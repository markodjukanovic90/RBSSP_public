"""
bbc_improve.py

Branch-and-Benders-Cut (BBC) za Minimum Red-Blue Separating Set Problem (MRBSSP).

Razlika u odnosu na ilp_improve.py (klasicni Set Cover ILP):
    - Master problem NE sadrzi sve covering constraint-e Sum_{v in D(r,b)} x_v >= 1
      unapred. Sadrzi samo objective (+ essential-node fixing + UB constraint iz
      GRASP warm starta).
    - Feasibility oracle (logicka provera, ne LP dualnost) je ugnjezdena u
      branch-and-bound preko CPLEX LazyConstraintCallback-a (docplex
      ConstraintCallbackMixin + cplex.callbacks.LazyConstraintCallback).
    - Svaki put kad CPLEX nadje novi integer incumbent, callback provjerava
      koji D(r,b) parovi (samo oni koji procure kroz kernelizaciju, tj. nisu
      vec pokriveni essential cvorom) nisu pokriveni izabranim separatorom.
      Iz skupa takvih parova bira se samo ogranicen broj (podrazumijevano 5)
      i dodaje se nekoliko lazy cut-ova po callback pozivu. Time prvi poziv
      callback-a ne ubacuje automatski cijelu preostalu ILP formulaciju.
      Samim dodavanjem cut-a koji kandidat narusava, CPLEX
      automatski tretira incumbent kao nevazeci (nema posebne reject()
      metode u CPLEX legacy LazyConstraintCallback API-ju).
    - Search tree se NE restartuje (za razliku od iterativnog LBBD) -- to je
      kljucna prednost BBC pristupa.

Kernelizacija (essential nodes, row/column dominance) i GRASP warm start su
preuzeti IDENTICNO iz ilp_improve.py, jer su nezavisni od toga da li se
covering constraint-i generisu unapred ili lenjo.

Pokretanje (na masini sa CPLEX + docplex instaliranim):
    python bbc_improve.py -i instance.csv -o rezultat.txt -t 60

Eksperiment sa zabranjenim cvorom/cvorovima (x_v = 0):
    python bbc_improve.py -i instance.csv -o rezultat_no_17.txt -t 3600 --forbid_nodes 17
    python bbc_improve.py -i instance.csv -o rezultat_no_17_42.txt -t 3600 --forbid_nodes 17 42
"""

import argparse
from time import time
from itertools import product
from os import remove

import networkx as nx

# ---------------------------------------------------------------------------
# CPLEX / docplex import (path isto kao u ilp_improve.py -- prilagodi po potrebi)
# ---------------------------------------------------------------------------
import os
import sys
import signal
import threading

# CPLEX Python API path. The environment variable takes precedence; the
# cluster path below is the user's current CPLEX installation.
_cplex_python_path = os.environ.get(
    "CPLEX_PYTHON_PATH",
    "/ceph/hpc/home/djukanovicm/VGLCS/cplex/python/3.10/x86-64_linux",
)
if _cplex_python_path and _cplex_python_path not in sys.path:
    sys.path.insert(0, _cplex_python_path)

import cplex
from docplex.mp.model import Model
from docplex.mp.callbacks.cb_mixin import ConstraintCallbackMixin
from cplex.callbacks import LazyConstraintCallback, MIPInfoCallback

# orloge is only needed by the optional log parser; BBC itself does not need it.
try:
    import orloge
except ImportError:
    orloge = None

from math import ceil
from grasp import grasp_mrbssp


# ---------------------------------------------------------------------------
# Graph & kernelization utilities -- IDENTICNO kao u ilp_improve.py
# ---------------------------------------------------------------------------

def closed_neighborhood(G, v):
    return set(G.neighbors(v)).union({v})


def compute_separation_sets(G, red_nodes, blue_nodes, selectable_nodes=None):
    """
    Za svaki red-blue par (r, b) racuna D(r,b) = skup SELEKTABILNIH cvorova v
    takvih da zatvoreno susedstvo N[v] sadrzi tacno jedan od {r, b}.

    ``selectable_nodes`` omogucava eksperiment sa zabranjenim cvorovima. Vazno:
    zabranjeni cvorovi se NE brisu iz grafa G, jer bi to promijenilo zatvorena
    susjedstva ostalih cvorova. Oni se samo iskljucuju kao moguci selektovani
    separatori, sto je ekvivalentno uslovu x_v = 0.
    """
    pair_sets = {}
    candidates = set(G.nodes()) if selectable_nodes is None else set(selectable_nodes)

    for r, b in product(red_nodes, blue_nodes):

        D_rb = set()

        for v in candidates:
            if len({r, b} & closed_neighborhood(G, v)) == 1:
                D_rb.add(v)

        pair_sets[(r, b)] = D_rb

    return pair_sets


def find_essential_nodes(pair_sets):
    """
    Ako je |D(r,b)| == 1, taj jedini cvor MORA biti u svakom validnom
    separatoru (jer je on jedini koji moze pokriti taj par).
    """
    essential = set()

    for D_rb in pair_sets.values():
        if len(D_rb) == 1:
            essential.update(D_rb)

    return essential


def dominated_rows(pair_sets):
    """
    Row dominance: ako D_i (strogo) podskup D_j, red j je dominiran i
    moze se ukloniti (njegovo pokrivanje je implicirano pokrivanjem D_i).
    """
    pairs = list(pair_sets.keys())
    to_remove = set()

    for i in range(len(pairs)):
        p1 = pairs[i]

        if p1 in to_remove:
            continue

        D1 = pair_sets[p1]

        for j in range(len(pairs)):

            if i == j:
                continue

            p2 = pairs[j]

            if p2 in to_remove:
                continue

            D2 = pair_sets[p2]

            if D1 < D2:
                to_remove.add(p2)

    reduced_pair_sets = {
        p: D
        for p, D in pair_sets.items()
        if p not in to_remove
    }

    return reduced_pair_sets, to_remove


def dominated_columns(pair_sets, nodes):
    """
    Column dominance: C(v) = skup parova koje v pokriva.
    Ako C(u) (strogo) podskup C(v), u je dominiran i moze se ukloniti
    (jer v pokriva sve sto i u, i moguce jos vise).
    """
    coverage = {v: set() for v in nodes}

    for pair, D in pair_sets.items():
        # Only surviving nodes are valid columns at this stage of kernelization.
        # Intersecting with ``nodes`` prevents stale vertices that were removed
        # in a previous dominance iteration from indexing ``coverage``.
        for v in D & nodes:
            coverage[v].add(pair)

    dominated = set()
    node_list = list(nodes)

    for i in range(len(node_list)):

        u = node_list[i]

        if u in dominated:
            continue

        Cu = coverage[u]

        for j in range(len(node_list)):

            if i == j:
                continue

            v = node_list[j]

            if v in dominated:
                continue

            Cv = coverage[v]

            if Cu < Cv:
                dominated.add(u)
                break

    return dominated


def kernelize(pair_sets, essential, all_nodes):
    """
    Iterativno primenjuje row i column dominance dok se nista vise ne moze
    ukloniti. Vraca redukovan pair_sets i redukovan skup cvorova.
    """
    changed = True

    removed_rows = 0
    removed_cols = 0

    while changed:

        changed = False

        new_pair_sets, removed = dominated_rows(pair_sets)

        if removed:
            removed_rows += len(removed)
            pair_sets = new_pair_sets
            changed = True

        dominated = dominated_columns(pair_sets, all_nodes - essential)

        if dominated:
            removed_cols += len(dominated)
            all_nodes -= dominated

            for pair in pair_sets:
                # Keep every separation set synchronized with the current
                # surviving-node set. This is stronger than subtracting only
                # the vertices removed in the present iteration and prevents
                # stale vertices from leaking into later dominance passes.
                pair_sets[pair] &= all_nodes

            changed = True

    print(f"Row dominance removed {removed_rows} constraints")
    print(f"Column dominance removed {removed_cols} variables")

    return pair_sets, all_nodes


def read_graph(path):
    """
    Cita graf iz CSV (DIMACS-like) fajla:
        n,<id>,<color>      -- color in {"red", "blue", ...}
        e,<u>,<v>
    """
    import csv

    G = nx.Graph()
    with open(path) as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == "n":
                G.add_node(int(row[1]), color=row[2])
            elif row[0] == "e":
                G.add_edge(int(row[1]), int(row[2]))

    red_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "red"}
    blue_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "blue"}
    return G, red_nodes, blue_nodes


def sanity_check(G, R, B, selectable_nodes=None):
    """
    Provjerava da je instanca rjesiva koristeci samo dozvoljene selektabilne
    cvorove. Ako su neki cvorovi zabranjeni, oni ostaju u grafu (i time u
    susjedstvima), ali se ne smiju koristiti kao separator.
    """
    pairs = product(R, B)
    candidates = set(G.nodes()) if selectable_nodes is None else set(selectable_nodes)

    for r, b in pairs:
        if not any(
            (r in closed_neighborhood(G, v)) != (b in closed_neighborhood(G, v))
            for v in candidates
        ):
            print(f"Unseparable pair under current selectable-node restriction: ({r}, {b})")
            return False

    return True


def checkFeasibility(solution, G, R, B):
    """
    Provjerava da dato rjesenje (skup cvorova) razdvaja SVE red-blue parove.
    """
    pairs = product(R, B)

    for r, b in pairs:
        uncovered = True
        for v in solution:
            if len(closed_neighborhood(G, v) & {r, b}) == 1:
                uncovered = False
                break
        if uncovered:
            return False

    return True


def retrieve_results_from_log(path_to_log):
    """
    Isto kao u ilp_improve.py -- parsira CPLEX log preko orloge.
    Napomena: kod lazy-constraint modela CPLEX log izgleda malo drugacije
    (broj generisanih "user cuts" se ne pojavljuje na isti nacin kao kod
    klasicnog B&B), ali best_bound / best_solution / status polja ostaju
    pouzdana.
    """
    if orloge is None:
        raise RuntimeError(
            "orloge is not installed. Install it only if retrieve_results_from_log() is used."
        )

    log_data = orloge.get_info_solver(path_to_log, "CPLEX")

    print("\n===== LOG DATA =====")
    print(log_data)
    print("====================\n")

    best_bound = log_data.get("best_bound")
    best_solution = log_data.get("best_solution")

    if best_bound is None:
        best_bound = 0

    if best_solution is None:
        best_solution = 0

    status = int(log_data.get("status") == "MIP - Integer optimal")

    remove(path_to_log)

    return best_solution, ceil(best_bound), status


# ---------------------------------------------------------------------------
# Crash-/SLURM-safe progress checkpointing
# ---------------------------------------------------------------------------

def _safe_float(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return default


class BBCProgressState:
    """Shared progress state written atomically to the requested output file.

    The output file is intentionally kept in the same simple ``Key: value``
    format as the final result file.  If SLURM kills the process, the latest
    checkpoint therefore remains parseable and still contains the best known
    valid UB/LB observed before termination.
    """

    def __init__(
        self,
        output_file,
        start_wall,
        initial_ub=None,
        initial_lb=None,
        initial_solution=None,
        forbidden_nodes=None,
        checkpoint_calls=50,
        checkpoint_nodes=1000,
        checkpoint_seconds=60.0,
    ):
        self.output_file = output_file
        self.start_wall = start_wall
        self.current_ub = _safe_float(initial_ub)
        self.current_lb = _safe_float(initial_lb)
        self.current_solution = (
            sorted(initial_solution) if initial_solution is not None else None
        )
        self.forbidden_nodes = sorted(forbidden_nodes or [])
        self.nb_cuts = 0
        self.nb_calls = 0
        self.max_violated_in_call = 0
        self.nodes = 0
        self.nodes_left = None
        self.last_checkpoint_call = 0
        self.last_checkpoint_nodes = 0
        self.last_checkpoint_wall = start_wall
        self.checkpoint_calls = max(0, int(checkpoint_calls))
        self.checkpoint_nodes = max(0, int(checkpoint_nodes))
        self.checkpoint_seconds = max(0.0, float(checkpoint_seconds))
        self.lock = threading.Lock()

    def _elapsed(self):
        return max(0.0, time() - self.start_wall)

    def update_bounds(self, lb=None, ub=None):
        # Minimization: LB should increase, UB should decrease.
        lb = _safe_float(lb)
        ub = _safe_float(ub)
        if lb is not None:
            if self.current_lb is None or lb > self.current_lb:
                self.current_lb = lb
        if ub is not None:
            if self.current_ub is None or ub < self.current_ub:
                self.current_ub = ub

    def update_incumbent(self, ub, solution):
        """Update the best known feasible incumbent and its vertex set atomically.

        The UB and stored solution must always correspond to the same feasible
        MRBSSP solution. This prevents checkpoints from containing a finite UB
        together with an empty or stale solution vector.
        """
        ub = _safe_float(ub)
        if ub is None or solution is None:
            return

        solution = sorted(solution)

        if self.current_ub is None or ub < self.current_ub - 1e-9:
            self.current_ub = ub
            self.current_solution = solution
        elif abs(ub - self.current_ub) <= 1e-9 and self.current_solution is None:
            self.current_solution = solution

    def checkpoint_due(self, calls=None, nodes=None):
        now = time()
        due_time = (
            self.checkpoint_seconds > 0
            and now - self.last_checkpoint_wall >= self.checkpoint_seconds
        )
        due_calls = (
            calls is not None
            and self.checkpoint_calls > 0
            and calls - self.last_checkpoint_call >= self.checkpoint_calls
        )
        due_nodes = (
            nodes is not None
            and self.checkpoint_nodes > 0
            and nodes - self.last_checkpoint_nodes >= self.checkpoint_nodes
        )
        return due_time or due_calls or due_nodes

    def write_checkpoint(self, reason="periodic", status="RUNNING", force=False):
        """Atomically replace the output file with the latest known progress."""
        with self.lock:
            if not force and not self.checkpoint_due(self.nb_calls, self.nodes):
                return

            directory = os.path.dirname(os.path.abspath(self.output_file)) or "."
            os.makedirs(directory, exist_ok=True)
            tmp_path = self.output_file + ".tmp"

            lb = self.current_lb if self.current_lb is not None else -1
            ub = self.current_ub if self.current_ub is not None else -1
            elapsed = self._elapsed()
            if ub is not None and ub != -1 and lb is not None and lb != -1 and ub > 0:
                gap = 100.0 * max(0.0, ub - lb) / ub
            else:
                gap = -1

            solution = (
                sorted(self.current_solution)
                if self.current_solution is not None
                else []
            )
            feasible_text = "True" if self.current_solution is not None else "unknown"

            with open(tmp_path, "w") as f:
                f.write(f"Solution: {solution}\n")
                f.write(f"LB: {lb}\n")
                f.write(f"UB: {ub}\n")
                f.write(f"Gap: {gap}\n")
                f.write(f"Status: {status}\n")
                f.write(f"Feasible: {feasible_text}\n")
                f.write(f"Time: {elapsed:.2f}\n")
                f.write(f"LazyCuts: {self.nb_cuts}\n")
                f.write(f"CallbackCalls: {self.nb_calls}\n")
                f.write(f"NodesProcessed: {self.nodes}\n")
                if self.nodes_left is not None:
                    f.write(f"NodesLeft: {self.nodes_left}\n")
                f.write(f"MaxViolatedPairsInOneCall: {self.max_violated_in_call}\n")
                f.write(f"ForbiddenNodes: {self.forbidden_nodes}\n")
                f.write("Checkpoint: True\n")
                f.write(f"CheckpointReason: {reason}\n")
                f.flush()
                os.fsync(f.fileno())

            os.replace(tmp_path, self.output_file)
            self.last_checkpoint_wall = time()
            self.last_checkpoint_call = self.nb_calls
            self.last_checkpoint_nodes = self.nodes


class BBCInfoCallback(MIPInfoCallback):
    """Periodically records global B&B information while CPLEX is running."""

    def __init__(self, env):
        MIPInfoCallback.__init__(self, env)
        self.progress_state = None

    def setup(self, progress_state):
        self.progress_state = progress_state

    def __call__(self):
        st = self.progress_state
        if st is None:
            return

        try:
            st.nodes = int(self.get_num_nodes())
        except Exception:
            pass
        try:
            st.nodes_left = int(self.get_num_remaining_nodes())
        except Exception:
            pass

        # Global best bound is a valid LB for the minimization master.
        try:
            st.update_bounds(lb=self.get_best_objective_value())
        except Exception:
            pass

        # Do not update the UB here. The MIP info callback exposes the
        # incumbent objective, but not the synchronized incumbent vertex set
        # needed by crash-safe checkpointing. The lazy callback updates UB and
        # solution together only after the candidate has passed the complete
        # MRBSSP feasibility check.

        if st.checkpoint_due(st.nb_calls, st.nodes):
            st.write_checkpoint(reason="MIP-progress", status="RUNNING", force=True)


# ---------------------------------------------------------------------------
# BBC -- Lazy Constraint Callback (logicki feasibility oracle)
# ---------------------------------------------------------------------------

class MRBSSPLazyCallback(ConstraintCallbackMixin, LazyConstraintCallback):
    """
    Implementira "Logic-Based" feasibility oracle iz PDF-a (sekcija 4/6):

        - Master problem (samo objective + essential/UB constraint-i) je
          vec rijesen do integer incumbent-a od strane CPLEX-a.
        - Ovaj callback se poziva SVAKI PUT kad B&B nadje novi integer
          incumbent (cak i pre nego sto je formalno prihvacen).
        - Mi pretvaramo trenutne vrednosti x_v u kandidatski separator S,
          pa logicki provjerimo da li S pokriva sve preostale (D_rb) parove
          (essential cvorovi su vec zagarantovano u svakom rjesenju zbog
          fixing constraint-a x_v == 1 u masteru, pa ih ne moramo provjeravati).
        - Ako postoje nepokriveni parovi, biramo najvise k njih i dodajemo
          odgovarajuce lazy cut-ove Sum_{v in D_remaining} x_v >= 1.
          Preostali redovi se ostavljaju za naredne callback pozive. Sam
          dodati cut, kojeg trenutni kandidat narusava, je dovoljan da CPLEX
          tretira incumbent kao nevazeci -- ne postoji posebna reject()
          metoda u CPLEX legacy LazyConstraintCallback API-ju. Search tree
          NIJE restartovan -- CPLEX jednostavno nastavlja B&B sa novim
          cut-ovima u poolu.
        - Ako su svi parovi pokriveni, incumbent je validan i prihvata se
          (ne radimo nista -- CPLEX ga automatski prihvata).
    """

    def __init__(self, env):
        LazyConstraintCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
        self.pair_sets = None          # dict[(r,b)] -> set(vertex) -- SAMO preostali parovi (post-kernelizacija)
        self.x = None                  # dict[vertex] -> docplex Var (samo "remaining_nodes")
        self.var_index = None          # dict[vertex] -> CPLEX-ov integer indeks varijable x_v
        self.remaining_nodes = None    # set preostalih (ne-essential, ne-dominiranih) cvorova
        self.essential = set()           # essential vertices present in every full solution
        self.nb_cuts = 0                # brojac generisanih lazy cut-ova (za statistiku)
        self.nb_calls = 0               # brojac poziva callback-a (za statistiku)
        self.max_cuts_per_call = 5      # najvise lazy cut-ova po callback pozivu
        self.cut_strategy = "smallest_set"
        self.added_pairs = set()        # sigurnost od ponovnog dodavanja istog reda
        self.max_violated_in_call = 0   # dijagnostika
        self.progress_state = None       # shared checkpoint state

    def setup(
        self,
        pair_sets,
        x,
        remaining_nodes,
        essential,
        max_cuts_per_call=5,
        cut_strategy="smallest_set",
        progress_state=None,
    ):
        self.pair_sets = pair_sets
        self.x = x
        self.remaining_nodes = remaining_nodes
        self.essential = set(essential)
        # docplex Var.index daje indeks kolone u CPLEX modelu -- to je ono
        # sto native get_values()/add() API ocekuje (ne docplex objekte).
        self.var_index = {v: x[v].index for v in remaining_nodes}
        self.max_cuts_per_call = max(1, int(max_cuts_per_call))
        self.cut_strategy = cut_strategy
        self.added_pairs = set()
        self.progress_state = progress_state

    def _select_violated_pairs(self, violated_pairs):
        """Bira najvise k trenutno narusenih redova.

        Kod integer kandidata svaki naruseni red ima istu numericku
        violaciju (lijeva strana je 0), pa je podrazumijevana strategija
        ``smallest_set``: prvo se dodaju kraci, tipicno restriktivniji redovi.
        """
        if self.cut_strategy == "smallest_set":
            ordered = sorted(
                violated_pairs,
                key=lambda item: (len(item[1]), item[0]),
            )
        elif self.cut_strategy == "largest_set":
            ordered = sorted(
                violated_pairs,
                key=lambda item: (-len(item[1]), item[0]),
            )
        elif self.cut_strategy == "lexicographic":
            ordered = sorted(violated_pairs, key=lambda item: item[0])
        else:
            raise ValueError(
                f"Unknown cut strategy: {self.cut_strategy}"
            )

        return ordered[: self.max_cuts_per_call]

    def __call__(self):
        """
        Pozvan od CPLEX-a za SVAKI integer-feasible candidate (incumbent)
        koji B&B pronadje, PRIJE nego sto je formalno prihvacen.

        Logicki feasibility oracle (PDF sekcija 4/6):
            1. Procitaj trenutne vrijednosti x_v (native CPLEX get_values).
            2. Rekonstruisi kandidatski separator S.
            3. Za svaki preostali (r,b) par provjeri da li S sece D(r,b).
            4. Izaberi najvise k nepokrivenih parova i dodaj njihove lazy
               cut-ove. Najmanje jedan rez je dovoljan da trenutni kandidat
               bude odbacen; ostale violacije se otkrivaju kasnije.
               Samo dodavanje cut-a koji kandidat narusava je dovoljno da
               CPLEX tretira incumbent kao nevazeci -- legacy
               LazyConstraintCallback API nema posebnu reject() metodu.
            Ako su svi parovi pokriveni, ne radi se nista -- CPLEX ga
            automatski prihvata kao validan incumbent.
        """
        self.nb_calls += 1

        indices = list(self.var_index.values())
        vertices = list(self.var_index.keys())

        # Native CPLEX legacy-callback metoda: vraca listu vrijednosti u
        # istom redosledu kao predati indeksi varijabli.
        values = self.get_values(indices)

        candidate_S = {
            v for v, val in zip(vertices, values) if val > 0.5
        }

        violated_pairs = [
            (pair, D_rb)
            for pair, D_rb in self.pair_sets.items()
            if pair not in self.added_pairs and not (D_rb & candidate_S)
        ]
        self.max_violated_in_call = max(
            self.max_violated_in_call, len(violated_pairs)
        )

        if self.progress_state is not None:
            st = self.progress_state
            st.nb_calls = self.nb_calls
            st.max_violated_in_call = self.max_violated_in_call
            try:
                st.update_bounds(lb=self.get_best_objective_value())
            except Exception:
                pass
            # Do NOT use the current candidate as UB until it has passed the
            # full lazy-feasibility check below.

        if violated_pairs:
            selected_pairs = self._select_violated_pairs(violated_pairs)

            # Dodavanje makar jednog narusenog lazy reda odbacuje trenutni
            # kandidat. Preostali naruseni redovi ce biti otkriveni pri
            # narednim integer kandidatima, bez restartovanja B&B stabla.
            for pair, D_rb in selected_pairs:
                self._add_lazy_cut(D_rb)
                self.added_pairs.add(pair)
                self.nb_cuts += 1

            if self.progress_state is not None:
                st = self.progress_state
                st.nb_cuts = self.nb_cuts
                st.nb_calls = self.nb_calls
                if st.checkpoint_due(st.nb_calls, st.nodes):
                    st.write_checkpoint(reason="lazy-callback", status="RUNNING", force=True)

            # Napomena: CPLEX-ov legacy LazyConstraintCallback NEMA metodu
            # reject()/accept() -- to nije deo API-ja (za razliku od neke
            # druge dokumentacije/jezika gdje takva metoda postoji). Sam
            # poziv self.add(...) sa cut-om koji trenutni kandidat narusava
            # je dovoljan: CPLEX automatski tretira kandidat kao nevazeci
            # (jer ne zadovoljava novododati cut) i nastavlja B&B sa njim
            # u poolu lazy constraint-a -- nista vise ne treba pozvati.

        # Ako violated_pairs je prazan -- kandidat je validan za ORIGINALNI
        # MRBSSP. Tek tada je njegova ciljna vrijednost siguran UB.
        if not violated_pairs and self.progress_state is not None:
            st = self.progress_state

            # candidate_S contains only surviving non-essential decision
            # vertices. Restore essential vertices before checkpointing.
            full_candidate = set(self.essential) | candidate_S

            try:
                candidate_obj = self.get_objective_value()
            except Exception:
                candidate_obj = len(full_candidate)

            # Store UB and its corresponding feasible solution atomically.
            st.update_incumbent(candidate_obj, full_candidate)
            st.nb_calls = self.nb_calls
            st.nb_cuts = self.nb_cuts
            if st.checkpoint_due(st.nb_calls, st.nodes):
                st.write_checkpoint(reason="valid-incumbent", status="RUNNING", force=True)

        # Ako violated_pairs je prazan -- incumbent je feasible (pokriva
        # sve preostale parove) i CPLEX ga automatski prihvata.

    def _add_lazy_cut(self, D_rb):
        """
        Dodaje lazy constraint Sum_{v in D_rb} x_v >= 1 preko native CPLEX
        legacy-callback API-ja (cplex.SparsePair), sto je standardni nacin
        dodavanja lazy cut-ova unutar LazyConstraintCallback-a.
        """
        indices = [self.var_index[v] for v in D_rb]
        coeffs = [1.0] * len(indices)

        self.add(
            constraint=cplex.SparsePair(ind=indices, val=coeffs),
            sense="G",   # >=
            rhs=1.0
        )


# ---------------------------------------------------------------------------
# BBC model generation
# ---------------------------------------------------------------------------

def generate_bbc_model(
    G,
    pair_sets,
    essential,
    warm_start,
    all_nodes,
    ub=None,
):
    """
    Build the BBC master problem.

    IMPORTANT: ``ub`` and ``warm_start`` are optional. They must only be
    supplied when the heuristic solution has been verified feasible for the
    ORIGINAL MRBSSP instance. An infeasible heuristic candidate is not a valid
    upper bound and must never restrict the BBC master.

    Covering constraints are still omitted here and generated lazily by the
    callback. Essential nodes are accounted for as a constant contribution to
    the objective.
    """
    mdl = Model(name="MRBSSP_BBC")

    remaining_nodes = all_nodes - essential

    x = {
        v: mdl.binary_var(name=f"x_{v}")
        for v in remaining_nodes
    }

    total_objective = len(essential) + mdl.sum(x[v] for v in remaining_nodes)
    mdl.minimize(total_objective)

    # A GRASP value is a valid upper bound ONLY when GRASP was verified
    # feasible for the complete original instance.
    if ub is not None:
        mdl.add_constraint(total_objective <= int(ub), ctname="grasp_upper_bound")

    # Add a MIP start only when a representable, verified-feasible warm start
    # is provided. (A valid UB may still be used even when the original warm
    # start contains columns removed by kernelization and therefore cannot be
    # represented directly in the reduced master.)
    if warm_start is not None:
        warm_start = set(warm_start)
        mip_start = mdl.new_solution()
        for v in remaining_nodes:
            mip_start.add_var_value(x[v], 1 if v in warm_start else 0)
        mdl.add_mip_start(mip_start)

    return mdl, x, remaining_nodes

def solve_bbc(
    mdl,
    x,
    remaining_nodes,
    pair_sets,
    essential,
    time_limit,
    threads=1,
    max_cuts_per_call=5,
    cut_strategy="smallest_set",
    progress_state=None,
):
    """
    Registruje LazyConstraintCallback i pokrece B&B sa lenjim cut-ovima.

    pair_sets: vec filtriran na D_remaining (presjek sa remaining_nodes) za
               sve parove koji NISU vec pokriveni essential cvorom -- ovo se
               racuna jednom unapred (van callback-a) radi efikasnosti, jer
               se D_rb skupovi ne mijenjaju tokom B&B, samo se njihova
               "pokrivenost" provjerava na svakom incumbentu.
    """
    mdl.parameters.timelimit = time_limit
    mdl.parameters.threads = threads

    # KLJUCNO za korektnost lazy constraint-a u CPLEX-u: CPX_PARAM_REDUCE
    # (preprocessing.reduce) MORA biti 0 (nikakve redukcije) ili 1 (samo
    # primalne redukcije) kad postoji lazy-constraint pool/callback --
    # dualne redukcije pretpostavljaju da je model kompletan, sto ovde
    # nije slucaj (vidi CPLEX dokumentaciju: "Adding User Cuts and Lazy
    # Constraints"). Vrednost 1 = CPX_PREREDUCE_PRIMALONLY.
    mdl.parameters.preprocessing.reduce = 1

    # Napomena: ne postoji posebаn parametar koji "ukljucuje" lazy-constraint
    # rezim -- CPLEX ga automatski aktivira kada se registruje
    # LazyConstraintCallback (sto radimo redom ispod).
    cb = mdl.register_callback(MRBSSPLazyCallback)
    cb.setup(
        pair_sets,
        x,
        remaining_nodes,
        essential,
        max_cuts_per_call=max_cuts_per_call,
        cut_strategy=cut_strategy,
        progress_state=progress_state,
    )

    # MIPInfoCallback supplies global LB / accepted incumbent UB / node counts
    # even when no new lazy incumbent arrives for a long time.
    info_cb = mdl.register_callback(BBCInfoCallback)
    info_cb.setup(progress_state)

    sol = mdl.solve(log_output=True)

    return sol, cb


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Branch-and-Benders-Cut (BBC) Red-Blue Separator sa CPLEX lazy constraints"
    )
    parser.add_argument("-i", "--i", required=True, help="Input CSV graph file (DIMACS-like format)")
    parser.add_argument("-o", "--o", required=True, help="Output file path for results")
    parser.add_argument("-t", "--time_limit", type=float, default=60.0, help="CPLEX time limit in seconds")
    parser.add_argument("-gk", "--grasp_k", type=int, default=5, help="k parametar za GRASP greedy")
    parser.add_argument("-gi", "--grasp_iter", type=int, default=10, help="Broj GRASP iteracija")
    parser.add_argument(
        "-f",
        "--forbid_nodes",
        "--forbid_node",
        nargs="+",
        type=int,
        default=[],
        metavar="NODE_ID",
        help=(
            "Jedan ili vise cvorova koji ne smiju biti u rjesenju. "
            "Primjer: --forbid_nodes 17 ili --forbid_nodes 17 42 103. "
            "Zabranjeni cvorovi se iskljucuju iz svih D(r,b) skupova, sto je "
            "ekvivalentno x_v=0, ali se NE brisu iz grafa."
        ),
    )
    parser.add_argument(
        "-kc",
        "--max_cuts_per_callback",
        type=int,
        default=5,
        help="Najveci broj narusenih lazy cut-ova dodatih u jednom callback pozivu (default: 5)",
    )
    parser.add_argument(
        "--cut_strategy",
        choices=("smallest_set", "largest_set", "lexicographic"),
        default="smallest_set",
        help="Strategija izbora narusenih parova kada ih ima vise od dozvoljenog broja",
    )
    parser.add_argument(
        "--checkpoint_calls",
        type=int,
        default=50,
        help="Write BBC progress every N lazy-callback calls; 0 disables (default: 50)",
    )
    parser.add_argument(
        "--checkpoint_nodes",
        type=int,
        default=1000,
        help="Write BBC progress every N processed B&B nodes; 0 disables (default: 1000)",
    )
    parser.add_argument(
        "--checkpoint_seconds",
        type=float,
        default=60.0,
        help="Write BBC progress at least every N seconds when a callback is reached; 0 disables (default: 60)",
    )
    args = parser.parse_args()

    if args.max_cuts_per_callback < 1:
        parser.error("--max_cuts_per_callback mora biti najmanje 1")
    if args.checkpoint_calls < 0 or args.checkpoint_nodes < 0 or args.checkpoint_seconds < 0:
        parser.error("Checkpoint intervals must be non-negative")

    input_file = args.i
    output_file = args.o
    time_limit = float(args.time_limit)

    G, R, B = read_graph(input_file)

    # -----------------------------------------------------------------------
    # Optional forbidden-node experiment.
    # -----------------------------------------------------------------------
    # A forbidden node remains in G so that the original closed neighborhoods
    # are unchanged, but it is removed from the selectable node universe. This
    # is exactly the intended x_v = 0 experiment and also means that the node is
    # absent from every D(r,b) set used by preprocessing and BBC.
    forbidden_nodes = set(args.forbid_nodes or [])
    unknown_forbidden = forbidden_nodes - set(G.nodes())
    if unknown_forbidden:
        parser.error(
            "Unknown forbidden node id(s): "
            + ", ".join(str(v) for v in sorted(unknown_forbidden))
        )

    selectable_nodes = set(G.nodes()) - forbidden_nodes
    if forbidden_nodes:
        print("Forbidden nodes (x_v = 0):", sorted(forbidden_nodes))
        print("Selectable nodes after prohibition:", len(selectable_nodes))

    # Check feasibility under the forbidden-node restriction BEFORE GRASP/BBC.
    sanity = sanity_check(G, R, B, selectable_nodes=selectable_nodes)

    heuristic_solution = None
    heuristic_value = None
    grasp_feasible = False

    if sanity:
        print(f"Running GRASP warm start with {args.grasp_iter} iterations and k={args.grasp_k}")
        try:
            
            from LS import build_pair_cover_map
            pair_sets_grasp =  build_pair_cover_map(G, R, B)
            candidate = grasp_mrbssp(G=G, R=R, B=B, k=int(args.grasp_k), max_iter=int(args.grasp_iter), seed=None, pair_cover_map=pair_sets_grasp) #, time_limit

        except Exception as exc:
            print(f"WARNING: GRASP failed ({type(exc).__name__}: {exc}).")
            print("BBC will continue without a warm start or GRASP upper bound.")
            candidate = None

        if candidate is not None:
            candidate = set(candidate)
            uses_forbidden = bool(candidate & forbidden_nodes)
            grasp_feasible = checkFeasibility(candidate, G, R, B) and not uses_forbidden
            print("GRASP feasible under forbidden-node restriction =", grasp_feasible)

            if uses_forbidden:
                print(
                    "WARNING: GRASP warm start uses forbidden node(s): "
                    f"{sorted(candidate & forbidden_nodes)}. "
                    "Disabling both the MIP start and the GRASP upper-bound constraint."
                )

            if grasp_feasible:
                heuristic_solution = list(candidate)
                heuristic_value = len(heuristic_solution)
                print(f"Valid GRASP upper bound = {heuristic_value}")
            else:
                if not uses_forbidden:
                    print(
                        "WARNING: GRASP returned an infeasible candidate. "
                        "Disabling both the MIP start and the GRASP upper-bound constraint."
                    )
        else:
            print("GRASP feasible = False")
    else:
        print("Instance failed sanity_check(); skipping GRASP and BBC solve.")

    total_time = 0.0
    nb_cuts = 0
    nb_calls = 0
    max_violated_in_call = 0

    if sanity:

        # Build D(r,b) only from selectable vertices. Forbidden vertices are
        # therefore not part of any separation set and can never become essential.
        pair_sets = compute_separation_sets(
            G, R, B, selectable_nodes=selectable_nodes
        )

        # Defensive check: banning a word can make the constrained variant
        # infeasible if some red-blue pair loses its last possible separator.
        empty_pairs = [pair for pair, D_rb in pair_sets.items() if not D_rb]
        if empty_pairs:
            raise RuntimeError(
                "Internal inconsistency: sanity_check passed but a constrained "
                f"D(r,b) set is empty; first pair: {empty_pairs[0]}"
            )

        essential = find_essential_nodes(pair_sets)
        print(f"Essential nodes under restriction: {len(essential)}")

        all_nodes = set(selectable_nodes)
        pair_sets, all_nodes = kernelize(pair_sets, essential, all_nodes)

        print("Remaining nodes after kernelization:", len(all_nodes))
        print("Remaining rows:", len(pair_sets))

        # Pripremi pair_sets ZA CALLBACK: ukloni essential cvorove iz D_rb i
        # iskljuci parove vec pokrivene essential cvorom -- identicna logika
        # kao filter u generate_docplex_model() iz ilp_improve.py, samo sto
        # se rezultat ne koristi za pravljenje constraint-a unapred, vec se
        # cuva za lenju provjeru u callback-u.
        remaining_nodes_for_cb = all_nodes - essential
        callback_pair_sets = {}

        for (r, b), D_rb in pair_sets.items():

            if D_rb & essential:
                continue

            D_remaining = D_rb & remaining_nodes_for_cb

            if D_remaining:
                callback_pair_sets[(r, b)] = D_remaining

        print("Pairs requiring lazy cuts (worst case):", len(callback_pair_sets))

        # A feasible GRASP solution always gives a valid UB for the original
        # problem. However, if kernelization removed a selected GRASP column,
        # that exact MIP start cannot be represented in the reduced master. In
        # that case keep the valid UB but omit the MIP start.
        warm_start_for_model = None
        if heuristic_solution is not None:
            selected_nonessential = set(heuristic_solution) - essential
            if selected_nonessential <= (all_nodes - essential):
                warm_start_for_model = heuristic_solution
            else:
                print(
                    "GRASP solution contains columns removed by kernelization; "
                    "keeping its valid UB but omitting the MIP start."
                )

        mdl, x_vars, remaining_nodes = generate_bbc_model(
            G,
            pair_sets,
            essential,
            warm_start_for_model,
            all_nodes,
            ub=heuristic_value,
        )

        start = time()
        progress_state = BBCProgressState(
            output_file=output_file,
            start_wall=start,
            initial_ub=heuristic_value,
            initial_lb=len(essential),
            initial_solution=heuristic_solution,
            forbidden_nodes=forbidden_nodes,
            checkpoint_calls=args.checkpoint_calls,
            checkpoint_nodes=args.checkpoint_nodes,
            checkpoint_seconds=args.checkpoint_seconds,
        )
        progress_state.write_checkpoint(
            reason="before-BBC-solve", status="RUNNING", force=True
        )

        # Best effort on graceful SLURM termination (SIGTERM). SIGKILL cannot
        # be caught, which is why periodic atomic checkpoints are still needed.
        def _termination_handler(signum, frame):
            try:
                progress_state.write_checkpoint(
                    reason=f"signal-{signum}", status="INTERRUPTED", force=True
                )
            finally:
                raise SystemExit(128 + signum)

        old_sigterm = signal.signal(signal.SIGTERM, _termination_handler)
        old_sigint = signal.signal(signal.SIGINT, _termination_handler)

        try:
            sol, cb = solve_bbc(
                mdl,
                x_vars,
                remaining_nodes,
                callback_pair_sets,
                essential,
                time_limit,
                max_cuts_per_call=args.max_cuts_per_callback,
                cut_strategy=args.cut_strategy,
                progress_state=progress_state,
            )
        finally:
            signal.signal(signal.SIGTERM, old_sigterm)
            signal.signal(signal.SIGINT, old_sigint)

        total_time = round(time() - start, 2)

        nb_cuts = cb.nb_cuts
        nb_calls = cb.nb_calls
        max_violated_in_call = cb.max_violated_in_call
        progress_state.nb_cuts = nb_cuts
        progress_state.nb_calls = nb_calls
        progress_state.max_violated_in_call = max_violated_in_call

        if sol:

            bbc_solution = list(essential)

            for v in remaining_nodes:
                if x_vars[v].solution_value > 0.5:
                    bbc_solution.append(v)

            ub = mdl.objective_value

            try:
                lb = mdl.solve_details.best_bound
            except Exception:
                lb = ub

            status = mdl.solve_details.status

            feas = checkFeasibility(bbc_solution, G, R, B)

        else:

            # Even if CPLEX has no accepted incumbent, a previously verified
            # GRASP solution is still a valid MRBSSP upper bound.
            bbc_solution = list(heuristic_solution) if heuristic_solution is not None else []
            ub = heuristic_value if heuristic_value is not None else -1

            try:
                lb = mdl.solve_details.best_bound
            except Exception:
                lb = progress_state.current_lb if progress_state.current_lb is not None else -1

            status = "No CPLEX incumbent"
            feas = bool(heuristic_solution is not None)

    else:

        bbc_solution = []
        ub = -1
        lb = -1
        status = (
            "Infeasible under forbidden-node restriction (sanity check failed)"
            if forbidden_nodes
            else "Infeasible instance (sanity check failed)"
        )
        feas = False

    # Final result replaces the running checkpoint atomically as well.
    tmp_output = output_file + ".tmp"
    if ub is not None and ub != -1 and lb is not None and lb != -1 and ub > 0:
        final_gap = 100.0 * max(0.0, float(ub) - float(lb)) / float(ub)
    else:
        final_gap = -1
    with open(tmp_output, "w") as f:
        f.write(f"Solution: {sorted(bbc_solution)}\n")
        f.write(f"LB: {lb}\n")
        f.write(f"UB: {ub}\n")
        f.write(f"Gap: {final_gap}\n")
        f.write(f"Status: {status}\n")
        f.write(f"Feasible: {feas}\n")
        f.write(f"Time: {total_time}\n")

        #initial (GRASP) solution stats:
        f.write(f"InitialSolutionValue: {heuristic_value}\n")
        f.write(f"InitialSolutionNodes: { sorted(heuristic_solution) if heuristic_solution != None else ''}\n")

        f.write(f"LazyCuts: {nb_cuts}\n")
        f.write(f"CallbackCalls: {nb_calls}\n")
        f.write(f"MaxCutsPerCallback: {args.max_cuts_per_callback}\n")
        f.write(f"CutStrategy: {args.cut_strategy}\n")
        f.write(f"MaxViolatedPairsInOneCall: {max_violated_in_call}\n")
        f.write(f"ForbiddenNodes: {sorted(forbidden_nodes)}\n")
        f.write("Checkpoint: False\n")
        f.write("CheckpointReason: final\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp_output, output_file)


if __name__ == "__main__":
    main()

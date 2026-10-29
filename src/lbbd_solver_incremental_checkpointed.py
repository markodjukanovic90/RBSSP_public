#!/usr/bin/env python3
from docplex.mp.model import Model

import argparse
import csv
from pathlib import Path
from time import time
from itertools import product
from collections import defaultdict
from os import name as os_name, remove
import os
import signal
import math


import networkx as nx
import pulp as pl
#import orloge  # pip install orloge
from math import ceil

from grasp import grasp_mrbssp

# Detect platform and set CPLEX path
import sys
#sys.path.append("/home/marko/Desktop/cplex/cplex/python/3.10/x86-64_linux")
#

sys.path.insert(
0,
"/ceph/hpc/home/djukanovicm/VGLCS/cplex/python/3.10/x86-64_linux"
)


# ---------------------- Graph & Model Utilities ---------------------------- #

def read_graph(path):
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

def closed_neighborhood(G, v):
    return set(G.neighbors(v)).union({v})
    
    
def sanity_check(G, R, B):

    pairs=product(R, B)
    
    for r, b in pairs:
        if not any((r in closed_neighborhood(G, v)) != (b in closed_neighborhood(G, v)) for v in G.nodes()):
            print(f"Unseparable pair: ({r}, {b})")
            return False
     
    return True 
    
def checkFeasibility(solution, G, R, B):

    pairs=product(R, B)
    
    for r, b in pairs:
        uncovered=True
        for v in solution:
            if len(closed_neighborhood(G, v) & {r, b}) == 1:
                uncovered=False
        if uncovered:
            return False 
    
    return True


def _relative_gap(lb, ub):
    """Return 100*(UB-LB)/UB when both bounds are finite and UB > 0."""
    try:
        if lb is None or ub is None:
            return None
        lb = float(lb) 
        ub = float(ub)
        if not math.isfinite(lb) or not math.isfinite(ub) or ub <= 0:
            return None
        return max(0.0, 100.0 * (ub - lb) / ub)
    except (TypeError, ValueError):
        return None


def _write_lbbd_checkpoint(path, stats, reason="periodic", checkpoint=True, other_stats=None):
    """Atomically persist the latest LBBD state.

    The temporary-file + os.replace pattern prevents SLURM termination from
    leaving a partially written result file.  UB is always the best solution
    already verified feasible for the *complete* MRBSSP; a restricted-master
    incumbent is never reported as a full-problem UB until it passes the oracle.
    """
    if not path:
        return

    lb = ceil(stats.get("best_lower_bound")) if stats.get("best_lower_bound") is not None else None     
    ub = stats.get("best_feasible_value")
    gap = _relative_gap(lb, ub)
    tmp_path = f"{path}.tmp.{os.getpid()}"

    lines = [
        f"status: {stats.get('status', 'running')}\n",
        f"LB: {lb}\n",
        f"UB: {ub}\n",
        f"gap: {gap}\n",
        f"initial_value: {stats.get('initial_value')}\n",
        f"initial_solution: {stats.get('initial_solution')}\n",
        f"best_lower_bound: {lb}\n",
        f"best_feasible_value: {ub}\n",
        f"best_feasible_solution: {stats.get('best_feasible_solution')}\n",
        f"last_master_objective: {stats.get('last_master_objective')}\n",
        f"last_master_solution: {stats.get('last_master_solution')}\n",
        f"iterations: {stats.get('iterations', 0)}\n",
        f"cuts_added: {stats.get('cuts_added', 0)}\n",
        f"elapsed: {float(stats.get('elapsed', 0.0)):.6f}\n",
        f"checkpoint: {checkpoint}\n",
        f"checkpoint_reason: {reason}\n",
        f"kernel time: {other_stats['kernelization_time']:.6f}\n" if other_stats and 'kernelization_time' in other_stats else "",  
        f"grasp time: {other_stats['grasp_time']:.6f}\n" if other_stats and 'grasp_time' in other_stats else ""
    ]
     
    with open(tmp_path, "w") as out_file:
        out_file.writelines(lines)
        out_file.flush()
        os.fsync(out_file.fileno())
    os.replace(tmp_path, path)


################### LBBD ####################################################

def find_violated_pairs(solution, pair_sets, already_added=None):
    """Return violated pair constraints that have not yet been added."""
    selected = set(solution)
    already_added = already_added or set()
    return [
        pair for pair, separating_nodes in pair_sets.items()
        if pair not in already_added and not (selected & separating_nodes)
    ]


def select_violated_pairs(violated, pair_sets, max_cuts=1, strategy="smallest_set"):
    """Select at most ``max_cuts`` violated pair constraints.

    Strategies
    ----------
    smallest_set:
        Add constraints with the smallest |D_rb| first. These are typically
        the most restrictive rows and are a natural proxy for strong cuts.
    largest_set:
        Add constraints with the largest |D_rb| first.
    lexicographic:
        Deterministic pair ordering.
    all:
        Add every currently violated constraint.
    """
    if strategy == "all" or max_cuts is None or max_cuts <= 0:
        return list(violated)

    if strategy == "smallest_set":
        ordered = sorted(violated, key=lambda p: (len(pair_sets[p]), p))
    elif strategy == "largest_set":
        ordered = sorted(violated, key=lambda p: (-len(pair_sets[p]), p))
    elif strategy == "lexicographic":
        ordered = sorted(violated)
    else:
        raise ValueError(f"Unknown cut-selection strategy: {strategy}")

    return ordered[:max_cuts]

def build_master(nodes, essential=None, warm_start=None, ub=None):
    mdl = Model("MRBSSP_LBBD")
    x = {v: mdl.binary_var(name=f"x_{v}") for v in nodes}
    mdl.minimize(mdl.sum(x.values()))

    essential = set(essential or [])
    for v in essential:
        if v in x:
            mdl.add_constraint(x[v] == 1, ctname=f"essential_{v}")

    if ub is not None:
        mdl.add_constraint(mdl.sum(x.values()) <= ub, ctname="upper_bound")

    if warm_start is not None:
        warm_start = set(warm_start)
        ms = mdl.new_solution()
        for v in nodes:
            ms.add_var_value(x[v], 1 if v in warm_start or v in essential else 0)
        mdl.add_mip_start(ms)

    return mdl, x


def solve_lbbd(
    nodes,
    pair_sets,
    essential=None,
    warm_start=None,
    ub=None,
    time_limit=60,
    max_cuts_per_iteration=1,
    cut_strategy="smallest_set",
    log_output=True,
    checkpoint_path=None,
    checkpoint_every_iterations=1,
    checkpoint_seconds=60.0,
    other_stats=None     
):
    """Solve MRBSSP by incremental LBBD row generation.

    Only a selected subset of currently violated pair constraints is added
    after each master solve. Besides the final solution, the routine records:

    * the supplied initial (GRASP) solution and its objective value;
    * the best valid lower bound obtained over all restricted masters;
    * the best complete feasible solution known during the run.

    The optimum of every restricted master is a valid lower bound for the
    complete formulation because the restricted master contains only a subset
    of the original covering constraints.
    """
    start_time = time()
    essential = set(essential or [])
    initial_solution = sorted(set(warm_start or []))
    initial_value = len(initial_solution) if warm_start is not None else None

    mdl, x = build_master(nodes, essential, warm_start, ub)
    mdl.parameters.threads = 1

    added_pairs = set()
    iteration = 0
    best_lower_bound = float("-inf")
    last_master_solution = None
    last_master_objective = None

    # A feasible warm start is already a valid full-problem incumbent.
    best_feasible_solution = list(initial_solution) if warm_start is not None else None
    best_feasible_value = initial_value

    last_checkpoint_time = start_time
    last_checkpoint_iteration = 0
    previous_sigterm_handler = signal.getsignal(signal.SIGTERM)

    def make_stats(status, other_stats=None):
        return {  
            "status": status,
            "iterations": iteration,
            "cuts_added": len(added_pairs),
            "elapsed": time() - start_time,
            "initial_solution": initial_solution,
            "initial_value": initial_value,
            "best_lower_bound": None if best_lower_bound == float("-inf") else ceil(best_lower_bound),
            "best_feasible_solution": best_feasible_solution,
            "best_feasible_value": best_feasible_value,
            "last_master_solution": last_master_solution,
            "last_master_objective": last_master_objective,
            "grasp_time": other_stats['grasp_time'] if other_stats and 'grasp_time' in other_stats else None,
            "kernelization_time": other_stats['kernelization_time'] if other_stats and 'kernelization_time' in other_stats else None 
        }

    def maybe_checkpoint(reason="periodic", force=False, status="running", other_stats=None):
        nonlocal last_checkpoint_time, last_checkpoint_iteration
        if not checkpoint_path:
            return
        now = time()
        by_iter = (
            checkpoint_every_iterations is not None
            and checkpoint_every_iterations > 0
            and iteration - last_checkpoint_iteration >= checkpoint_every_iterations
        )
        by_time = (
            checkpoint_seconds is not None
            and checkpoint_seconds > 0
            and now - last_checkpoint_time >= checkpoint_seconds
        )
        if force or by_iter or by_time:
            _write_lbbd_checkpoint(
                checkpoint_path, make_stats(status), reason=reason, checkpoint=True, other_stats=other_stats
            )
            last_checkpoint_time = now
            last_checkpoint_iteration = iteration

    def _sigterm_handler(signum, frame):  
        # SLURM commonly sends SIGTERM shortly before SIGKILL.  Persist the
        # latest *certified* LB/full-feasible UB state before termination.
        try:
            maybe_checkpoint(reason=f"signal-{signum}", force=True, status="interrupted")
        finally:
            os._exit(128 + int(signum))

    signal.signal(signal.SIGTERM, _sigterm_handler)
    maybe_checkpoint(reason="initial", force=True, status="running")

    while True:
        elapsed = time() - start_time
        remaining = time_limit - elapsed
        if remaining <= 0:
            print("LBBD global time limit reached.")
            maybe_checkpoint(reason="global-time-limit", force=True, status="time_limit")
            signal.signal(signal.SIGTERM, previous_sigterm_handler)
            return best_feasible_solution, mdl, make_stats("time_limit", other_stats=other_stats)

        mdl.parameters.timelimit = remaining
        iteration += 1
        print(f"LBBD iteration {iteration}; active pair cuts: {len(added_pairs)}")

        sol = mdl.solve(log_output=log_output)
        solve_details = mdl.solve_details

        # CPLEX best_bound is the valid lower bound for the current restricted
        # master. Since the complete problem has all these constraints and
        # possibly more, this is also a valid lower bound for the full problem.
        current_lb = getattr(solve_details, "best_bound", None)
        if current_lb is not None:
            try:

                current_lb = ceil(float(current_lb)) # Round up to nearest integer for valid LB on integer objective.   
                if current_lb > best_lower_bound:
                    best_lower_bound = current_lb #update best lower bound if current LB is better  
            except (TypeError, ValueError):
                pass

        if sol is None:
            print("No master solution was returned.")
            status = "time_limit" if "time" in str(solve_details.status).lower() else "no_solution"
            maybe_checkpoint(reason="master-no-solution", force=True, status=status, other_stats=other_stats)
            signal.signal(signal.SIGTERM, previous_sigterm_handler)
            return best_feasible_solution, mdl, make_stats(status, other_stats=other_stats)

        chosen = [v for v in nodes if x[v].solution_value > 0.5]
        last_master_solution = sorted(chosen)
        last_master_objective = float(mdl.objective_value)

        # When a restricted master is solved to optimality, its optimum itself
        # is a valid lower bound and is often stronger/more stable than the
        # reported floating-point best bound.
        status_text = str(solve_details.status).lower()

        if "optimal" in status_text:
            best_lower_bound = max(best_lower_bound, last_master_objective)

        violated = find_violated_pairs(chosen, pair_sets, added_pairs)

        print(
            f"  objective: {last_master_objective:g}, "
            f"current_LB: {current_lb}, best_LB: {None if best_lower_bound == float('-inf') else best_lower_bound:g}, "
            f"violated_unadded: {len(violated)}, elapsed: {time() - start_time:.2f}s"
        )

        # Periodic checkpoint after every requested number of completed LBBD
        # iterations (or elapsed checkpoint interval).  At this point the LB is
        # certified; UB remains the best solution already verified feasible for
        # the full problem.
        maybe_checkpoint(reason="iteration", status="running")

        if not violated:
            # The current master solution satisfies all generated rows and no
            # omitted row is violated. It is therefore feasible for the full
            # formulation. Its objective also matches a valid restricted-master
            # lower bound, proving global optimality.
            best_feasible_solution = sorted(chosen)
            best_feasible_value = int(round(last_master_objective))
            best_lower_bound = max(best_lower_bound, last_master_objective)
            maybe_checkpoint(reason="optimal", force=True, status="optimal", other_stats=other_stats)
            signal.signal(signal.SIGTERM, previous_sigterm_handler)
            return best_feasible_solution, mdl, make_stats("optimal", other_stats=other_stats)

        selected_pairs = select_violated_pairs(
            violated,
            pair_sets,
            max_cuts=max_cuts_per_iteration,
            strategy=cut_strategy,
        )

        for pair in selected_pairs:
            r, b = pair
            separating_nodes = pair_sets[pair]
            mdl.add_constraint(
                mdl.sum(x[v] for v in separating_nodes) >= 1,
                ctname=f"pair_{r}_{b}",
            )
            added_pairs.add(pair)

        print(f"  added {len(selected_pairs)} cut(s): {selected_pairs}")

    
    
def compute_separation_sets(G, red_nodes, blue_nodes):

    pair_sets = {}

    for r, b in product(red_nodes, blue_nodes):

        D_rb = set()
        for v in G.nodes():

            if len({r, b} & closed_neighborhood(G, v)) == 1:
                D_rb.add(v)

        pair_sets[(r, b)] = D_rb
    return pair_sets
    
            

def find_essential_nodes(pair_sets):

    essential = set()
    
    for D_rb in pair_sets.values():

        if len(D_rb) == 1:
            essential.update(D_rb)
    
    return essential
    
    


#second pre-solving steps:
def dominated_rows(pair_sets):
    """
    Row dominance:
    Ako D_i ⊆ D_j, onda je red j dominiran i može se ukloniti.
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

            # D1 subset D2 => remove D2
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
    Column dominance.

    C(v) = skup redova koje pokriva čvor v.

    Ako C(u) ⊆ C(v),
    tada se u može ukloniti.
    """

    coverage = {
        v: set()
        for v in nodes
    }

    for pair, D in pair_sets.items():

        # Only surviving candidate vertices may contribute to coverage.
        # This prevents stale vertices left in D(r,b) from causing KeyError
        # after previous column-dominance reductions.
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

    changed = True

    removed_rows = 0
    removed_cols = 0

    while changed:

        changed = False

        # row dominance
        new_pair_sets, removed = dominated_rows(
            pair_sets
        )

        if removed:

            removed_rows += len(removed)

            pair_sets = new_pair_sets

            changed = True

        # column dominance
        dominated = dominated_columns(
            pair_sets,
            all_nodes - essential
        )

        if dominated:

            removed_cols += len(dominated)

            all_nodes -= dominated

            # Keep every separation set synchronized with the surviving
            # node set.  Using intersection is stronger than subtracting only
            # the nodes removed in the current iteration and also cleans up
            # any stale vertices left by earlier reductions.
            for pair in pair_sets:

                pair_sets[pair] &= all_nodes

            changed = True

    print(
        f"Row dominance removed {removed_rows} constraints"
    )

    print(
        f"Column dominance removed {removed_cols} variables"
    )

    return pair_sets, all_nodes
            
def main():
    
    parser = argparse.ArgumentParser(description="ILP Red-Blue Separator with CPLEX and log parsing")
    parser.add_argument("-i", "--i", required=True, help="Input CSV graph file in DIMACS format")
    parser.add_argument("-o", "--o", required=True, help="Output file path for results")
    parser.add_argument("-t", "--time_limit", type=float, default=60.0, help="CPLEX time limit in seconds")
    parser.add_argument("-warm_start", "--warm_start", type=int, default=1, help="CPLEX time limit in seconds")
    parser.add_argument(
        "--max_cuts_per_iteration",
        type=int,
        default=1,
        help="Maximum number of violated pair constraints added per LBBD iteration; use 0 for all.",
    )
    parser.add_argument(
        "--cut_strategy",
        choices=["smallest_set", "largest_set", "lexicographic", "all"],
        default="smallest_set",
        help="Rule used to select violated pair constraints.",
    )
    parser.add_argument(
        "--checkpoint_iterations",
        type=int,
        default=1,
        help="Write LB/UB/gap checkpoint every N completed LBBD iterations (default: 1; 0 disables iteration-based checkpoints).",
    )
    parser.add_argument(
        "--checkpoint_seconds",
        type=float,
        default=60.0,
        help="Also checkpoint when this many seconds have elapsed since the previous checkpoint (checked between LBBD master solves).",
    )
    
    args = parser.parse_args()

    input_file = args.i
    output_file = args.o
    time_limit = float(args.time_limit)
    warm_start = int(args.warm_start)
    
    # Read input graph
    G, R, B = read_graph(input_file)

    # IMPORTANT: detect an infeasible instance before running GRASP.
    # Otherwise the heuristic can fail on an empty candidate list even though
    # the correct outcome is simply that no red-blue separator exists.
    if not sanity_check(G, R, B):
        print("Instance is infeasible: at least one red-blue pair is unseparable.")
        with open(output_file, "w") as out_file:
            out_file.write("status: infeasible\n")
            out_file.write("initial_value: None\n")
            out_file.write("initial_solution: []\n")
            out_file.write("best_lower_bound: None\n")
            out_file.write("best_feasible_value: None\n")
            out_file.write("best_feasible_solution: None\n")
            out_file.write("last_master_objective: None\n")
            out_file.write("last_master_solution: None\n")
            out_file.write("iterations: 0\n")
            out_file.write("cuts_added: 0\n")
            out_file.write("elapsed: 0.000000\n")
        return

    # run preprocessing steps: compute separation sets, find essential nodes, and kernelize the instance     
    time_of_kernelization_start = time()   
    pair_sets = compute_separation_sets(G, R, B)
    essential = find_essential_nodes(pair_sets)
    print(f"Essential nodes: {len(essential)}")
    all_nodes = set(G.nodes())  
    pair_sets, all_nodes = kernelize(pair_sets, essential, all_nodes)
    kernelization_time_overall = time() - time_of_kernelization_start

    heuristic_solution = None
    heuristic_value = None

    if warm_start == 1:
        
        print("Running GRASP warm start with k=5, #numiter=10...")
        from LS import build_pair_cover_map
        pair_sets_grasp =  build_pair_cover_map(G, R, B)
        time_of_grasp_start = time()
        candidate = grasp_mrbssp(G=G, R=R, B=B, k=5, max_iter=10, seed=None, pair_cover_map=pair_sets_grasp, forbidden=None)  #, time_limit)  # max_iter=50
        grasp_time_overall = time() - time_of_grasp_start   

        #print(candidate)  
        grasp_feasible = checkFeasibility(candidate, G, R, B)
        print("GRASP feasible =", grasp_feasible)
    
        # A GRASP solution is a valid incumbent/upper bound ONLY when it is
        # feasible for the complete MRBSSP instance.  Passing an infeasible
        # candidate as ub would incorrectly cut off the true optimum.
        if grasp_feasible:
            heuristic_solution = list(candidate)
            heuristic_value = len(heuristic_solution)
            print(f"Valid warm-start UB = {heuristic_value}")
        else:
            print(
                "WARNING: GRASP returned an infeasible candidate. "
                "Disabling both the MIP start and the GRASP upper bound."
            )

    # Run incremental LBBD.
    result = solve_lbbd(
        all_nodes,
        pair_sets,
        essential=essential,
        warm_start=heuristic_solution,
        ub=heuristic_value,
        time_limit=time_limit,
        max_cuts_per_iteration=args.max_cuts_per_iteration,
        cut_strategy=args.cut_strategy,
        checkpoint_path=output_file,
        checkpoint_every_iterations=args.checkpoint_iterations,
        checkpoint_seconds=args.checkpoint_seconds,
        other_stats={
            "kernelization_time": kernelization_time_overall,
            "grasp_time": grasp_time_overall if warm_start == 1 else None,
        }   
    )

    solution, model, stats = result

    # Final safety guard: an "optimal" LBBD result must also satisfy the
    # original, unkernelized graph instance.  This catches any future
    # preprocessing/model inconsistency before it is reported as optimal.
    final_feasible = (
        solution is not None
        and checkFeasibility(solution, G, R, B)
    )
    if stats["status"] == "optimal" and not final_feasible:
        print(
            "WARNING: LBBD returned an 'optimal' master solution that is "
            "infeasible for the original graph. Marking result as invalid."
        )
        stats["status"] = "invalid"

    final_gap = _relative_gap(stats["best_lower_bound"], stats["best_feasible_value"])
    with open(output_file, "w") as out_file:
        out_file.write(f"status: {stats['status']}\n")
        out_file.write(f"LB: {stats['best_lower_bound']}\n")
        out_file.write(f"UB: {stats['best_feasible_value']}\n")
        out_file.write(f"gap: {final_gap}\n")
        out_file.write(f"initial_value: {stats['initial_value']}\n")
        out_file.write(f"initial_solution: {stats['initial_solution']}\n")
        out_file.write(f"best_lower_bound: {stats['best_lower_bound']}\n")
        out_file.write(f"best_feasible_value: {stats['best_feasible_value']}\n")
        out_file.write(f"best_feasible_solution: {stats['best_feasible_solution']}\n")
        out_file.write(f"last_master_objective: {stats['last_master_objective']}\n")
        out_file.write(f"last_master_solution: {stats['last_master_solution']}\n")
        out_file.write(f"iterations: {stats['iterations']}\n")
        out_file.write(f"cuts_added: {stats['cuts_added']}\n")
        out_file.write(f"elapsed: {stats['elapsed']:.6f}\n")
        out_file.write(f"final_feasible: {final_feasible}\n")
        out_file.write("checkpoint: False\n")
        out_file.write("checkpoint_reason: final\n")
        out_file.write(f"kernelization_time: {stats['kernelization_time']}\n")
        out_file.write(f"grasp_time: {stats['grasp_time']}\n")
    

if __name__=="__main__":
    print("LBBD solver scaffold. Integrate read_graph(), pair_sets and kernelization from your project.")
    
    main()
    

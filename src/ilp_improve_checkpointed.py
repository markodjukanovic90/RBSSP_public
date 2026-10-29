
import argparse
import csv
from pathlib import Path
from time import time
from itertools import product
from collections import defaultdict
from os import name as os_name, remove

import networkx as nx
import pulp as pl
#import orloge  # pip install orloge
from math import ceil
from grasp import grasp_mrbssp

# Detect platform and set CPLEX path
import sys
#sys.path.append(r'/ceph/hpc/home/djukanovicm/VGLCS/cplex/bin/x86-64_linux/cplex')
#"/home/marko/Desktop/cplex/cplex/python/3.10/x86-64_linux")

import sys

sys.path.insert(
    0,
    "/ceph/hpc/home/djukanovicm/VGLCS/cplex/python/3.10/x86-64_linux"
)

import cplex
import os
import signal
import math

print("CPLEX:", cplex.__version__)
print("Loaded from:", cplex.__file__)

## cplex link 

# ---------------------- Graph & Model Utilities ---------------------------- #

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

            # ukloni dominirane čvorove iz svih redova
            for pair in pair_sets:

                pair_sets[pair] -= dominated

            changed = True

    print(
        f"Row dominance removed {removed_rows} constraints"
    )

    print(
        f"Column dominance removed {removed_cols} variables"
    )

    return pair_sets, all_nodes
    


    
    

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

def retrieve_results_from_log1(path_to_log):
    log_data = orloge.get_info_solver(path_to_log, "CPLEX")
    best_bound = ceil(log_data.get("best_bound", 0))
    best_solution = log_data.get("best_solution", 0)
    status = int(log_data.get("status") == "MIP - Integer optimal")
    remove(path_to_log)  # Clean up log file
    return best_solution, best_bound, status
    

def retrieve_results_from_log(path_to_log):

    log_data = orloge.get_info_solver(
        path_to_log,
        "CPLEX"
    )

    print("\n===== LOG DATA =====")
    print(log_data)
    print("====================\n")

    best_bound = log_data.get("best_bound")
    best_solution = log_data.get("best_solution")

    if best_bound is None:
        best_bound = 0

    if best_solution is None:
        best_solution = 0

    status = int(
        log_data.get("status")
        == "MIP - Integer optimal"
    )

    remove(path_to_log)

    return (
        best_solution,
        ceil(best_bound),
        status
    )

# ---------------------- Main Execution Block ------------------------------ #
########################################## ILP model generation #####################################################
def generate_ilp_model(
        G,
        red_nodes,
        blue_nodes,
        pair_sets,
        essential, 
        sol_warm_start):

    remaining_nodes = set(G.nodes()) - essential
    print(sol_warm_start)
    prob = pl.LpProblem(
        "MinRedBlueSeparatingSet",
        pl.LpMinimize
    )

    x = {
        v: pl.LpVariable(
            f"x_{v}",
            cat="Binary"
        )
        for v in remaining_nodes
    }
    
    for v in remaining_nodes:
    
        if v in sol_warm_start:
            x[v].setInitialValue(1)
        else:
            x[v].setInitialValue(0)
    
    for v in essential:
        prob += x[v] == 1
    
        
    prob += (
        len(essential)    +
        pl.lpSum(x[v] for v in remaining_nodes)
    ), "Minimize_Separator_Size"
    
    prob += (
    pl.lpSum(
        x[v]
        for v in remaining_nodes
    )  <=  len(sol_warm_start) )#- len(essential))
    
    for D_rb in pair_sets.values():

        # već pokriven essential čvorom
        if len(D_rb & essential) > 0:
            continue

        D_remaining = D_rb & remaining_nodes

        if len(D_remaining) > 0:

            prob += (
                pl.lpSum(x[v] for v in D_remaining)
                >= 1
            )
    
    #for name, c in list(prob.constraints.items())[:10]:
    #    print(c)

    return prob, x, remaining_nodes
    

from docplex.mp.model import Model


def generate_docplex_model(
        G,
        pair_sets,
        essential,
        warm_start_sol, all_nodes):

    mdl = Model(name="MRBSSP")

    remaining_nodes = all_nodes - essential
    #remaining_nodes = set(G.nodes()) - essential

    x = {
        v: mdl.binary_var(
            name=f"x_{v}"
        )
        for v in remaining_nodes
    }

    mdl.minimize(
        len(essential)
        +
        mdl.sum(
            x[v]
            for v in remaining_nodes
        )
    )
    
    # obj UB constraint
    if warm_start_sol != None:
        mdl.add_constraint(mdl.sum( x[v] for v in remaining_nodes) <= len(warm_start_sol) - len(essential))
   

    for D_rb in pair_sets.values():

        if D_rb & essential:
            continue

        D_remaining = D_rb & remaining_nodes

        if D_remaining:

            mdl.add_constraint(
                mdl.sum(
                    x[v]
                    for v in D_remaining
                ) >= 1
            )
    
    #presolve -- provide a warm start 
    mip_start = mdl.new_solution()

    for v in remaining_nodes:

        if warm_start_sol != None and v in warm_start_sol:

            mip_start.add_var_value(
                x[v], 1
            )

        else:

            mip_start.add_var_value(
                x[v],
                0
            )

    mdl.add_mip_start(mip_start)

    return mdl, x, remaining_nodes
    
    
# ---------------------- Checkpoint support --------------------------------- #

def _safe_float(value, default=None):
    try:
        if value is None:
            return default
        value = float(value)
        if math.isnan(value) or math.isinf(value):
            return default
        return value
    except (TypeError, ValueError, OverflowError):
        return default


def _relative_gap(ub, lb):
    ub = _safe_float(ub)
    lb = _safe_float(lb)
    if ub is None or lb is None:
        return None
    if abs(ub) <= 1e-12:
        return 0.0 if abs(ub - lb) <= 1e-9 else None
    return max(0.0, 100.0 * (ub - lb) / abs(ub))


def write_checkpoint(output_file, state, checkpoint=True, reason="periodic", other_stats=None ):
    """Atomically write the latest valid ILP state.

    Atomic replacement means that a scheduler kill during a write leaves either
    the previous complete checkpoint or the new complete checkpoint, rather than
    a partially written result file.
    """
    output_file = str(output_file)
    parent = os.path.dirname(os.path.abspath(output_file))
    if parent:
        os.makedirs(parent, exist_ok=True)

    tmp_file = f"{output_file}.tmp.{os.getpid()}"
    ub = state.get("UB")
    lb = state.get("LB")
    gap = _relative_gap(ub, lb)

    with open(tmp_file, "w") as f:
        f.write(f"Solution: {state.get('Solution', [])}\n")
        f.write(f"LB: {ceil(lb) if lb is not None else -1}\n")
        f.write(f"UB: {ub if ub is not None else -1}\n")
        f.write(f"Gap: {gap if gap is not None else -1}\n")
        f.write(f"Status: {state.get('Status', 'RUNNING')}\n")
        f.write(f"Feasible: {state.get('Feasible', False)}\n")
        f.write(f"Time: {state.get('Time', 0.0):.6f}\n")
        f.write(f"Initial: {state.get('Initial', -1)}\n")

        #GRASP solution (for warmstart) stats:
        f.write(f"InitialSolutionValue: {state.get('InitialSolutionValue', -1)}\n")
        f.write(f"InitialSolutionNodes: {state.get('InitialSolutionNodes', [])}\n")     

        f.write(f"NodesProcessed: {state.get('NodesProcessed', -1)}\n")
        f.write(f"NodesLeft: {state.get('NodesLeft', -1)}\n")
        f.write(f"Checkpoint: {checkpoint}\n")
        f.write(f"CheckpointReason: {reason}\n")

        # preprocessing stats (GRASP, kernelization, model generation): TODO 
        if other_stats:
            for key, value in other_stats.items():
                f.write(f"{key}: {value}\n")            
                

    os.replace(tmp_file, output_file)


class ILPCheckpointCallback(cplex.callbacks.MIPInfoCallback):
    """Periodic MIP progress checkpoint for the compact ILP.

    The compact ILP contains the complete constraint set, so every CPLEX
    incumbent is a valid primal upper bound for MRBSSP and CPLEX's best
    objective bound is a valid lower bound.  No incumbent is fabricated from
    fractional or otherwise incomplete information.
    """

    def __call__(self):  
        now = time()
        elapsed = now - self.start_wall

        try:
            nodes = int(self.get_num_nodes())
        except Exception:
            nodes = -1
        try:
            nodes_left = int(self.get_num_remaining_nodes())
        except Exception:
            nodes_left = -1

        try:
            lb = _safe_float(self.get_best_objective_value())
        except Exception:
            lb = None
        if lb is not None:
            old_lb = _safe_float(self.state.get("LB"))
            if old_lb is None or lb > old_lb:
                self.state["LB"] = lb

        try:
            has_incumbent = bool(self.has_incumbent())
        except Exception:
            has_incumbent = False

        if has_incumbent:
            try:
                ub = _safe_float(self.get_incumbent_objective_value())
            except Exception:
                ub = None
            if ub is not None:
                old_ub = _safe_float(self.state.get("UB"))
                if old_ub is None or ub < old_ub:
                    self.state["UB"] = ub
                # An incumbent of the full compact ILP is feasible by definition.
                self.state["Feasible"] = True

        self.state["Status"] = "RUNNING"
        self.state["Time"] = elapsed
        self.state["NodesProcessed"] = nodes
        self.state["NodesLeft"] = nodes_left

        due_time = (
            self.checkpoint_seconds > 0
            and elapsed - self.last_checkpoint_time >= self.checkpoint_seconds
        )
        due_nodes = (
            self.checkpoint_nodes > 0
            and nodes >= 0
            and (self.last_checkpoint_nodes < 0 or
                 nodes - self.last_checkpoint_nodes >= self.checkpoint_nodes)
        )

        if due_time or due_nodes:
            reason = "seconds" if due_time else "nodes"
            try:
                write_checkpoint(self.output_file, self.state, True, reason)
                self.last_checkpoint_time = elapsed
                self.last_checkpoint_nodes = nodes
            except Exception as exc:
                # A checkpoint failure must never abort the optimization.
                print(f"WARNING: could not write ILP checkpoint: {exc}", flush=True)

    def setup(
        self,
        output_file,
        state,
        start_wall,
        checkpoint_seconds=60.0,
        checkpoint_nodes=1000,
    ):
        self.output_file = output_file
        self.state = state
        self.start_wall = start_wall
        self.checkpoint_seconds = max(0.0, float(checkpoint_seconds))
        self.checkpoint_nodes = max(0, int(checkpoint_nodes))
        self.last_checkpoint_time = 0.0
        self.last_checkpoint_nodes = -1


# ---------------------- Main Execution Block ------------------------------ #

def main():
    parser = argparse.ArgumentParser(
        description="ILP Red-Blue Separator with periodic UB/LB checkpoints"
    )
    parser.add_argument("-i", "--i", required=True, help="Input CSV graph file in DIMACS format")
    parser.add_argument("-o", "--o", required=True, help="Output file path for results/checkpoints")
    parser.add_argument("-t", "--time_limit", type=float, default=1200.0, help="CPLEX time limit in seconds")
    parser.add_argument("-warm", "--warm_start", type=int, default=1, help="Use GRASP warm start (1=yes, 0=no)")
    parser.add_argument(
        "--checkpoint_seconds",
        type=float,
        default=60.0,
        help="Write current UB/LB to -o at least approximately every N seconds (0 disables).",
    )
    parser.add_argument(
        "--checkpoint_nodes",
        type=int,
        default=1000,
        help="Also checkpoint after approximately N additional B&B nodes (0 disables).",
    )
    args = parser.parse_args()

    input_file = args.i
    output_file = args.o
    time_limit = float(args.time_limit)
    warm_start = int(args.warm_start)
    overall_start = time()
    solve_start = None #time when the algorithm starts solving the ILP (after preprocessing and warm start) 

    state = {
        "Solution": [],
        "LB": None,
        "UB": None,
        "Status": "INITIALIZING",
        "Feasible": False,
        "Time": 0.0,
        "Initial": -1,
        "NodesProcessed": -1,
        "NodesLeft": -1,
    }

    # Best-effort final snapshot when SLURM sends SIGTERM before SIGKILL.
    def _handle_signal(signum, frame, other_stats=None):
        state["Status"] = "INTERRUPTED"
        state["TimeOverall"] = time() - overall_start
        state["Time"] = time() - solve_start

        try:
            write_checkpoint(output_file, state, True, f"signal-{signum}", other_stats=other_stats  )
        finally:
            raise SystemExit(128 + signum)

    try:
        signal.signal(signal.SIGTERM, _handle_signal)
        signal.signal(signal.SIGINT, _handle_signal)
    except Exception:
        pass

    G, R, B = read_graph(input_file)  

    # Check structural feasibility before GRASP, so infeasible instances still
    # leave a meaningful output file rather than failing inside the heuristic.
    sanity = sanity_check(G, R, B)
    if not sanity:
        state.update({
            "Status": "Infeasible instance (sanity check failed)",
            "LB": -1,
            "UB": -1,
            "Time": time() - overall_start,
        })
        write_checkpoint(output_file, state, False, "final-infeasible")
        return

    heuristic_solution = None
    heuristic_value = None
    grasp_time = kernel_time = model_time = 0.0

    if warm_start == 1:
        print("############ Running GRASP warm start ... #################  ", flush=True)
        try:
            from LS import build_pair_cover_map
            pair_sets_grasp =  build_pair_cover_map(G, R, B)
            grasp_start = time()
            candidate = grasp_mrbssp(G=G, R=R, B=B, k=5, max_iter=10, seed=None, pair_cover_map=pair_sets_grasp, forbidden=None ) #, time_limit)
            grasp_time = time() - grasp_start   

        except Exception as exc:
            print(f"WARNING: GRASP failed ({type(exc).__name__}: {exc}). Continuing without warm start.", flush=True)
            candidate = None

        if candidate is not None and checkFeasibility(candidate, G, R, B):
            heuristic_solution = list(candidate)
            heuristic_value = len(heuristic_solution)
            state["UB"] = float(heuristic_value)
            state["Solution"] = sorted(heuristic_solution)
            state["Feasible"] = True
            state["Initial"] = heuristic_value
            print("GRASP feasible = True", flush=True)
            print(f"Warm start solution size = {heuristic_value}", flush=True)
        else:
            print("GRASP feasible = False; warm start/UB disabled.", flush=True)

    #kernelization and essential node detection
    time_kernel_start = time()      
    pair_sets = compute_separation_sets(G, R, B)
    essential = find_essential_nodes(pair_sets)
    print(f"###############Essential nodes: {len(essential)}", flush=True)

    all_nodes = set(G.nodes())
    pair_sets, all_nodes = kernelize(pair_sets, essential, all_nodes)
    kernel_time = time() - time_kernel_start

    print("Remaining nodes after kernelization:", len(all_nodes), flush=True)
    print("Remaining rows:", len(pair_sets), flush=True)

    # A trivial, always-valid lower bound after essential-node detection.
    state["LB"] = float(len(essential))
    state["Status"] = "READY"
    state["Time"] = time() - overall_start

    state["InitialSolutionValue"] = heuristic_value if heuristic_value is not None else -1
    state["InitialSolutionNodes"] = sorted(heuristic_solution) if heuristic_solution is not None else []    

    write_checkpoint(output_file, state, True, "before-solve")

    model_time_start = time()
    mdl, x_vars, remaining_nodes = generate_docplex_model(
        G, pair_sets, essential, heuristic_solution, all_nodes
    )
    model_time = time() - model_time_start  

    # model parameters for CPLEX    
    mdl.parameters.timelimit = time_limit
    mdl.parameters.threads = 1

    #stats of preprocessinf steps (GRASP, kernelization, model generation):
    other_stats = {
        "grasp_time": grasp_time,
        "kernel_time": kernel_time,
        "model_time": model_time
    }

    #Model solving with periodic checkpointing of UB/LB and other stats (GRASP, kernelization, model generation times)
    solve_start = time()  
    #overall_start = time()
    cb = mdl.register_callback(ILPCheckpointCallback)
    cb.setup(
        output_file=output_file,
        state=state,
        start_wall=solve_start,
        checkpoint_seconds=args.checkpoint_seconds,
        checkpoint_nodes=args.checkpoint_nodes,
    )


    #cplex solve:
    #time_cplex_start = time()
    sol = mdl.solve(log_output=True)
    total_time = round(time() - solve_start, 2)
    #time_ilp_approach_run = round(time() - time_cplex_start, 2)
    details = mdl.solve_details

    # Always preserve the strongest final CPLEX lower bound if available.
    final_lb = _safe_float(getattr(details, "best_bound", None))

    if final_lb is not None:
        final_lb = ceil(final_lb)
        old_lb = _safe_float(state.get("LB"))
        if old_lb is None or final_lb > old_lb:
            state["LB"] = final_lb

    if sol:

        ilp_solution = list(essential)
        for v in remaining_nodes:
            if x_vars[v].solution_value > 0.5:
                ilp_solution.append(v)

        ub = _safe_float(mdl.objective_value)
        # The independently verified GRASP solution remains a valid fallback UB.
        if heuristic_value is not None and (ub is None or heuristic_value < ub):
            ub = float(heuristic_value)
            ilp_solution = list(heuristic_solution)

        state["UB"] = ub
        state["Solution"] = sorted(ilp_solution)
        #GRASP solution for warm-start stats:
        state["InitialSolutionValue"] = heuristic_value if heuristic_value is not None else -1
        state["InitialSolutionNodes"] = sorted(heuristic_solution) if heuristic_solution is not None else []    
    
        state["Feasible"] = checkFeasibility(ilp_solution, G, R, B)
        state["Status"] = str(details.status)
    else:
        # Even if mdl.solve() returns no solution, a periodic callback may have
        # already recorded a valid incumbent and a useful certified lower bound.
        state["Status"] = str(getattr(details, "status", "No solution"))

    state["Time"] = total_time

    try:
        state["NodesProcessed"] = int(getattr(details, "nb_nodes_processed", state["NodesProcessed"]))
        #state["TimeTotal"] = time_ilp_approach_run
    except Exception:
        pass

    write_checkpoint(output_file, state, False, "final", other_stats)


# ---------------------------- Entry Point --------------------------------- #

if __name__ == "__main__":
    main()

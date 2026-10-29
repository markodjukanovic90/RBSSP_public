import networkx as nx
import random
from greedy import randomized_greedy_separator, load_graph_from_dimacs_csv
from LS import local_search, build_pair_cover_map 
import argparse
from time import time 

# -----------------------------------
# Full GRASP Algorithm  
# -----------------------------------

def grasp_mrbssp(G, R, B, k=5, max_iter=10, seed=None, tx=1000, pair_cover_map=None):#@k does not need to be fixed, but varied according to @alpha param. and min and max greedy value
    """
    GRASP for Minimum Red-Blue Neighborhood Separator.
    """
    if seed is not None:
        random.seed(seed)
    
    #pair_cover_map = build_pair_cover_map(G, R, B)
    #print(pair_cover_map)
    best_solution = None
    best_cost = float('inf') 
    start_time = time()

    for _ in range(max_iter): 
    
        S_greedy = randomized_greedy_separator(G, R, B, k)
        S_local = local_search(G, R, B, S_greedy, pair_cover_map) #try to improve S_greedy 

        if len(S_local) < best_cost: 
            best_solution = S_local
            best_cost = len(S_local)
        
        end_time = time()

        if end_time - start_time >= tx:
            break

    return best_solution

# -----------------------------------
# Example Usage
# -----------------------------------

if __name__ == "__main__":  #grasp-k-3-maxiter-150


    # terminal parameters:
    parser = argparse.ArgumentParser(description="Optimized greedy red-blue separator")
    parser.add_argument("-i", "--i", help="Path to DIMACS-style CSV input graph")
    parser.add_argument("-o", "--o", help="Path to output file where results will be stored")
    parser.add_argument("-maxiter", "--maxiter", default="1000", help="Number of interations in GRASP") #default = 10
    parser.add_argument("-k", "--k", default="3", help="k in the randomized Greedy procedure")
    parser.add_argument("-t", "--t", default="1200", help="k in the randomized Greedy procedure")

    args = parser.parse_args()

    # Load graph
    G = load_graph_from_dimacs_csv(args.i)
    red_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "red"}
    blue_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "blue"}
    
    print("Red: ", red_nodes)
    print("Blue: ", blue_nodes)
    print("Max. iterations:", args.maxiter.strip())
    # Run GRASP heuristic
    start = time()
    #
    pair_cover_map = build_pair_cover_map(G, red_nodes, blue_nodes)
    final_solution = grasp_mrbssp(G, red_nodes, blue_nodes, int(args.k), int(args.maxiter.strip()), None, int(args.t), pair_cover_map)
    end = time()
    duration = round(end - start, 2)

    # Output results
    with open(args.o, "w") as f:
        f.write(f"Solution: {final_solution}\n")         # Solution
        f.write(f"Quality: {len(final_solution)}\n")    # Solution quality (size)
        f.write(f"Time: {duration}\n")          # Execution time

    print(f"Separator written to {args.o}")
    



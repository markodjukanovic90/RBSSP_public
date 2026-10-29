import networkx as nx
import csv
import argparse
from itertools import product
from collections import defaultdict
from time import time
import random 
#from draw_graph import draw_graph

def closed_neighborhood_opt(G, v):
    return set(G.neighbors(v)).union({v})
    
def sanity_check(G, pairs):

    for r, b in pairs:
        if not any((r in closed_neighborhood_opt(G, v)) != (b in closed_neighborhood_opt(G, v)) for v in G.nodes()):
            print(f"Unseparable pair: ({r}, {b})")
            return False
    
    return True 

####################### optimized GREEDY version ##############################################################


def optimized_greedy_separator(G, red, blue):


    #draw the graph
    #draw_graph(G, red, blue)

    pairs = set(product(red, blue))
    separator = set()
    node_separates = defaultdict(set)
    
    #sanity_check (verify if at least one feasible solution is possible):
    separable=sanity_check(G, pairs)  
    if not separable:
        return [] # no fasible solution

    for v in G.nodes():
        N_v = closed_neighborhood_opt(G, v)
        for (r, b) in pairs:
            if (r in N_v) != (b in N_v):
                node_separates[v].add((r, b))

    while pairs:
        best_node = max(node_separates.keys(), key=lambda v: len(node_separates[v] & pairs))
        covered = node_separates[best_node] & pairs
        separator.add(best_node)
        pairs -= covered

        del node_separates[best_node]
        for v in node_separates:
            node_separates[v] -= covered
  
    return separator

def load_graph_from_dimacs_csv(path):
    G = nx.Graph()
    with open(path, "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == "n":
                G.add_node(int(row[1]), color=row[2])
            elif row[0] == "e":
                G.add_edge(int(row[1]), int(row[2]))
    return G

    
############################# RANDOMIZED greedy approach ##############################################


def randomized_greedy_separator(G, red, blue, k=3, seed=None):

    #random.seed(seed)
    pairs = set(product(red, blue))
    separator = set()
    node_separates = defaultdict(set)

    # Precompute which node can separate which red-blue pairs
    for v in G.nodes():
        N_v = closed_neighborhood_opt(G, v)
        for (r, b) in pairs:
            if (r in N_v) != (b in N_v):  # XOR: one in, one out
                node_separates[v].add((r, b))
    
    # While there are unresolved red-blue pairs
    while pairs:
        cover_counts = []
        for v in node_separates:
            cover = node_separates[v] & pairs
            cover_counts.append((v, len(cover)))
        
        #if len(cover_counts) == 0:
        #    return []
        
        # Sort by cover size in descending order
        cover_counts.sort(key=lambda x: x[1], reverse=True)

        # Build Candidate List with top-k nodes
        candidate_list = [v for (v, _) in cover_counts[:min(k, len(cover_counts)) ]]
        
        # Randomly select one candidate to add to separator
        chosen_v = random.choice(candidate_list) #chosen_v selected
        covered = node_separates[chosen_v] & pairs
        separator.add(chosen_v)
        pairs -= covered

        # Clean up: remove selected node and update remaining covers
        del node_separates[chosen_v]
        for v in node_separates:
            node_separates[v] -= covered

    return separator


######################################################################################################

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Optimized greedy red-blue separator")
    parser.add_argument("-i", "--i", help="Path to DIMACS-style CSV input graph")
    parser.add_argument("-o", "--o", help="Path to output file where results will be stored")
    args = parser.parse_args()

    # Load graph
    G = load_graph_from_dimacs_csv(args.i)
    red_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "red"}
    blue_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "blue"}

    # Run heuristic
    start = time()
    #separator = optimized_greedy_separator(G, red_nodes, blue_nodes)
    separator=randomized_greedy_separator(G, red_nodes, blue_nodes, 10)
    end = time()
    duration = round(end - start, 2)
    
    # Output results
    with open(args.o, "w") as f:
        f.write(f"Solution: {separator}\n")         # Solution
        f.write(f"Quality: {len(separator)}\n")    # Solution quality (size)
        f.write(f"Time: {duration}\n")          # Execution time

    print(f"Separator written to {args.o}")
    
    

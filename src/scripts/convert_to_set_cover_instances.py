import networkx as nx
import csv
import argparse
from itertools import product
from collections import defaultdict
from time import time
import random 
from draw_graph import draw_graph

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

    
############################# Instance transformator ##############################################

def convert_pairs_into_num(pair, num_vertices):
    return pair[0]*num_vertices + pair[1]

def create_set_cover_instance_from_mrbssp(G, red, blue):

    pairs = set(product(red, blue))
    pairs_enumerate=list(enumerate(pairs))

    separator = set()
    node_separates = defaultdict(set)

    # Precompute which node can separate which red-blue pairs
    for index_of_pair, (r, b) in pairs_enumerate:
        for v in G.nodes():
            N_v = closed_neighborhood_opt(G, v)
            if (r in N_v) != (b in N_v):  # XOR: one in, one out
                node_separates[ v ].add( index_of_pair+1 ) 
                
    return node_separates, len(pairs), len(G.nodes())
    
def create_set_cover_instance_from_mrbssp1(G, red, blue): # S_i: pairs covered by i 

    pairs = set(product(red, blue))
    pairs_enumerate=list(enumerate(pairs))

    separator = set()
    node_separates = defaultdict(set)

    # Precompute which node can separate which red-blue pairs
    for index_of_pair, (r, b) in pairs_enumerate:
        for v in G.nodes():
            N_v = closed_neighborhood_opt(G, v)
            if (r in N_v) != (b in N_v):  # XOR: one in, one out
                node_separates[ index_of_pair+1 ].add(v+1)  #   vertices v  => list of pairs covered by v
                  
    return node_separates, len(pairs), len(G.nodes()) # |S|=|V(G)|
    
def save_set_cover_instance_into_file(file_path, node_separates, num_U, num_S): #num_S, num_U
    # Output results
    with open(file_path, "w") as f:
        f.write(f" {num_U} {num_S} \n")      # first line {num_S} {num_U}
        #the second line, all 1s
        for _ in range(num_S): #num_
            f.write(" 1") 

        f.write(" \n")
        #write S_i and its elements:
        for u, S_u in node_separates.items(): #u-pair
            f.write(" " + str(len(S_u)) + "\n")   
            for e in S_u:
                f.write(" " + str(e) )
                
            f.write(" \n")
        
         
    print(f"Separator written to {args.o}")
    
######################################################################################################

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Optimized greedy red-blue separator")
    parser.add_argument("-i", "--i", help="Path to DIMACS-style CSV input graph")
    parser.add_argument("-o", "--o", help="Path to the output set cover instance to be stored")
    args = parser.parse_args()

    # Load graph
    G = load_graph_from_dimacs_csv(args.i)
    red_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "red"}
    blue_nodes = {v for v, d in G.nodes(data=True) if d.get("color") == "blue"}

    separator, len_U, len_S = create_set_cover_instance_from_mrbssp1(G, red_nodes, blue_nodes)

    #save into a file:
    save_set_cover_instance_into_file(args.o, separator, len_U, len_S)
    

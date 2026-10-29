from collections import defaultdict
from itertools import product

## we do this just once  (r, b)==> set of all vertices that resolves this pair
def build_pair_cover_map(G, red_nodes, blue_nodes):
    """Precompute which nodes can separate each red-blue pair."""
    pair_cover_map = defaultdict(set)
    for v in G.nodes():
        N_v = set(G.neighbors(v)).union({v})
        for r in red_nodes:
            for b in blue_nodes:
                if (r in N_v) != (b in N_v):  # exclusive OR
                    pair_cover_map[(r, b)].add(v)
    return pair_cover_map

## cover for solution @sol
def get_initial_covered_pairs(sol, pair_cover_map):
    """Return set of covered (r,b) pairs for initial solution."""
    covered = set()
    for (r, b), nodes in pair_cover_map.items():
        if nodes & sol:
            covered.add((r, b))
    return covered

def is_complete(covered_pairs, total_pairs):
    return len(covered_pairs) == len(total_pairs)

#LS procedure with two neighborhood strucutre (1-remove and 1-swap); Comment: 1-swap now removed (can be used to convert from infeasible to feasible solution -- VNS?)
def local_search(G, red_nodes, blue_nodes, solution, pair_cover_map=None):
    """
    Local search with 1-removal and 1-swap using incremental coverage tracking.
    Input: G, sets R, B, and initial feasible set `solution`
    Returns: improved feasible set
    """
    #pair_cover_map = build_pair_cover_map(G, red_nodes, blue_nodes) #do not need to be always called: TODO: optimize
    total_pairs = set(product(red_nodes, blue_nodes))
    solution = set(solution)  # ensure it's a set

    # Compute current coverage state
    covered_pairs = get_initial_covered_pairs(solution, pair_cover_map)

    improvement = True
    while improvement:
        
        improvement = False

        # 1-Removal step
        for v in list(solution):
            new_covered = set(covered_pairs)
            for (r, b) in total_pairs:
                if v in pair_cover_map[(r, b)]:
                    # Check if v is the only one covering it
                    if solution.intersection(pair_cover_map[(r, b)]) == {v}:
                        new_covered.discard((r, b))

            if is_complete(new_covered, total_pairs):
                solution.remove(v)
                covered_pairs = new_covered
                improvement = True
                break  # restart local search

        #if improvement:
        #    continue  # restart after improvement

        # 1-Swap step
        '''
        for v_out in list(solution):
            for v_in in set(G.nodes()) - solution:
                new_solution = (solution - {v_out}) | {v_in}
                new_covered = set()

                for (r, b), nodes in pair_cover_map.items():
                    if new_solution & nodes:
                        new_covered.add((r, b))

                if is_complete(new_covered, total_pairs) and len(new_solution) <= len(solution):
                    solution = new_solution
                    covered_pairs = new_covered
                    improvement = True
                    break
            if improvement:
                break
        '''
    return solution
    


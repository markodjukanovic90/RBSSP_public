from itertools import product
from collections import defaultdict
import networkx as nx
import random
from greedy import randomized_greedy_separator, load_graph_from_dimacs_csv
from LS import local_search 
import argparse
from time import time 


# ---------------------------------------------------
# Closed neighborhood
# ---------------------------------------------------

def closed_neighborhood_opt(G, v):
    return set(G.neighbors(v)).union({v})

# ---------------------------------------------------
# Improved randomized greedy
# ---------------------------------------------------

def randomized_greedy_separator(
        G,
        red,
        blue,
        alpha=None):

    pairs = set(product(red, blue))

    separator = set()

    node_separates = defaultdict(set)

    # which node separates which pair
    for v in G.nodes():

        Nv = closed_neighborhood_opt(G, v)

        for r, b in pairs:

            if (r in Nv) != (b in Nv):

                node_separates[v].add((r, b))

    # pair frequencies
    print("######## Pair frequencies calculated ... ###################################### ")   
    pair_frequency = defaultdict(int)

    for v in node_separates:

        for p in node_separates[v]:

            pair_frequency[p] += 1

    while pairs:

        scores = []

        for v in node_separates:

            covered = node_separates[v] & pairs

            if len(covered) == 0:
                continue

            # rare pair bonus
            rare_bonus = sum(
                1.0 / pair_frequency[p]
                for p in covered
            )

            # objective
            score = (
                10.0 * len(covered)
                + rare_bonus
                - 0.05 * G.degree(v)
            )

            scores.append((v, score))

        if len(scores) == 0:
            break

        scores.sort(
            key=lambda x: x[1],
            reverse=True
        )

        if alpha is None:

            alpha = random.uniform(
                0.05,
                0.35
            )

        best_score = scores[0][1]
        worst_score = scores[-1][1]

        threshold = best_score - alpha * (
            best_score - worst_score
        )

        RCL = [
            v
            for v, score in scores
            if score >= threshold
        ]

        chosen = random.choice(RCL)

        covered = node_separates[chosen] & pairs

        separator.add(chosen)

        pairs -= covered

        del node_separates[chosen]

        for v in node_separates:

            node_separates[v] -= covered

    return separator

# ---------------------------------------------------
# Fast local search
# ---------------------------------------------------

def local_search_improved(
        G,
        red,
        blue,
        solution):

    current = set(solution)

    improved = True

    while improved:

        improved = False

        # 1-delete neighborhood

        for v in list(current):

            candidate = current - {v}

            if is_feasible(
                    G,
                    red,
                    blue,
                    candidate):

                current = candidate

                improved = True

                break

    return current

# ---------------------------------------------------
# Feasibility check
# ---------------------------------------------------

def is_feasible(
        G,
        red,
        blue,
        separator):

    for r in red:

        for b in blue:

            separated = False

            for v in separator:

                Nv = closed_neighborhood_opt(G, v)

                if (r in Nv) != (b in Nv):

                    separated = True

                    break

            if not separated:

                return False

    return True

# ---------------------------------------------------
# GRASP
# ---------------------------------------------------

def grasp_mrbssp(
        G,
        red,
        blue,
        max_iter=300,
        seed=None):

    if seed is not None:

        random.seed(seed)

    best_solution = None

    best_cost = float("inf")

    elite_pool = []

    for it in range(max_iter):

        alpha = random.uniform(
            0.05,
            0.35
        )

        S0 = randomized_greedy_separator(
            G,
            red,
            blue,
            alpha
        )

        S = local_search_improved(
            G,
            red,
            blue,
            S0
        )

        cost = len(S)

        if cost < best_cost:

            best_solution = S.copy()

            best_cost = cost

        elite_pool.append(S.copy())

        elite_pool.sort(key=len)

        if len(elite_pool) > 10:

            elite_pool.pop()

    return best_solution
    
    
if __name__ == "__main__":


    import argparse
    from time import time

    parser = argparse.ArgumentParser(
        description="GRASP for Minimum Red-Blue Separating Set Problem"
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="DIMACS CSV graph"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file"
    )

    parser.add_argument(
        "-maxiter",
        "--maxiter",
        type=int,
        default=300,
        help="Number of GRASP iterations"
    )

    parser.add_argument(
        "-seed",
        "--seed",
        type=int,
        default=None,
        help="Random seed"
    )

    args = parser.parse_args()

    # ----------------------------------
    # Load graph
    # ----------------------------------

    G = load_graph_from_dimacs_csv(
        args.input
    )

    red_nodes = {
        v
        for v, d in G.nodes(data=True)
        if d.get("color") == "red"
    }

    blue_nodes = {
        v
        for v, d in G.nodes(data=True)
        if d.get("color") == "blue"
    }

    print(
        f"Graph loaded: "
        f"{G.number_of_nodes()} nodes, "
        f"{G.number_of_edges()} edges"
    )

    print(
        f"Red={len(red_nodes)} "
        f"Blue={len(blue_nodes)}"
    )

    # ----------------------------------
    # Run GRASP
    # ----------------------------------

    start = time()

    best_solution = grasp_mrbssp(
        G,
        red_nodes,
        blue_nodes,
        max_iter=args.maxiter,
        seed=args.seed
    )

    end = time()

    runtime = round(
        end - start,
        4
    )

    # ----------------------------------
    # Output
    # ----------------------------------

    quality = len(best_solution)

    print("\n===== BEST SOLUTION =====")
    print("Separator size:", quality)
    print("Separator:", sorted(best_solution))
    print("Runtime:", runtime)

    with open(args.output, "w") as f:

        f.write(
            f"Quality: {quality}\n"
        )

        f.write(
            f"Time: {runtime}\n"
        )

        f.write(
            "Solution: "
            + " ".join(
                map(str,
                    sorted(best_solution))
            )
            + "\n"
        )

    print(
        f"\nResults written to "
        f"{args.output}"
    )

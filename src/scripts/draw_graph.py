import networkx as nx
import matplotlib.pyplot as plt

def draw_graph(G, R, B):

    # Define color mapping
    color_map = {
        'R': 'red',
        'B': 'blue'
    }
    
    # Extract node colors
    node_colors = [None] * len(G.nodes()) #[color_map[G.nodes[n]['color']] for n in G.nodes()]
    
    for r in R:   
        node_colors[r] = 'red'
    for b in B:
        node_colors[b] = 'blue'
    
    # Draw
    pos = nx.spring_layout(G, seed=42)  # You can use other layouts too

    plt.figure(figsize=(6, 6))
    nx.draw(G, pos, with_labels=True, labels={n: n for n in G.nodes()}, 
            node_color=node_colors, node_size=1000, font_color='white', font_weight='bold')
    plt.title("Graph with Colored Nodes (R = red, B = blue)")
    plt.axis('off')
    plt.show()
    
    


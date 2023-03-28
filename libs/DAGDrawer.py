from itertools import count

# Third-party libraries
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

class DagDrawer:
    def __init__(self,  out_path: str = "",
                        title: str = "",
                        VERTICAL_SPACE: int = 14,
                        HORIZONTAL_SPACE: int = 14,
                        X_MARGIN: int = 5,
                        Y_MARGIN: int = 0) -> None:
        self.out_path = out_path if out_path != "" else ("./" + title)
        self.title = title
        self.VERTICAL_SPACE = VERTICAL_SPACE
        self.HORIZONTAL_SPACE = HORIZONTAL_SPACE
        self.X_MARGIN = X_MARGIN
        self.Y_MARGIN = Y_MARGIN

    def arrange_pos_nodes(self, G):
        pos_nodes = {}
        
        topology = [e for e in nx.topological_generations(G)]

        num_generations = len(topology)
        x_offset = self.HORIZONTAL_SPACE / (num_generations - 1)

        for idx_level, level in enumerate(topology):
            num_nodes_in_level = len(level)
            y_offset = self.VERTICAL_SPACE / (num_nodes_in_level + 1)
            for idx_node, node in enumerate(level):
                noise = 0
                if len(level) == 1: # Else labels will be overlapping
                    noise = np.random.uniform(-0.1, 0.1)
                y = self.Y_MARGIN + ((idx_node + 1) * y_offset) + noise
                x = self.X_MARGIN + ((idx_level + 1) * x_offset) + noise
                pos_tuple = (x, y)
                pos_nodes[node] = pos_tuple
        
        return pos_nodes
            
    def construct_attr_dict(self, G):
        attr_nodes = {}
        topology = [e for e in nx.topological_generations(G)]
        for idx_level, level in enumerate(topology):
            for idx_node, node in enumerate(level):
                node_class = None
                if (len(nx.ancestors(G, node)) == 0):
                    node_class = "Input"
                elif (len(nx.descendants(G, node)) == 0):
                    node_class = "Output"
                else:
                    node_class = "Step"

                attr_nodes[node] = {"class": node_class }
        
        return attr_nodes

    def arrange_pos_labels(self, G, pos_nodes):
        pos_labels = { k: (x, y) for k, (x, y) in pos_nodes.items() }

        for node, (x, y) in pos_nodes.items():
            dice = np.random.uniform(0, 1)
            if dice > 0.5:
                pos_labels[node] = (x, y + 0.25)
            else:
                pos_labels[node] = (x, y - 0.25)
            # if y > self.VERTICAL_SPACE / 2:
            #     pos_labels[node] = (x, y + 0.2)
            # else:
            #     pos_labels[node] = (x, y - 0.2)
            
        return pos_labels


    def draw(self, G):
        plt.figure(figsize=(self.HORIZONTAL_SPACE, 
                            self.VERTICAL_SPACE)) 

        # Construct position of each node
        pos_nodes = self.arrange_pos_nodes(G)    

        # Construct position of each label
        pos_labels = self.arrange_pos_labels(G, pos_nodes)

        # Construct attribute of each node
        attr_nodes = self.construct_attr_dict(G)
        nx.set_node_attributes(G, attr_nodes)

        # Construct color map
        groups = set(nx.get_node_attributes(G,'class').values())
        mapping = dict(zip(sorted(groups), count()))
        nodes = G.nodes()
        colors = [mapping[G.nodes[n]['class']] for n in nodes]
        cmap = plt.cm.jet

        # Draw just the nodes
        nx.draw_networkx_nodes(G,  
                            node_color = colors, 
                            cmap=cmap,
                            pos=pos_nodes
        )

        # Draw node labels
        nx.draw_networkx_labels(G, pos = pos_labels,
                                    font_color='black')
        # Draw edges 
        nx.draw_networkx_edges(G, pos = pos_labels,
                            arrowsize=40,
                            connectionstyle='arc3',
                            edge_color='grey')

        # Create scatter plot with correct color legend for each group
        for k, v in mapping.items():
            plt.scatter([],[], c=[cmap(v/np.max(colors))], label='{}'.format(k))

        plt.title(self.title)
        plt.legend(loc='upper left')
        plt.savefig(self.out_path + '.png')

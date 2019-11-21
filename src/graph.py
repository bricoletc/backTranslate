import networkx as nx
from sampler import DNASample


class GraphPurifier:
    """
    'Purify' here means resolving the graph down to an edgeless graph
    with the greatest possible number of nodes (DNASamples)

    The _purify functions are ordered by how well they (seem to) perform, which
    is also how long they take to run.
    """

    def __init__(self, graph: nx.Graph):
        self.graph = graph
        initial_ccs = self.get_connected_components(self.graph)
        self._disjoint_samples = None
        # self._purify_degreeCentrality()
        self._purify_maxIndependentSet()
        # self._purify_minCut(initial_ccs)
        assert (
            len(self.graph.edges) == 0
        ), "ERROR: the graph should be edgeless after call to ._purify()"

    @property
    def disjoint_samples(self):
        if self._disjoint_samples is None:
            self._disjoint_samples = list(self.graph.nodes)
        return self._disjoint_samples

    @staticmethod
    def get_connected_components(graph: nx.Graph):
        """
        .subgraph() returns an immutable *view* on input graph
        """
        connected_components = list(nx.connected_components(graph))
        return [graph.subgraph(cc) for cc in connected_components]

    def _purify_degreeCentrality(self):
        degrees = nx.degree_centrality(self.graph)
        while self.graph.number_of_edges() > 0:
            max_node = max(degrees, key=degrees.get)
            degrees.pop(max_node)
            self.graph.remove_node(max_node)

    def _purify_maxIndependentSet(self):
        nodes = nx.maximal_independent_set(self.graph)
        self.graph = self.graph.subgraph(nodes).copy()

    def _purify_minCut(self, connected_components):
        for cc in connected_components:
            if len(cc.nodes) == 1:
                continue
            # spanning_tree = nx.minimum_spanning_tree(cc)
            # min_cut = list(nx.minimum_node_cut(spanning_tree))
            min_cut = list(nx.minimum_node_cut(cc, flow_func=shortest_augmenting_path))
            self.graph.remove_nodes_from(min_cut)
            new_ccs = self.get_connected_components(cc)
            self._purify_minCut(new_ccs)

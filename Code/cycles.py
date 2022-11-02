from itertools import combinations
import pandas as pd
import os

# This script implements an efficient modified breadth-first search algorithm
# to find the chordless cycles of a zeolite framework graph. It then reduces
# the chordless cycle set to a ring set (a ring is a chorldess cycle that
# cannot be constructed by an XOR operation on the edge sets of smaller rings
# belonging to the graph).

# Please cite the crystallographic data source if you wish to use the data:
# Ch. Baerlocher and L.B. McCusker
# Database of Zeolite Structures: http://www.iza-structure.org/databases/


class Node:
    """
    A class to represent a graph's node with attributes needed for the cycle
    search in the Graph class.

    ...

    Attributes
    ----------
    name : str or int
        Original name of the node, passed from an adjacency list.
    degree : int
        Node degree, or number of connected neighbors.
    label : int
        An integer label 1 to N, where N is number of nodes in graph.
        Generated with algorithm described in Graph class.
    blocked : int
        An integer indicating extent of blocking.
        Generated with algorithm described in Graph class.
    neighbors : list(names)
        A list of names of nodes that are diectly connected to the node
    """

    def __init__(self, name):
        """Initializes a node object with the given name and default values
        for other attributes.

        Args:
            name (str or int): Original name of the node, passed from an
            adjacency list.
        """
        self.name = name
        self.degree = 0
        self.label = 0
        self.blocked = 0
        self.neighbors = []


class Cycle:
    """
    A class to represent cycles. Currently, does not use Node or Graph classes,
    as their attributes are not necessary for Cycle operations.

    ...

    Attributes
    ----------
    nodes : list(names)
        A list of the names of nodes that are sequentially connceted in a cycle
    edges : list(list(names))
        A list of edges that form the cycle, which themselves are a list of two
        nodes. The nodes are not necessarily seqeuential.

    Methods
    -------
    edges_from_nodes():
        Generates the cycle's edges, if cycle initialized with list of nodes
    xor(other):
        Generates a cycle sum/difference defined as an XOR operation on edges
        of the given cycle and the "other"cycle
    check_ring(valid_rings):
        Checks if a cycle can be composed of smaller cycles from basis set

    """

    def edges_from_nodes(self):
        """Generates the Cycle's edges, if Cycle initialized with nodes, by
        appending consecutive nodes as edges, ordering the pairs by node name.
        This removes directional ambiguity when comparing edges.
        """
        for i, node in enumerate(self.nodes):
            if node > self.nodes[i-1]:
                self.edges.append([self.nodes[i-1], node])
            else:
                self.edges.append([node, self.nodes[i-1]])

    def __init__(self, constructor, type='nodes'):
        """Initializes a Cycle using data passed in the "constructor"
        list. Can be constructed using either a list of nodes or list
        of edges.

        Args:
            constructor (list): A sequence of nodes that form a cycle or a list
            of edges.
            type (str, optional): The constructor type to clarify how the cycle
            is initialized. Defaults to 'nodes'.

        Raises:
            Exception: If an improper type is specified, the cycle cannot be
            initialized.
        """
        self.nodes = []
        self.edges = []
        if type == 'nodes':
            self.nodes = constructor
            self.edges_from_nodes()
        elif type == 'edges':
            self.edges = constructor
        else:
            raise Exception("Cycle can only be constructed from sequential "
                            "lists of nodes or edges ")

    def xor(self, other):
        """Generates a cycle sum/difference defined as an XOR operation on
        edges of the given cycle and the "other" cycle.

        Args:
            other (Cycle): The other Cycle object with which the XOR is
            performed.

        Raises:
            ValueError: Raised if the XOR operation generates two cycles.

        Returns:
            (Cycle): Resulting Cycle object of the XOR operation.
        """
        # Initialize common and exclussive nodes and edges.
        common_edges = []
        common_nodes = []
        self_xor_b = self.edges[:]
        # For each edge in the "other" cycle, check if node is also in self.
        for edge in other.edges:
            if edge in self.edges:
                common_edges.append(edge)
                common_nodes += edge
                self_xor_b.remove(edge)
            else:
                self_xor_b.append(edge)
        # Check if XOR yields a single cycle. Raise error if not
        if len(set(common_nodes)) > len(common_edges)+1:
            raise ValueError("The cycle was split into two cycles")
        return Cycle(self_xor_b, 'edges')

    def check_ring(self, valid_rings):
        """Checks if a cycle can be composed of smaller cycles from basis list
        of rings, which is passed as "valid_rings" arg.

        Args:
            valid_rings (list(Cycle)): A basis list of smaller Cycle objects
            that have been found to be rings (i.e., Cycles that cannot be
            composed through XOR operation of two or more other rings) of the
            graph.

        Returns:
            (bool): True if cycle is a ring.
        """
        # Simplify naming scheme such that c = a XOR b, where a is original
        # cycle, b is a ring from the basis set and c is the XOR result.
        a_edge_count = len(self.edges)
        # Each ring in the valid_rings list is checked.
        for i, b in enumerate(valid_rings):
            try:
                c = self.xor(b)
            except ValueError:
                # The larger cycle was split into two cycles
                return False
            c_edge_count = len(c.edges)
            # Collapse the list of lists to get list of nodes
            c_node_count = len(
                set([node for edge in c.edges for node in edge]))
            if c_edge_count > a_edge_count:
                # The cycle c is larger than cycle a, so a may still be a ring
                continue
            if c_edge_count != c_node_count:
                # Failed because nodes not equal to edges, and c is not a cycle
                return False
            if c_edge_count == a_edge_count:
                if len(valid_rings) == 1:
                    # If only one basis ring was passed, cycle a could not be
                    # constructed by an XOR sum of smaller rings, and is thus a
                    # ring itself
                    return True
                # If additional rings are left to check, recursively call
                # the check_ring() method, passing a reduced basis, excluding
                # the already checked rings.
                ring = c.check_ring(valid_rings[i+1:])
                if not ring:
                    # If further recursion indicates that c was not a ring, a
                    # is also not a ring
                    return False
            elif c_edge_count < a_edge_count:
                # Found decomposition to smaller rings, so a was not a ring
                return False
        return True


class Graph:
    """
    A class to represent a graph to be used for a cycle search.
    Only nodes implemented here because edge operations were not necessary.
    The cycle search algorithm is a slightly modified implementation of
    the chordless cycle search algorithm described in:
    https://doi.org/10.48550/arXiv.1309.1051
    This is a modified breadth-first search with an O(M+N) time complexity,
    where M is the number of edges, and N is the number of nodes. For zeolites,
    the Si-O-Si connectivity leads to approximately 4 edges/node (excluding
    edges of crystal structure and under-coordinated sites), leading to
    approximately O(5N) complexity, which reduces to O(N).

    Though much of the pseudocode is directly implemented, some modifications
    were made to make the code more Pythonic. An option for a cutoff cycle
    length is added to terminate the search of any cycles longer than a cutoff.
    This improves performance, if the maximum cycle size of interest is known.
    ...

    Attributes
    ----------
    nodes : dictionary
        A dictionary that stores user-defined adjacency list.
    triplets : list
        A list of triplets that serve as origin of cycle search.
    cycles : list
        A list of chordless cycles found by the cycle search.

    Methods
    -------
    degree_labeling():
        Generates node labels based on node degrees (see reference for details)
    triplets_search():
        Finds valid triplets of nodes that are the starting sequences for the
        cycle search method, cc_visit()
    block_neighbors(v):
        Increments the "blocked" label of neighbors of node v.
    unblock_neighbors(v):
        Decrements the "blocked" label of neighbors of node v.
    cc_visit(p, cutoff):
        Recursively elongates node sequence until cycle found or cutoff size is
        reached
    chorldess_cycles(cutoff=None):
        Finds all chordless cycles of the graph that are smaller in size than
        specified cutoff. The default cutoff is None, resulting in no
        restrictions on cycle size.
    """

    def __init__(self, adjacency):
        """Initializes a Graph object given an adjacency dictionary.

        Args:
            adjacency (dict): A user-provided adjacency dictionary with nodes
            for keys, and a list of the node's neighbors for values
        """
        # Initialize blank graph attributes
        self.nodes = {}
        self.triplets = []
        self.cycles = []

        # Convert the user-provided adjacency dictionary to Graph object.
        for node, node_neighbors in adjacency.items():
            # Use only unique neighbors
            node_neighbors = list(set(node_neighbors))
            if node not in self.nodes:
                node_i = Node(node)  # initiate a Node object
                # Get the neighbors directly indicated in adjacency matrix
                node_i.neighbors = node_neighbors
                self.nodes[node] = node_i  # assign node to graph object
            else:
                self.nodes[node].neighbors = list(
                    set(self.nodes[node].neighbors + node_neighbors))

            # Initiate or update the nodes corresponding to neighbors
            # This makes sure that nodes that are pointed to, but are not
            # keyed in the adjacency dict are included
            for neighbor in node_neighbors:
                if neighbor not in self.nodes:
                    neighbor_node = Node(neighbor)
                    neighbor_node.neighbors = [node]
                    self.nodes[neighbor] = neighbor_node
                else:
                    if node not in self.nodes[neighbor].neighbors:
                        self.nodes[neighbor].neighbors.append(node)

        # Calculate the degree for each node
        for node in self.nodes:
            self.nodes[node].degree = len(self.nodes[node].neighbors)

    def degree_labeling(self):
        """Generates node labels based on node degrees (see
         https://doi.org/10.48550/arXiv.1309.1051 for details)
        """
        # Make a temp degree dict to keep all unlabeled nodes
        degrees = {node: self.nodes[node].degree for node in self.nodes}
        for i in range(len(self.nodes)):
            # Set v to an unlabeled node with minimum degree
            v = min(degrees, key=degrees.get)
            # Label the v node with integer "i" and pop v after it is labeled
            self.nodes[v].label = i
            degrees.pop(v)
            # Decrement the degree of each uncolored neighbor of v
            for u in self.nodes[v].neighbors:
                if u in degrees:
                    degrees[u] -= 1

    def triplets_search(self):
        """Finds valid triplets of nodes that are the starting sequences for
        cycle search method, cc_visit(). Triplets of form (x,u,y) must satisfy
        the condition that: the y.label > x.label > u.label.
        """
        for u in self.nodes:
            neighbors = [
                i for i in self.nodes[u].neighbors if
                self.nodes[i].label > self.nodes[u].label]
            # Form triplet from all pair combinations of neighbors
            for x, y in list(combinations(neighbors, 2)):
                # Ensure that y label is > x label
                if self.nodes[x].label > self.nodes[y].label:
                    x, y = y, x
                # If x and y are connected, found a 3-ring
                if y in self.nodes[x].neighbors:
                    self.cycles.append([x, u, y])
                # Else, a feasible path origin triplet was found
                else:
                    self.triplets.append([x, u, y])

    def block_neighbors(self, v):
        """Increments the "blocked" label of neighbors of node v.

        Args:
            v (int or str): A label of node whose neighbors must be blocked
        """
        for u in self.nodes[v].neighbors:
            self.nodes[u].blocked += 1

    def unblock_neighbors(self, v):
        """Decrements the "blocked" label of neighbors of node v.

        Args:
            v (int or str): A label of node whose neighbors must be unblocked
        """
        for u in self.nodes[v].neighbors:
            if self.nodes[u].blocked > 0:
                self.nodes[u].blocked -= 1

    def cc_visit(self, p, cutoff):
        """Recursively elongates node sequence until cycle found or cutoff
        size is reached

        Args:
            p (list): A sequence of nodes that forms a path that could become
            a cycle
            cutoff (int): A limit (inclusive) of cycle size to consider.
        """
        # Early stop if path length is above a cutoff for ring size
        if len(p) > cutoff-1:
            return

        # Block the neighbors of the last element in path
        self.block_neighbors(p[-1])

        # Try to elongate path with the neighbors of the last element in path
        for v in self.nodes[p[-1]].neighbors:
            # If a neighbor's label is greater than the path's 2nd element, and
            # its blocked_counter is 1, append it to path
            if (self.nodes[v].label > self.nodes[p[1]].label
                    and self.nodes[v].blocked == 1):
                p_temp = p + [v]
                # If a cycle is formed, add to cycles list
                if p[0] in self.nodes[v].neighbors:
                    self.cycles.append(p_temp)
                # Otherwise, call cc_visit, recursively
                else:
                    self.cc_visit(p_temp, cutoff)

        # Unblock the neighbors that were blocked in this function call
        self.unblock_neighbors(p[-1])

    def chordless_cycles(self, cutoff=None):
        """Finds all chordless cycles of the graph that are smaller in size
        than specified cutoff. The default cutoff is None, resulting in no
        restrictions on cycle size.

        Args:
            cutoff (int, optional): A limit (inclusive) of cycle size to
            consider. Defaults to None, in which case all possible cycles
            are found.
        """
        # If no cutoff supplied, explore cycles sizes up to number of nodes
        if cutoff is None:
            cutoff = len(self.nodes)
        # Use algorithm-specific node labeling to initialize cycle search
        self.degree_labeling()
        # Initialize the searchable path origins and find all 3-membered rings
        self.triplets_search()
        # Use every triplet origin only once to search for cycles
        while len(self.triplets) > 0:
            # Remove an element from triplets array and search
            p = self.triplets.pop()
            # Increment the block counter for neighbors of the triplet center
            self.block_neighbors(p[1])
            # Visit unblocked neighbors of the last node and check for a cycle
            self.cc_visit(p, cutoff)
            # Decrement the block counter for neighbors of the triplet center
            self.unblock_neighbors(p[1])


def df_to_graph(positions_df):
    """
    Turns a pandas dataframe to a graph object. Specifically implemented for
    Si-O connectivity in mind, with parsing logic for externally-defined data
    formatting. Easily extended to other graph adjacency types, once data-
    specific parsing is known.

    Args:
        positions_df (pandas.Dataframe): A pandas dataframe that contains all
        Si and O connectivity.

    Returns:
        (Graph): A graph object that contains methods to perform a search for
        cordless cycles.
    """
    # Clean up dataframe and get index of oxygens
    positions_df.fillna('', inplace=True)
    o_index = positions_df['Atom'] == 'O'
    graph_dict = {}
    # For all oxygens with more than one bonded silicon atom, append the graph
    # adjacency dictionary
    for o in o_index[o_index].index:
        if positions_df.loc[o, 'Bonded Atoms']:
            si_bond = positions_df.loc[o, 'Bonded Atoms'].split(';')[:-1]
            si_bond = [int(x) for x in si_bond]
            if len(si_bond) > 1:
                if si_bond[0] not in graph_dict:
                    graph_dict[si_bond[0]] = [si_bond[1]]
                else:
                    graph_dict[si_bond[0]] += [si_bond[1]]
    return Graph(graph_dict)


# Process all files generated by bond_search.py script to find all rings
# with 6 or fewer nodes.
path = "./FrameworkBonds"
dir_list = os.listdir(path)
for file in dir_list:
    positions_df = pd.read_csv(f"./FrameworkBonds/{file}")
    # Transform the dataframe to a graph object
    zeolite_graph = df_to_graph(positions_df)
    # Search for chordless cycles, using the largest ring size of 6 as a cutoff
    zeolite_graph.chordless_cycles(6)

    # Create a dictionary to hold possible rings and populate with candidate
    # cycles generated by the search for chorldess cycles
    ring_sizes = [3, 4, 5, 6]
    possible_rings = dict([(i, []) for i in ring_sizes])
    for A in zeolite_graph.cycles:
        if len(A) in ring_sizes:
            possible_rings[len(A)].append(A)

    # Because 3, 4, and 5 membered chorldess cycles are rings, and 6 membered
    # chorldess cycles cannot be constructed from 3 or 5 membered rings, only
    # need to check if the discovered 6 membered chorldess cycles can be
    # constructed from the 4 membered rings to see if they are 6 membered rings
    i_pop = []
    for i, A in enumerate(possible_rings[6]):
        A = Cycle(A)
        ring = A.check_ring([Cycle(ring) for ring in possible_rings[4]])
        if not ring:
            i_pop.append(i)
    # Remove all 6-membered cycles that are not rings
    for i in i_pop[::-1]:
        possible_rings[6].pop(i)

    # Convert indeces to positions
    rings_xyz = []
    for size in possible_rings.keys():
        for ring in possible_rings[size]:
            # Interlace oxygens between Si atoms
            o_neighbors = positions_df.loc[ring, 'Bonded Atoms']
            interlace = []
            for i, si in enumerate(ring):
                o_1 = o_neighbors.iloc[i-1].split(';')
                o_2 = o_neighbors.iloc[i].split(';')
                o_common = set(o_1).intersection(o_2).pop()
                interlace.append(int(o_common))
                interlace.append(si)
            xyzs = positions_df.loc[interlace, [
                'x', 'y', 'z']].values.flatten().tolist()
            rings_xyz.append(xyzs)

    # Save the ring atom coordinates as a series of flat lists of form:
    # [x0,y0,z0,...xN,yN,zN]. These coordinates can be used to visualize
    # the rings.
    rings_series = pd.Series(rings_xyz)
    rings_series.to_csv(
        f"./Rings/{file[:-4]}_rings.csv", index=False)

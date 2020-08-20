"""Methods to create, filter, and export moelcular network
"""

import networkx as nx

def make_network(nodes,edges):
    """Takes a list of Spectrum objects and a dictionary of spectrum matches and returns a NetworkX graph object
    """
    #make network
    network=nx.Graph(edges)

    #add metadata to nodes
    attributes={}
    for N in nodes:
        network.add_node(N)
        attributes[N]={}
        if hasattr(N,'parameters'):
            for P in N.parameters:
                attributes[N][P]=str(N.parameters[P])
        
        if hasattr(N, 'library_parameters'):
            network.add_node(N,library_match=True)
            for L in N.library_parameters:
                attributes[N]["".join(["library_",L])]=N.library_parameters[L]
        else:
            network.add_node(N,library_match=False)
          
    nx.set_node_attributes(network,attributes)

    return network

def filter_neighbors(graph,M):
    """Takes networkx graph object and an int as the maximum numbers of connections a node can have
    Edges are removed from each node from the lowest cosine score, until the node has below the threshold number.
    Returns the filtered network
    """
    for node in nx.nodes(graph):

        #sort edges connecting to query node by descending cosine score
        edges = sorted(graph.edges(node,data=True), key=lambda x: x[2]['cosine'],reverse=True)

        if len(edges)>M: #if too many neighbours i.e. connecting edges
            graph.remove_edges_from(edges[M:]) #remove edges from lowest cosine score end of list
    
    return graph

def filter_family(graph, M):
    """Takes a networkx graph object and an int for threshold molecular family size in the network.
    Removes edges from families that are above the threshold, from the lowest cosine score edges, until the family is small enough.
    Returns the filtered network.
    """

    used=[]
    for node in nx.nodes(graph):

        #don't recheck any nodes in a family that has already been filtered
        if node in used:
            continue

        while True:
            family=nx.node_connected_component(graph,node)

            if len(family)>M:
                #sort from highest to lowest cosine
                edges = sorted(graph.edges(family,data=True), key=lambda x: x[2]['cosine'],reverse=True)
                while True:
                    #get edge to remove
                    remove=edges[-1]
                    #remove that edge from the lsit of existing edges
                    edges=edges[:-1]
                    #remove edge from the network
                    graph.remove_edge(remove[0],remove[1])

                    #check if there is a path from the two nodes that the edge was removed from
                    #if a path exists then the family was not broken and another edge should be removed
                    if not nx.has_path(graph,remove[0],remove[1]):
                        break
                
                continue
                
            for F in family:
                used.append(F)
            break

    return graph

def write_graphml(graph, file_name):
    """Takes a networkx graph object and a file name and writes the network to a graphml file
    """
    nx.write_graphml(graph, file_name)
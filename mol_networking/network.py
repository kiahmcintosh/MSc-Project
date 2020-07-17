import networkx as nx

def make_network(nodes,edges):

    network=nx.Graph(edges)
    attributes={}
    for N in nodes:
        network.add_node(N)
        attributes[N]={}
        if hasattr(N,'parameters'):
            for P in N.parameters:
                attributes[N][P]=N.parameters[P]
        
        if hasattr(N, 'library_parameters'):
            network.add_node(N,library_match=True)
            for L in N.library_parameters:
                attributes[N]["".join(["library_",L])]=N.library_parameters[L]
        else:
            network.add_node(N,library_match=False)
          
    nx.set_node_attributes(network,attributes)

    return network

def filter_neighbors(graph,M):
    for node in nx.nodes(graph):

        #sort edges connecting to query node by descending cosine score
        edges = sorted(graph.edges(node,data=True), key=lambda x: x[2]['cosine'],reverse=True)

        if len(edges)>M: #if too many neighbours i.e. connecting edges
            graph.remove_edges_from(edges[M:]) #remove edges from lowest cosine score end of list
    
    return graph

def filter_family(graph, M):
    used=[]
    for node in nx.nodes(graph):
        if node in used:
            # print(f"\nnode {node} already used")
            continue

        while True:
            family=nx.node_connected_component(graph,node)

            if len(family)>M:
                edges = sorted(graph.edges(family,data=True), key=lambda x: x[2]['cosine'],reverse=True)
                # print(edges)
                while True:
                    remove=edges[-1]
                    edges=edges[:-1]
                    graph.remove_edge(remove[0],remove[1])
                    
                    if not nx.has_path(graph,remove[0],remove[1]):
                        break
                
                continue
                
            for F in family:
                used.append(F)
            break

    return graph

def write_graphml(graph, file_name):
    nx.write_graphml(graph, file_name)
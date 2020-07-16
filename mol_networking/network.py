import networkx as nx

def make_network(nodes,edges):

    graph=nx.Graph(edges)
    for N in nodes:
        if hasattr(N, 'name'):
            graph.add_node(N,colour=1,ID=N.name)
        else:
            graph.add_node(N,colour=0)

    return graph

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
            # print("\n",family)
            # print(f"\n{node}\t{family}\t{len(family)}")

            if len(family)>M:
                edges = sorted(graph.edges(family,data=True), key=lambda x: x[2]['cosine'],reverse=True)
                # print(edges)
                while True:
                    remove=edges[-1]
                    edges=edges[:-1]
                    # print("remove: ",remove)
                    # family.remove(remove[1])
                    # print(family)
                    graph.remove_edge(remove[0],remove[1])
                    
                    # print(remove[0],"\t",remove[1],"\t",nx.has_path(graph,remove[0],remove[1]))
                    if not nx.has_path(graph,remove[0],remove[1]):
                        break
                
                continue
                
            for F in family:
                used.append(F)
            break

    return graph

def write_graphml(graph, file_name):
    nx.write_graphml(graph, file_name)
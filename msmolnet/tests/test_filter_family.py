from msmolnet import network
import networkx as nx

def test_filter_family():
    nodes=[1,2,3,4,5,6,7]
    edges=[(1,2,{'cosine':0.9}),(1,3,{'cosine':0.5}),(1,4,{'cosine':0.2}),(2,5,{'cosine':0.4}),(5,6,{'cosine':0.6})]
    test_network=network.make_network(nodes,edges)

    filtered_network=network.filter_family(test_network,3)

    family=list(nx.node_connected_component(filtered_network,2))

    assert family == [1,2,3], "Incorrect family"


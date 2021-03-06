"""Method to test n neighbours filtering
"""

from msmolnet import network
import networkx as nx

def test_filter_neighbors():
    """Tests network n neighbours filtering method.
    Throws an assertion error if the filtering is not correct
    """
    nodes=[1,2,3,4,5,6,7,8]
    edges=[(1,2,{'cosine':0.5}),(1,3,{'cosine':0.6}),(1,4,{'cosine':0.9}),(1,5,{'cosine':0.75}),(1,6,{'cosine':0.2}),(1,7,{'cosine':0.5}),]
    test_network=network.make_network(nodes,edges)

    filtered_network=network.filter_neighbors(test_network,3)
    neighbors=list(filtered_network.edges(1))

    assert len(neighbors) == 3, "Incorrect number of neighbors"
    assert neighbors == [(1,3),(1,4),(1,5)], "Expected different neighbors"
    
   
import time
import random
from pathlib import Path
import networkx as nx
from graphs import bidir_dijkstra


def fetch_graph(file_name):
    with open(file_name) as input_file:
        line = input_file.readline().rstrip('\n')
        n_nodes, n_vertices = map(int, line.split(' '))
        edges = []
        for _ in range(n_vertices):
            line = input_file.readline().rstrip('\n')
            vertex1, vertex2, weight = map(int, line.split(' '))
            edges.append(((vertex1, vertex2, weight)))
        line = input_file.readline().rstrip('\n')
        n_queries = int(line)
        queries = []
        for _ in range(n_queries):
            line = input_file.readline().rstrip('\n')
            vertex1, vertex2 = map(int, line.split(' '))
            queries.append(((vertex1, vertex2)))

    return edges, queries


def fetch_social_media_combined(file_name):
    with open(file_name) as input_file:
        line = input_file.readline().rstrip('\n')
        edges = []
        while line:
            v1, v2 = map(int, line.split(' '))
            edges.append((v1, v2))
            line = input_file.readline().rstrip('\n')

    graph = nx.DiGraph()
    graph.add_edges_from([(v1, v2, {'weight': 1}) for (v1, v2) in edges])

    return graph


def build_graph(edges) -> nx.DiGraph:
    graph = nx.DiGraph()
    weighted_edges = [(v1, v2, {'weight': w}) for (v1, v2, w) in edges]
    graph.add_edges_from(weighted_edges)
    return graph


def test_bidir_dijkstra():
    edges, queries = fetch_graph(Path('../test/test03.txt'))
    graph = build_graph(edges)
    for source, sink in queries:
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == 5
        assert path == [0, 1, 4, 5]

    edges, queries = fetch_graph(Path('../test/test04.txt'))
    graph = build_graph(edges)
    for source, sink in queries:
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == 5
        assert path == [0, 1, 3, 4]

    edges, queries = fetch_graph(Path('../test/test05.txt'))
    graph = build_graph(edges)
    for source, sink in queries:
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == 1
        assert path == [0, 1]

    edges, queries = fetch_graph(Path('../test/test06.txt'))
    graph = build_graph(edges)
    for source, sink in queries:
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == 3
        assert path == [0, 1, 2]

    edges, queries = fetch_graph(Path('../test/test01.txt'))
    graph = build_graph(edges)
    for (source, sink), expected in zip(queries, [0, 0, 1, -1]):
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == expected

    edges, queries = fetch_graph(Path('../test/test02.txt'))
    graph = build_graph(edges)
    for source, sink in queries:
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == 3
        assert path == [1, 2, 3]

    graph = fetch_social_media_combined(Path('../test/facebook_combined.txt'))
    '''    try:
        expected = len(nx.shortest_path(graph, source=449, target=626)) - 1
    except nx.NetworkXNoPath:
        expected = -1
    distance, path = bidir_dijkstra(graph, 449, 626)
    assert distance == expected

    return '''

    """
    random.seed(42)
    for i in range(100):
        source = random.randint(0, 4031)
        sink = random.randint(0, 4031)
        # print(i, source, sink)
        try:
            expected = len(nx.shortest_path(graph, source=source, target=sink)) - 1
        except nx.NetworkXNoPath:
            expected = -1
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == expected
    """

    graph = fetch_social_media_combined(Path('../test/twitter_combined.txt'))
    nodes = list(graph.nodes)
    random.seed(42)
    print()
    for i in range(10):
        source = random.choice(nodes)
        sink = random.choice(nodes)
        print(i, source, sink)
        t1 = time.process_time()
        try:
            expected = len(nx.shortest_path(graph, source=source, target=sink)) - 1
        except nx.NetworkXNoPath:
            expected = -1
        t2 = time.process_time()
        distance, path = bidir_dijkstra(graph, source, sink)
        t3 = time.process_time()
        assert distance == expected
        print(distance, 'nx=',t2-t1, 'bidir-dijkstra=',t3-t2)

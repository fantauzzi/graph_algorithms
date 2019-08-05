import os
import urllib.request
import gzip
import random
from pathlib import Path
import networkx as nx
from graphs import bidir_dijkstra


def fetch_gzip(file_name, url):
    if not os.path.exists(Path(file_name)):
        if not os.path.exists(Path(file_name + '.gz')):
            print('\nDownloading file', Path(file_name + '.gz'))
            urllib.request.urlretrieve(url, Path(file_name + '.gz'))
        print('Unzipping file', Path(file_name + '.gz'))
        with gzip.open(Path(file_name + '.gz'), 'rb') as gzip_file, open(Path(file_name), 'wb') as out_file:
            out_file.writelines(gzip_file)
        os.remove(Path(file_name + '.gz'))


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
    graph.add_edges_from([(v2, v1, {'weight': 1}) for (v1, v2) in edges])

    return graph


def build_graph(edges) -> nx.DiGraph:
    graph = nx.DiGraph()
    weighted_edges = [(v1, v2, {'weight': w}) for (v1, v2, w) in edges]
    graph.add_edges_from(weighted_edges)
    return graph


def test_bidir_dijkstra():
    '''
    Downloads datasets from "SNAP Datasets: Stanford Large Network Dataset Collection
    Jure Leskovec and Andrej Krevl
    http://snap.stanford.edu/data
    June, 2014
    '''

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

    fb_filename = '../test/facebook_combined.txt'
    fetch_gzip(fb_filename, 'https://snap.stanford.edu/data/facebook_combined.txt.gz')
    graph = fetch_social_media_combined(Path('../test/facebook_combined.txt'))

    random.seed(42)
    for i in range(1):
        source = random.randint(0, 4031)
        sink = random.randint(0, 4031)
        try:
            expected = len(nx.shortest_path(graph, source=source, target=sink)) - 1
        except nx.NetworkXNoPath:
            expected = -1
        distance, path = bidir_dijkstra(graph, source, sink)
        assert distance == expected

    twitter_filename = '../test/twitter_combined.txt'
    fetch_gzip(twitter_filename, 'https://snap.stanford.edu/data/twitter_combined.txt.gz')
    graph = fetch_social_media_combined(twitter_filename)
    nodes = list(graph.nodes)

    random.seed(42)
    with open(Path('../test/twitter_tcs.txt'), 'wt') as output_file:
        for i in range(10):
            source = random.choice(nodes)
            sink = random.choice(nodes)
            try:
                expected = len(nx.shortest_path(graph, source=source, target=sink)) - 1
            except nx.NetworkXNoPath:
                expected = -1
            distance, path = bidir_dijkstra(graph, source, sink)
            assert distance == expected
            print(source, sink, expected, file = output_file)

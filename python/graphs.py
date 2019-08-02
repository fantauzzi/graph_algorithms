import networkx as nx
from sortedcontainers import SortedSet


def bidir_dijkstra(graph, source, sink):
    """
    :param graph:
    :param source:
    :param sink:
    :return:
    """

    if source == sink:
        return 0, []

    ''' Use two graphs: graph_1 is the same as the one passed to the function, graph_2 is obtained from graph_1
    by flipping the orientation of every edge, and maintaining the same distances (weights). '''
    graph_1 = graph
    graph_2 = nx.DiGraph()
    graph_2.add_edges_from([(v2, v1) for (v1, v2) in graph_1.edges])
    for v1, v2 in graph_1.edges:
        graph_2.edges[v2, v1]['weight'] = graph_1.edges[v1, v2]['weight']

    ''' Set of vertices yet to be processed, each with the current estimate of its shortest distance from the
    vertex of its graph. Note that graph_2 source is graph_1 sink, and graph_2 sink is graph_1 source.'''
    pending_1 = SortedSet([(0, source)])
    pending_2 = SortedSet([(0, sink)])

    ''' Associate each vertex with the current estimate of its shortest distance from the source in its graph.
    This information is redundant, as it is also contained in pending_1 and pending_2; however, having it in maps
    like these allows to look it up in constant time. Consequently, d_1 must be kept consistent with pending_1,
    and d_2 with pending_2: when updating the estimate for a given vertex, it must be updated in both the SortedSet
    and the map.'''
    d_1 = {source: 0}
    d_2 = {sink: 0}

    ''' Set of vertices already processed. Which implies, its estimate of shortest distance from the source is
    the exact measurement.'''
    processed_1 = set()
    processed_2 = set()

    ''' Will store information useful to backtrack the shortest path to source/sink.'''
    pred_1 = {source: None}
    pred_2 = {sink: None}

    '''for vertex2 in graph_1.adj[sink]:
        d_2[vertex2] = float('inf')
        pending_2.add((float('inf'), sink))'''

    shortest_so_far = float('inf')
    shortest_path_found = False
    step = 0
    vertex1, vertex2, best_vertex2 = None, None, None  # Useful for debugging. At the end of the while loop, both must be not None
    while not shortest_path_found:
        # Do one step of Dijkstra in graph_1
        if not pending_1:
            return -1, []
        d_f1, vertex1 = pending_1.pop(0)  # Pop the minimum.
        assert d_f1 == d_1[vertex1]  # Check the invariant (pending_1 and d_1 must be consistent)
        processed_1.add(vertex1)
        # Add vertices adjacent to vertex1, that have not been processed, to the pending vertices (unless already there)
        for vertex2 in graph_1.adj[vertex1]:
            if d_1.get(vertex2) is None:
                assert vertex2 not in processed_1
                d_1[vertex2] = float('inf')
                pending_1.add((float('inf'), vertex2))
        best_vertex2 = None
        for vertex2 in graph_1.adj[vertex1]:
            # Length of the shortest path from the source to vertex2, going through vertex1.
            d_via_vertex1 = d_f1 + graph_1.edges[vertex1, vertex2]['weight']
            # Relax edge (vertex1, vertex2)
            if vertex2 not in processed_1:
                d_f2 = d_1[vertex2]
                if d_via_vertex1 < d_f2:
                    d_1[vertex2] = d_via_vertex1
                    pending_1.remove((d_f2, vertex2))
                    pending_1.add((d_via_vertex1, vertex2))
                    pred_1[vertex2] = vertex1
            ''' If vertex2 has already been processed backward, then you have found a path
            from source to sink; update information on the shortest path found so far as necessary. '''
            if vertex2 in processed_2:
                length = d_via_vertex1 + d_2[vertex2]
                shortest_so_far = min(shortest_so_far, length)
                best_vertex2 = vertex2
        # Check termination condition.
        if pending_1 and pending_2:
            l1, _ = pending_1[0]
            l2, _ = pending_2[0]
            if l1 + l2 >= shortest_so_far:
                shortest_path_found = True
        ''' Trade graph_1 with grap_2 before the next step of Dijkstra. Therefore alternate one step in graph_1
        with one step in graph_2. '''
        graph_1, graph_2 = graph_2, graph_1
        pending_1, pending_2 = pending_2, pending_1
        d_1, d_2 = d_2, d_1
        processed_1, processed_2 = processed_2, processed_1
        pred_1, pred_2 = pred_2, pred_1
        step += 1

    # A shortest path from source to sink was found, and it must go through vertex1 and then best_vertex2
    assert vertex1 is not None
    assert best_vertex2 is not None

    # Ensure that graph_1 is the same as graph, swapping graph_1 and graph_2 one more time if necessary.
    if step % 2 == 1:
        graph_1, graph_2 = graph_2, graph_1
        pending_1, pending_2 = pending_2, pending_1
        d_1, d_2 = d_2, d_1
        processed_1, processed_2 = processed_2, processed_1
        pred_1, pred_2 = pred_2, pred_1
    else:
        vertex1, best_vertex2 = best_vertex2, vertex1

    ''' Stitch together the shortest path as 
    source -> ... -> vertex1 -> best_vertex2 -> ... -> sink'''

    def backtrack_path(from_vertex, pred):
        path = [from_vertex]
        previous = pred[from_vertex]
        while previous is not None:
            path.append(previous)
            previous = pred[previous]
        return path

    sub_path_1 = backtrack_path(vertex1, pred_1)  # This is reversed
    sub_path_2 = backtrack_path(best_vertex2, pred_2)
    path = sub_path_1[::-1] + sub_path_2

    # Compute the length of the shortest path, as the sum of the lengths of the stiched sub-paths.
    distance = d_1[best_vertex2] + d_2[best_vertex2] if step % 2 == 1 else d_1[vertex1] + d_2[vertex1]

    return distance, path

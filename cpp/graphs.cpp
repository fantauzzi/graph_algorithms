#include <utility>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <limits>

#include "../include/catch.hpp"

#include "graphs.h"

using std::string;
using std::vector;
using std::pair;
using std::tuple;
using std::make_tuple;
using std::make_pair;
using std::endl;
using std::unordered_map;
using std::map;
using std::multimap;
using std::unordered_set;
using std::set;
using std::swap;
using std::numeric_limits;

typedef vector<pair<int, int>> AdjList;
typedef unordered_map<int, AdjList> Graph;

void append_adj(Graph &graph, int v1, int v2, int w)
/**
 * Add an edge to a directed graph, with a given weight. It is possible to have multiple edges, with the same
 * orientation, between the same two vertices.
 * @param graph the graph.
 * @param v1 the vertex from where the edge originates.
 * @param v2 the vertex where the edge goes.
 * @param w the weight for the edge.
 */
{
    if (graph.find(v1) == graph.end())
        graph[v1] = {{v2, w}};
    else
        graph[v1].emplace_back(make_pair(v2, w));
}


template<typename Container, typename Value>
bool in(Value value, const Container &container)
/**
 * Verifies if a given values is in a given container.
 * @tparam Container the container type. It must implement method find().
 * @tparam Value the value type.
 * @param value the given value.
 * @param container the given container.
 * @return true if the given value is found in the container, false otherwise.
 */
{
    return container.find(value) != container.end();
}


const auto None = numeric_limits<int>::min();

vector<int> backtrack_path(int from_vertex, const map<int, int> &pred) {
    vector<int> path = {from_vertex};
    auto previous = pred.at(from_vertex);
    while (previous != None) {
        path.emplace_back(previous);
        previous = pred.at(previous);
    }

    return path;
}


pair<int, vector<int>> bidir_dijkstra(const Graph &graph, int source, int sink)
/**
 * Returns the length of the shortest path, and the shortest path, in a weighted, directed graph. They are computed
 * using bi-directional Dijkstra. Weights are assumed to be non-negative.
 * @param graph the directed, weighted graph.
 * @param source the source vertex in the graph, where the shortest path must begin.
 * @param sink the sink (destination) vertex in the graph, where the shortest path must end.
 * @return a pair with the length of and the sequence of vertices along it. If there are multiple shortest paths,
 * then one of them is returned in the vector. If sink is not reachable from source, that is there isn't any path in
 * the graph from source to sink, then returns the pair -1, {}.
 */
{
    if (source == sink)
        return make_pair<int, vector<int>>(0, {});

    Graph const *graph_1 = &graph;
    Graph reverse_graph;
    for (const auto &item: graph) {
        auto v1 = item.first;
        auto adj_list = item.second;
        for (const auto &adj: adj_list) {
            auto v2 = adj.first;
            auto w = adj.second;
            append_adj(reverse_graph, v2, v1, w);
        }
    }
    Graph const *graph_2 = &reverse_graph;

    /* Vertices yet to be processed, for either direction. Each vertex is stored with the current estimate of its
     * shortest distance from its source, as a pair <distance, vertex>. It allows accessing the vertex with the lowest
     * distance in constant time. */
    set<pair<int, int>> pending_1 = {{0, source}};
    set<pair<int, int>> pending_2 = {{0, sink}};

    /* Associates each vertex with the current estimate of its shortest distance from its source, the reciprocal of
     * pending_1 and pending_2. It allows finding the current distance estimate of a given vertex in constant time. */
    map<int, int> d_1 = {{source, 0}};
    map<int, int> d_2 = {{sink, 0}};

    /* Vertices already processed for the respective directions. If a vertex is in here, then its current distance
     * estimate from its source is the actual shortest distance */
    unordered_set<int> processed_1;
    unordered_set<int> processed_2;

    // Information useful to backtrack the shortest path to source/sink
    map<int, int> pred_1 = {{source, None}};
    map<int, int> pred_2 = {{sink, None}};

    bool shortest_path_found = false;
    int shortest_so_far = numeric_limits<int>::max();
    int counter_at_shortest_path = None;
    int best_vertex1 = None;
    int best_vertex2 = None;
    int step_counter = -1;
    while (!shortest_path_found) {
        ++step_counter;
        // One step of Dijkstra, along on of the two directions
        if (pending_1.empty()) // No more pending vertices => sink is not reachable from source
            return make_pair<int, vector<int>>(-1, {});
        /*  Pop the vertex with current minimum estimate of its shortest distance from its source, and mark it as
         * processed */
        auto[df_1, vertex1] = *pending_1.cbegin();
        assert(df_1 == d_1[vertex1]);  // Invariant: pending_1 and d_1 must be consistent
        pending_1.extract(pending_1.begin());
        processed_1.emplace(vertex1);
        // Add un-processed vertices adjacent to vertex1 to the pending vertices (unless already there).
        if (in(vertex1, *graph_1)) { // if vertex1 is not in the adj. list *graph1, then it has no outgoing edges
            for (const auto[vertex2, weight]: graph_1->at(vertex1))
                if (!in(vertex2, d_1)) {
                    assert(!in(vertex2, processed_1));
                    d_1[vertex2] = numeric_limits<int>::max();
                    pending_1.emplace(make_pair(numeric_limits<int>::max(), vertex2));
                }
            /* Relax edges outgoing from vertex1, whose other end-point hasn't been processed yet, and check if any of
             * them might belong to a shortest path from source to sink */
            for (const auto[vertex2, weight]: graph_1->at(vertex1)) {
                auto d_via_vertex1 = df_1 + weight;
                if (!in(vertex2, processed_1)) {
                    auto df_2 = d_1[vertex2];
                    if (d_via_vertex1 < df_2) {
                        d_1[vertex2] = d_via_vertex1;
                        pending_1.extract(make_pair(df_2, vertex2));
                        pending_1.emplace(make_pair(d_via_vertex1, vertex2));
                        pred_1[vertex2] = vertex1;
                    }
                }
                if (in(vertex2, processed_2)) {
                    auto length = d_via_vertex1 + d_2[vertex2];
                    if (length < shortest_so_far) {
                        shortest_so_far = length;
                        // Edge (best_vertex1, best_vertex2) is along the shortest path found so far from source to sink.
                        best_vertex1 = vertex1;
                        best_vertex2 = vertex2;
                        counter_at_shortest_path = step_counter;
                    }
                }
            }
        }
        // Check termination condition
        if (best_vertex2 != None && !pending_1.empty() && !pending_2.empty()) {
            const auto[l1, dummy1] = *pending_1.cbegin();
            const auto[l2, dummy2] = *pending_2.cbegin();
            if (l1 + l2 >= shortest_so_far)
                shortest_path_found = true;
        }
        /* Swap the information related to the two directions, as steps of the Dijkstra algorithms will alternate
         * between them. */
        if (!shortest_path_found) {
            swap(graph_1, graph_2);
            pending_1.swap(pending_2);
            d_1.swap(d_2);
            processed_1.swap(processed_2);
            pred_1.swap(pred_2);
        }
    }

    /* If you got here, a shortest path from source to sink was found, and it goes through edges best_vertex1 and
     * best_vertex2 (but not necessarily in that order). */
    assert(best_vertex1 != None);
    assert(best_vertex2 != None);

    /* Stitch together the shortest path as:
     * source -> ... -> vertex1 -> best_vertex2 -> ... -> sink''' */
    vector<int> sub_path_1, sub_path_2;
    if (step_counter % 2 == counter_at_shortest_path % 2) {
        sub_path_1 = backtrack_path(best_vertex1, pred_1); // This is reversed.
        sub_path_2 = backtrack_path(best_vertex2, pred_2);
    } else {
        sub_path_1 = backtrack_path(best_vertex2, pred_1); // This is reversed.
        sub_path_2 = backtrack_path(best_vertex1, pred_2);
    }
    std::reverse(sub_path_1.begin(), sub_path_1.end());
    vector<int> shortest_path = sub_path_1;
    shortest_path.insert(shortest_path.end(), sub_path_2.begin(), sub_path_2.end());
    if (step_counter % 2 == 1)
        std::reverse(shortest_path.begin(), shortest_path.end());

    // Compute the length of the shortest path, as the sum of the lengths of the two sub-paths.
    int distance = d_1[best_vertex2] + d_2[best_vertex2];
    assert(distance == shortest_so_far);

    const auto res = make_pair(distance, shortest_path);
    return res;
}





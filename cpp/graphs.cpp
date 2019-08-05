#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <stdexcept>
#include <limits>
#include <random>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#define CATCH_CONFIG_MAIN

#include "../include/catch.hpp"

using std::string;
using std::vector;
using std::ifstream;
using std::istream;
using std::pair;
using std::tuple;
using std::make_tuple;
using std::make_pair;
using std::cout;
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
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_weight_t, int>> BoostGraph;
typedef boost::graph_traits<BoostGraph>::vertex_descriptor Vertex;
typedef boost::property_map<BoostGraph, boost::vertex_index_t>::type IndexMap;

BoostGraph make_boost_graph(Graph &graph) {
    typedef pair<int, int> Edge;
    vector<Edge> edges;
    int n_edges = 0;
    unordered_set<int> vertices;
    for (const auto &[v1, adj_list]: graph)
        for (const auto &[v2, w]: adj_list) {
            ++n_edges;
            edges.emplace_back(make_pair(v1, v2));
            vertices.insert(v1);
            vertices.insert(v2);
        }
    vector<int> weights(n_edges, 1);
    const auto n_vertices{vertices.size()};
    BoostGraph boost_graph(edges.begin(), edges.end(), weights.begin(), n_vertices);
    return boost_graph;
}


struct TestCase {
    TestCase(Graph graph, vector<pair<int, int>> test_cases) : graph(std::move(graph)),
                                                               queries(std::move(test_cases)) {}

    Graph graph;
    vector<pair<int, int>> queries;
};


void append_adj(Graph &graph, int v1, int v2, int w) {
    if (graph.find(v1) == graph.end())
        graph[v1] = {{v2, w}};
    else
        graph[v1].emplace_back(make_pair(v2, w));
}


TestCase fetch_graph_test_case(const string &file_name) {
    ifstream input_file(file_name);
    if (!input_file.is_open()) {
        cout << "File not found: " << file_name << endl;
        throw std::invalid_argument("File not found: " + file_name);
    }

    int n_vertices, n_edges;
    input_file >> n_vertices >> n_edges;
    Graph graph;
    for (int i = 0; i < n_edges; ++i) {
        int v1, v2, w;
        input_file >> v1 >> v2 >> w;
        append_adj(graph, v1, v2, w);
    }
    int n_test_cases;
    input_file >> n_test_cases;
    vector<pair<int, int>> test_cases;
    for (int i = 0; i < n_test_cases; ++i) {
        int source, sink;
        input_file >> source;
        input_file >> sink;
        test_cases.emplace_back(make_pair(source, sink));
    }
    auto res = TestCase(graph, test_cases);

    return res;
}

vector<tuple<int, int, int>> fetch_twitter_test_case(const string & file_name) {
    ifstream input_file(file_name);
    if (!input_file.is_open()) {
        cout << "File not found: " << file_name << endl;
        throw std::invalid_argument("File not found: " + file_name);
    }

    vector<tuple<int, int, int>> res;
    while (!input_file.eof()) {
        int v1, v2, d;
        input_file >> v1 >> v2 >> d;
        res.emplace_back(make_tuple(v1, v2, d));
    }
    return res;
}


Graph fetch_social_media_combined(const string &file_name, const bool make_bidirectional = false) {
    ifstream input_file(file_name);
    if (!input_file.is_open()) {
        cout << "File not found: " << file_name << endl;
        throw std::invalid_argument("File not found: " + file_name);
    }

    Graph graph;
    while (!input_file.eof()) {
        int v1, v2;
        input_file >> v1 >> v2;
        append_adj(graph, v1, v2, 1);
        if (make_bidirectional)
            append_adj(graph, v2, v1, 1);
    }

    return graph;
}


template<typename Container, typename Value>
bool in(Value value, const Container &container) {
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


pair<int, vector<int>> bidir_dijkstra(const Graph &graph, int source, int sink) {
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
    int best_vertex1 = None;
    int best_vertex2 = None;
    int step_count;
    for (step_count = 0; !shortest_path_found; ++step_count) {
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
        if (in(vertex1, *graph_1)) {
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
        /* Trade the information related to the two directions, as steps of the Dijkstra algorithms will alternate
         * between them */
        swap(graph_1, graph_2);
        pending_1.swap(pending_2);
        d_1.swap(d_2);
        processed_1.swap(processed_2);
        pred_1.swap(pred_2);
    }

    /* If you got here, a shortest path from source to sink was found, and it goes through edge (best_vertex1,
     * best_vertex2). */
    assert(best_vertex1 != None);
    assert(best_vertex2 != None);

    if (step_count % 2 == 1) {
        swap(graph_1, graph_2);
        pending_1.swap(pending_2);
        d_1.swap(d_2);
        processed_1.swap(processed_2);
        pred_1.swap(pred_2);
    } else
        swap(best_vertex1, best_vertex2);

    /* Stitch together the shortest path as:
     * source -> ... -> vertex1 -> best_vertex2 -> ... -> sink''' */
    auto sub_path_1 = backtrack_path(best_vertex1, pred_1); // This is reversed.
    auto sub_path_2 = backtrack_path(best_vertex2, pred_2);
    std::reverse(sub_path_1.begin(), sub_path_1.end());
    vector<int> shortest_path = sub_path_1;
    shortest_path.insert(shortest_path.end(), sub_path_2.begin(), sub_path_2.end());

    // Compute the length of the shortest path, as the sum of the lengths of the two sub-paths.
    int distance = (step_count % 2 == 1) ? d_1[best_vertex2] + d_2[best_vertex2] : d_1[best_vertex1] +
                                                                                   d_2[best_vertex1];
    assert(distance == shortest_so_far);

    const auto res = make_pair(distance, shortest_path);
    return res;
}


TEST_CASE("in") {
    unordered_set<int> s = {1, 3, 5, 7, 9};
    REQUIRE(in(1, s));
    REQUIRE(in(3, s));
    REQUIRE(in(5, s));
    REQUIRE(in(7, s));
    REQUIRE(in(9, s));
    REQUIRE(!in(0, s));
    REQUIRE(!in(4, s));
}


TEST_CASE("bidirectional_dijkstra") {

    SECTION("small input") {
        auto test_case = fetch_graph_test_case("../../test/test03.txt");
        for (auto query: test_case.queries) {
            auto[distance, path] = bidir_dijkstra(test_case.graph, query.first, query.second);
            REQUIRE(distance == 5);
            REQUIRE(path == vector<int>({0, 1, 4, 5}));
        }

        test_case = fetch_graph_test_case("../../test/test04.txt");
        for (auto query: test_case.queries) {
            auto[distance, path] = bidir_dijkstra(test_case.graph, query.first, query.second);
            REQUIRE(distance == 5);
            REQUIRE(path == vector<int>({0, 1, 3, 4}));
        }

        test_case = fetch_graph_test_case("../../test/test05.txt");
        for (auto query: test_case.queries) {
            auto[distance, path] = bidir_dijkstra(test_case.graph, query.first, query.second);
            REQUIRE(distance == 1);
            REQUIRE(path == vector<int>({0, 1}));
        }

        test_case = fetch_graph_test_case("../../test/test06.txt");
        for (auto query: test_case.queries) {
            auto[distance, path] = bidir_dijkstra(test_case.graph, query.first, query.second);
            REQUIRE(distance == 3);
            REQUIRE(path == vector<int>({0, 1, 2}));
        }

        test_case = fetch_graph_test_case("../../test/test01.txt");
        vector<int> expected = {0, 0, 1, -1};
        for (int i = 0; i < test_case.queries.size(); ++i) {
            auto query = test_case.queries[i];
            auto[distance, path] = bidir_dijkstra(test_case.graph, query.first, query.second);
            REQUIRE(distance == expected[i]);
        }

        test_case = fetch_graph_test_case("../../test/test02.txt");
        for (auto query: test_case.queries) {
            auto[distance, path] = bidir_dijkstra(test_case.graph, query.first, query.second);
            REQUIRE(distance == 3);
            REQUIRE(path == vector<int>({1, 2, 3}));
        }
    }

    SECTION("facebook") {

        string file_name = "../../test/facebook_combined.txt";
        auto graph = fetch_social_media_combined(file_name, true);
        std::random_device rd;
        std::mt19937 generator(rd());
        generator.seed(42);
        std::uniform_int_distribution<> distribution(0, 4031);
        auto boost_graph = make_boost_graph(graph);
        std::vector<int> d(num_vertices(boost_graph));

        int count = 0;
        for (int i = 0; i < 1; ++i) {
            auto source = distribution(generator);
            auto sink = distribution(generator);
            auto boost_source = vertex(source, boost_graph);
            auto boost_sink = vertex(sink, boost_graph);
            dijkstra_shortest_paths(boost_graph, boost_source, boost::distance_map(&d[0]));
            auto boost_distance = d[boost_sink];
            if (boost_distance == numeric_limits<int>::max()) {
                boost_distance = -1;
                ++count;
            }
            auto[distance, path] = bidir_dijkstra(graph, source, sink);
            REQUIRE(distance == boost_distance);
            cout << distance << " ";
        }
        cout << "\ncount is " << count << endl;
    }

    SECTION("twitter") {

        string file_name = "../../test/twitter_combined.txt";
        auto graph = fetch_social_media_combined(file_name);
        auto boost_graph = make_boost_graph(graph);

        std::random_device rd;
        std::mt19937 generator(rd());
        generator.seed(42);
        auto n_vertices = num_vertices(boost_graph);
        std::uniform_int_distribution<> distribution2(0, n_vertices - 1);
        std::vector<int> d2(n_vertices);
        std::vector<int> vertices(n_vertices);
        IndexMap index = get(boost::vertex_index, boost_graph);
        int pos = 0;
        for (auto vi = boost::vertices(boost_graph).first; vi != boost::vertices(boost_graph).second; ++vi, ++pos)
            vertices[pos] = index(*vi);
        assert(pos + 1 == n_vertices);

        int count = 0;
        for (int i = 0; i < 1; ++i) {
            auto source = vertices[distribution2(generator)];
            auto sink = vertices[distribution2(generator)];
            auto boost_source = vertex(source, boost_graph);
            auto boost_sink = vertex(sink, boost_graph);
            dijkstra_shortest_paths(boost_graph, boost_source, boost::distance_map(&d2[0]));
            auto boost_distance = d2[boost_sink];
            if (boost_distance == numeric_limits<int>::max()) {
                boost_distance = -1;
                ++count;
            }
            auto[distance, path] = bidir_dijkstra(graph, source, sink);
            REQUIRE(distance == boost_distance);
            cout << distance << " ";
        }
        cout << "\ncount is " << count << endl;
    }

}

/* TODO
 * Add twitter test cases
 * add download of social media datasets
 * add in-line documentation
 * check the timing (for fun)
 * add portable paths, using Boost Filesystem
 * */

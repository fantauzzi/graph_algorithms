#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <limits>
#include <random>
// c++ version 7.4.0 still has "filesystem" under "experimental", even when compiling for C++17
#include <experimental/filesystem>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <curl/curl.h>
#include <cstdio>
#include <ctime>

#define CATCH_CONFIG_MAIN

#include "../include/catch.hpp"

#include "graphs.h"

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
using std::experimental::filesystem::path;
using std::experimental::filesystem::exists;
using std::experimental::filesystem::remove;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_weight_t, int>> BoostGraph;

struct TestCase {
    TestCase(Graph graph, vector<pair<int, int>> test_cases) : graph(std::move(graph)),
                                                               queries(std::move(test_cases)) {}

    Graph graph;
    vector<pair<int, int>> queries;
};


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


TestCase fetch_graph_test_case(const string &file_name) {
    ifstream input_file(path(file_name).native());
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


vector<tuple<int, int, int>> fetch_twitter_test_case(const string &file_name) {
    ifstream input_file(path(file_name).native());
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
    ifstream input_file(path(file_name).native());
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


// Useful to download a file via curl.
static size_t write_data(void *ptr, size_t size, size_t nmemb, void *stream) {
    size_t written = fwrite(ptr, size, nmemb, (FILE *) stream);
    return written;
}


bool fetch_as_needed(const string &file_name, const string &url) {
    auto native_name = path(file_name).native();
    auto gz_native_name = path(file_name + ".gz").native();
    if (!exists(native_name)) {
        if (!exists(gz_native_name)) {
            cout << "Donwloading " << gz_native_name << endl;
            curl_global_init(CURL_GLOBAL_ALL);
            auto curl_handle = curl_easy_init();
            curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str());
            // curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 1L);
            curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);

            auto pagefile = std::fopen(gz_native_name.c_str(), "wb");
            if (pagefile) {
                curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, pagefile);
                curl_easy_perform(curl_handle);
                fclose(pagefile);
            }
            curl_easy_cleanup(curl_handle);
            curl_global_cleanup();
        }
        if (exists(gz_native_name)) {
            cout << "Unzipping " << gz_native_name << " into " << native_name << endl;
            std::ofstream output_file;
            output_file.open(native_name);
            ifstream file(gz_native_name, std::ios_base::in | std::ios_base::binary);
            boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
            in.push(boost::iostreams::gzip_decompressor());
            in.push(file);
            boost::iostreams::copy(in, output_file);
            output_file.close();
            if (exists(gz_native_name))
                remove(gz_native_name);
        }
    }
    return exists(native_name);
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
    /* Downloads datasets from "SNAP Datasets: Stanford Large Network Dataset Collection
     * Jure Leskovec and Andrej Krevl
     * http://snap.stanford.edu/data
     * June, 2014
     */

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
        bool ok = fetch_as_needed(file_name, "https://snap.stanford.edu/data/facebook_combined.txt.gz");
        REQUIRE(ok);
        auto graph = fetch_social_media_combined(file_name, true);
        std::random_device rd;
        std::mt19937 generator(rd());
        generator.seed(42);
        std::uniform_int_distribution<> distribution(0, 4031);
        auto boost_graph = make_boost_graph(graph);
        std::vector<int> d(num_vertices(boost_graph));

        for (int i = 0; i < 100; ++i) {
            auto source = distribution(generator);
            auto sink = distribution(generator);
            auto boost_source = vertex(source, boost_graph);
            auto boost_sink = vertex(sink, boost_graph);
            dijkstra_shortest_paths(boost_graph, boost_source, boost::distance_map(&d[0]));
            auto boost_distance = d[boost_sink];
            if (boost_distance == numeric_limits<int>::max()) {
                boost_distance = -1;
            }
            auto[distance, path] = bidir_dijkstra(graph, source, sink);
            REQUIRE(distance == (path.empty() ? -1 : path.size() - 1));
            REQUIRE(distance == boost_distance);
        }
    }

    SECTION("twitter") {
        string file_name = "../../test/twitter_combined.txt";
        bool ok = fetch_as_needed(file_name, "https://snap.stanford.edu/data/twitter_combined.txt.gz");
        REQUIRE(ok);
        auto graph = fetch_social_media_combined(file_name);
        auto test_cases = fetch_twitter_test_case("../../test/twitter_tcs.txt");

        for (const auto[source, sink, expected]: test_cases) {
            auto[distance, path] = bidir_dijkstra(graph, source, sink);
            REQUIRE(distance == (path.empty() ? -1 : path.size() - 1));
            REQUIRE(distance == expected);
        }
    }
}

/* TODO
 * check the timing (for fun)
  * */

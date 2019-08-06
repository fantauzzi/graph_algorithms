#pragma once

#include <vector>
#include <unordered_map>
#include <utility>

typedef std::vector<std::pair<int, int>> AdjList;
typedef std::unordered_map<int, AdjList> Graph;

void append_adj(Graph &graph, int v1, int v2, int w);

template<typename Container, typename Value>
bool in(Value value, const Container &container);

std::pair<int, std::vector<int>> bidir_dijkstra(const Graph &graph, int source, int sink);

#include <Graph.hpp>
#include "dijkstra.hpp"
#include <list>

template <typename Graph,typename Vector>
Dijkstra<Graph,Vector>::Dijkstra() {
    _M = nullptr;
    _R = nullptr;
    _U = nullptr;
    _graph = nullptr;
    _shortestPath = nullptr;
}

template <typename Graph, typename Vector>
Dijkstra<Graph,Vector>::Dijkstra(Graph graph, int startIdx) {
    _graph = graph;
    _M.push_back(startIdx);
    _R = _graph.getNeighbour();

}
#ifndef DIJKSTRA_HPP
#define DIJKSTRA_HPP

#include <list>

template <typename Graph, typename Vector>
class Dijkstra{
private:
    std::list<int> _M;
    std::list<int> _R;
    std::list<int> _U;
    Graph _graph;
    Vector _shortestPath;

public:
    //default constructor
    Dijkstra();

    //constructor
    Dijkstra(Graph graph, int startIdx);




};



#endif
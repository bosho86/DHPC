//
// Created by Claudio Cannizzaro on 2019-10-22.
//

#ifndef DHPC_GRAPH_HPP
#define DHPC_GRAPH_HPP

#include <vector>

class Graph {
private:
    std::vector<int> adjMat;
    int numVertices;

public:
    Graph(int numVertex) {
        adjMat = std::vector<int>(numVertex * numVertex);
        numVertices = numVertex;
    }

    //default constructor, creates a predefined Graph, directed, only positive edge weights
    Graph() {
        adjMat = std::vector<int>(6 * 6);
        adjMat[0 * 6 + 1] = 2;
        adjMat[0 * 6 + 2] = 3;
        adjMat[1 * 6 + 2] = 2;
        adjMat[1 * 6 + 3] = 6;
        adjMat[2 * 6 + 0] = 3;
        adjMat[2 * 6 + 1] = 2;
        adjMat[2 * 6 + 4] = 1;
        adjMat[3 * 6 + 1] = 6;
        adjMat[3 * 6 + 4] = 3;
        adjMat[3 * 6 + 5] = 1;
        adjMat[4 * 6 + 2] = 1;
        adjMat[4 * 6 + 3] = 3;
        adjMat[4 * 6 + 5] = 1;
        adjMat[5 * 6 + 3] = 1;

        numVertices = 6;
    }

    std::vector<int> getAdjMat() {
        return adjMat;
    }

    int getNumVertex(){
        return numVertices;
    }
};

#endif //DHPC_GRAPH_HPP

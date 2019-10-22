//
// Created by Claudio Cannizzaro on 2019-10-22.
//

#include "Graph.hpp"
#include <iostream>

int main() {
    Graph a = Graph();

    int i = 0;
    int size = a.getNumVertex();
    for (auto x:a.getAdjMat()) {
        std::cout << x << " ";
        if ((i++ % size) == 5) {
            std::cout << std::endl;
        }
    }
}
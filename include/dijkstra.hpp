#ifndef DIJKSTRA_HPP
#define DIJKSTRA_HPP

#include <list>
#include <vector>
#include <set>
#include <algorithm>
#include <typeinfo>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <numeric>
#include <fstream>
#include <complex>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstring>
#include <new>
#include <cmath>
#include <functional>

#include <read_input_file.h>
#include <Graph.hpp>

template <typename T>
class Dijkstra{
private:
    Graph<T> graph;
    std::vector<T> AdjMat;  
    std::vector< std::vector< std::pair<T ,int> > > Vertex_cost_and_NN;
    std::vector<T> d;
    std::vector<int> pred;
    int Nvertices;
public:
    //default constructor
    Dijkstra();

    //constructor
    Dijkstra(Graph<T> &G){
        graph = G;        
        AdjMat = graph.getAdjMat();
        Nvertices = graph.getNumVertex();
        Vertex_cost_and_NN = graph.getVertex_cost_and_NN();
        d.resize(Nvertices);
        pred.resize(Nvertices);
        
        for(int u=0; u<Nvertices; u++){
            d[u] = infinity;
            pred[u] = -1;
        }

    }
    
    void Find_path_Dijkstra(const int &source, const int &target, std::vector<int> &path){
    
    
        d[source] = 0;
        pred[source] = -1;
        std::set< std::pair<T, int> > Q_queue;
        Q_queue.insert ( std::make_pair(d[source], source) );
        
        while(!Q_queue.empty()){
        
            typename std::set<std::pair<T, int>>::iterator iter = Q_queue.begin();
            int u = ( *iter ).second;
            Q_queue.erase( Q_queue.begin());
            
            
            int ns = 0;
            //for(int ic = 0; ic<Nvertices; ic++){
                
            //    if(AdjMat[u * Nvertices + ic] < infinity)
            //        ns++;
            
            //}
            //std::cout << ns << std::endl;
            ns = Vertex_cost_and_NN[u].size();
            for(int v=0; v< ns; v++ ){
                
                T weight = Vertex_cost_and_NN[u][v].first; //AdjMat[u * Nvertices + v];
                int ver = Vertex_cost_and_NN[u][v].second;
                
                if( d[u] + weight < d[ver] ){
                    
                    Q_queue.erase(std::make_pair(d[ver], ver));
                    
                    d[ver] = d[u] + weight;
                    pred[ver] = u;
                    Q_queue.insert(std::make_pair(d[ver], ver));
                    
                }
                
                
            }
        }
        
        for (int vertex = 0;  vertex <Nvertices; vertex++)
            if(pred[vertex]>-1)
                path.push_back(pred[vertex]);
        
        printSolution(d, source, pred);

    }
    
    // Function to print shortest 
    // path from source to j 
    // using parent array 
    void printPath(std::vector<int> &parent, int j) 
    { 

        // Base Case : If j is source 
        if (parent[j] == -1) 
            return; 
        
        int i = parent[j];
        printPath(parent, i); 

        printf("%d ", j); 
    } 

    // A utility function to print  
    // the constructed distance 
    // array 
    void printSolution(std::vector<T> & dist, int src, std::vector<int> & parent) 
    { 
        //int src = 1; 
        printf("Vertex\t Distance\tPath"); 
        for (int i = 0; i < Nvertices; i++) 
        { 
            if(i != src){
                printf("\n%d -> %d \t\t %g\t\t%d ", 
                              src, i, dist[i], src); 
                printPath(parent, i); 
            }
        } 
    } 



};



#endif
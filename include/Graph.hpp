//
// Created by Claudio Cannizzaro on 2019-10-22.
//

#ifndef DHPC_GRAPH_HPP
#define DHPC_GRAPH_HPP

#include <vector>
#include <valarray>
#define infinity 1E6

template <typename T>
class Graph {
private:
    std::vector<T> adjMat;
    int numVertices;
    std::vector< std::pair<T ,T> > vertex_X_Y;
    std::vector< std::vector< std::pair<T ,int> > > vertex_cost_and_NN;
public:
//    Graph(int numVertex) {
//        adjMat = std::vector<int>(numVertex * numVertex);
//        numVertices = numVertex;
//    }

    //default constructor, creates a predefined Graph, directed, only positive edge weights
    Graph() {
        adjMat = std::vector<T>(6 * 6);
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

    std::vector<T> getAdjMat() {
        return adjMat;
    }

    int getNumVertex(){
        return numVertices;
    }

    std::vector< std::vector< std::pair<T ,int> > > getVertex_cost_and_NN() {
        return vertex_cost_and_NN;
    }    
    
    std::vector< std::pair<T ,T> > getVertex_X_Y() {
        return vertex_X_Y;
    }
    
    std::pair<T ,T> getX_Y_at_Vertex(int vertex) {
        return vertex_X_Y[vertex];
    }    
    
    //get index of neighbours of vertex at idx
    std::vector<T> getNeighbour(int idx){
        std::vector<T> neighbours;
        for (int i = 0; i < numVertices; ++i) {
            if (adjMat[idx * numVertices + i] != 0){
                neighbours.push_back(i);
            }
        }
        return neighbours;
    }
    
    // added by Hamilton

    Graph(int numVertex, bool read_map, std::vector< std::vector<T> > MAP) {
        
        if (!read_map){
        
            map_generation(numVertex);
        
        }else{
            
            reading_map_from_file(MAP);
            
        }
       
        
    }
    
    void reading_map_from_file( std::vector< std::vector<T> > MAP){

        numVertices = MAP.size();
        adjMat.resize(numVertices * numVertices);

        // Assuming is a 2D map
        // randomly generate vertices position (x,y)


        std::vector<std::vector<int> > NN;
        NN.resize(numVertices);

        for (int i = 0; i < numVertices; i++) {

            int ns = MAP[i].size() - 3;
            NN[i].resize(ns);
            
            //std::cout<<  ns << " " << MAP[i].size() <<std::endl;

            T ix = MAP[i][1];
            T iy = MAP[i][2];
            vertex_X_Y.push_back( std::make_pair(ix, iy) );
            //std::cout << ix << " " << iy << " " << i << std::endl;

            for (int in = 0; in < ns; in++){
                NN[i][in] = (int) (MAP[i][3 + in]);
                //std::cout<<  i << " " <<  MAP[i][3 + in] << std::endl;
            }
        }

        // computing the cost or weight. They are defined as the 
        // Euclidean distance between nodes. Then, store in adjMat

        for (int i = 0; i < numVertices; i++) {
            for (int j = 0; j < numVertices; j++) {
                adjMat[i * numVertices + j] = T(infinity);
            }
        }
        
        for (int i = 0; i < numVertices; i++) {

            
            int ns = NN[i].size();
            std::vector< std::pair<T ,int> > temp;

            for (int j = 0; j < ns; j++) {

                int nn = NN[i][j];
                
                if(nn > -1 ){
                    
                    T dx = vertex_X_Y[i].first - vertex_X_Y[nn].first;
                    T dy = vertex_X_Y[i].second - vertex_X_Y[nn].second;

                    T dist = T(std::sqrt(dx * dx + dy * dy));
                    adjMat[i * numVertices + nn ] = dist;
                    temp.push_back(std::make_pair(dist,  nn));                
                
                }else{
                
                    T dist = T(1e6);
                    adjMat[i * numVertices + nn ] = dist;
                    temp.push_back(std::make_pair(dist,  nn));                            
                
                }
                
            }
            
            vertex_cost_and_NN.push_back(temp);

        }
        
        
    }

    // not tested but could be useful.
    // Random generation of a map    
    void map_generation(int numVertex){
        
        numVertices = numVertex;
        adjMat.resize(numVertex*numVertex);       
        
        // Nx and Ny are the boundary of the Map
        // Arbitrarily assume max distance between nodes to be some
        // value. In this case is 30.
        int Nx =100; int Ny = 100;
        int max_dist_nodes = 30;
               

        
        std::srand((unsigned)time(NULL)); 
        
        // Assuming is a 2D map
        // randomly generate vertices position (x,y)
        
        for(int i = 0; i < numVertex; i++){
            
            int ix = std::rand()%(Nx);
            int iy = std::rand()%(Ny);
            vertex_X_Y.push_back(std::make_pair(ix, iy));
            
        }
        
        // computing the cost or weight. They are defined as the 
        // Euclidean distance between nodes. Then, store in adjMat
        
        for( int i = 0; i<numVertex; i++){
            
            for (int j= 0; j<numVertex; j++){
                
                adjMat[i*numVertex + j] = 0;
                
                if( j != i ){
                    
                    int dx = vertex_X_Y[i].first -  vertex_X_Y[j].first;
                    int dy = vertex_X_Y[i].second -  vertex_X_Y[j].second;
                    int dist = std::floor(std::sqrt( dx*dx + dy*dy  ));
                    
                    if (dist <= max_dist_nodes){
                        
                        adjMat[i*numVertex + j] = dist;
                    
                    }
                    
                    
                    
                }
                
            }
                    
        }        
    
    }
    
};

#endif //DHPC_GRAPH_HPP

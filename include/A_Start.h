#include <list>
#include <algorithm>
#include <typeinfo>
#include <iostream>
#include <stdio.h>
#include <vector>
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
#include "Graph.hpp"

// initializing class point
template <typename T>
class point {
public:
    point(){x = T(0); y = T(0);};
    point( T a, T b ) { x = a; y = b; }
    bool operator ==( const point<T>& o ) { return o.x == x && o.y == y; }
    point operator +( const point<T>& o ) { return point<T>( o.x + x, o.y + y ); }
    T x, y;
};

template <typename T>
class map {
private:     
    Graph<T> graph;
    std::vector<T> AdjMat;
public:
    map(Graph<T> &G){
        graph = G;        
        AdjMat = graph.getAdjMat();
        w = graph.getNumVertex();
        h = graph.getNumVertex();
        
        for( int r = 0; r < h; r++ ){
            for( int s = 0; s < w; s++ ){
                m[s][r] = 0;        
            }
        }
        
    }
    map() {
        char t[8][8] = {
            {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 1, 1, 1, 0}, {0, 0, 1, 0, 0, 0, 1, 0},
            {0, 0, 1, 0, 0, 0, 1, 0}, {0, 0, 1, 1, 1, 1, 1, 0},
            {0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}
        };

        w = h = 8;
        for( int r = 0; r < h; r++ )
            for( int s = 0; s < w; s++ )
                m[s][r] = t[r][s];
         
    }
    int operator() ( int ix, int iy ) { return m[ix][iy]; }
    char m[8][8];
    int w, h;
    
    std::vector<T> getMap_Graph(){ 
        
        w = graph.getNumVertex();
        h = graph.getNumVertex();
        return graph.getAdjMat();        
                
    };
};

template <typename T>
class node {
public:
    bool operator == (const node<T>& o ) { return pos == o.pos; }
    bool operator == (const point<T>& o ) { return pos == o; }
    bool operator < (const node<T>& o ) { return dist + cost < o.dist + o.cost; }
    point<T> pos, parent;
    int dist, cost;
    int isnode;
};

template <typename T>
class aStar {
private:     
    Graph<T> graph;
    std::vector<T> vertex_x_y;
    std::vector< std::vector< std::pair<T ,int> > > Vertex_cost_and_NN;
public:
    aStar(Graph<T> &G){
    
        graph = G;
        Vertex_cost_and_NN = graph.getVertex_cost_and_NN();
        vertex_x_y  =  graph.getVertex_X_Y();
    
    }
    aStar() {
        
        neighbours.resize(8);
        neighbours[0] = point<T>( -1, -1 ); neighbours[1] = point<T>(  1, -1 );
        neighbours[2] = point<T>( -1,  1 ); neighbours[3] = point<T>(  1,  1 );
        neighbours[4] = point<T>(  0, -1 ); neighbours[5] = point<T>( -1,  0 );
        neighbours[6] = point<T>(  0,  1 ); neighbours[7] = point<T>(  1,  0 );
        
    }
 
    
    
    int calcDist( point<T>& p ){
        // need a better heuristic
        int x = end.x - p.x, y = end.y - p.y;
        return( x * x + y * y );
    }
 
    bool isValid( point<T>& p ) {
        return ( p.x >-1 && p.y > -1 && p.x < m.w && p.y < m.h );
    }
 
    bool existPoint( point<T>& p, T cost ) {
        typename std::list<node<T>>::iterator i;
        i = std::find( closed.begin(), closed.end(), p );
        
        if( i != closed.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost ) return true;
            else { closed.erase( i ); return false; }
        }
        
        i = std::find( open.begin(), open.end(), p );
        if( i != open.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost ) return true;
            else { open.erase( i ); return false; }
        }
        return false;
    }
 
    bool fillOpen( node<T>& n ) {
        
        int stepCost, nc, dist;
        point<T> neighbour;
 
        for( int x = 0; x < 8; x++ ) {
            // one can make diagonals have different cost
            stepCost = 1; // x < 4 ? 1 : 1; // here cost
            //inn = Vertex_cost_and_NN[x]
            //T p_x = vertex_x_y[x].first;
            //point<T> NN();
            neighbour = n.pos + neighbours[x];
            if( neighbour == end ) return true;
 
            if( isValid( neighbour ) && m( neighbour.x, neighbour.y ) != 1 ) {
                nc = stepCost + n.cost;
                dist = calcDist( neighbour );
                if( !existPoint( neighbour, nc + dist ) ) {
                    node<T> m;
                    m.cost = nc; m.dist = dist;
                    m.pos = neighbour; 
                    m.parent = n.pos;
                    open.push_back( m );
                }
            }
        }
        return false;
    }
 
    bool search( point<T>& s, point<T>& e, map<T>& mp ) {
        point<T> Zero;
        node<T> n; end = e; start = s; m = mp;
        
        n.cost = 0; n.pos = s; n.parent = Zero; n.dist = calcDist( s ); 
        
        open.push_back( n );
        while( !open.empty() ) {
            //open.sort();
            node<T> n = open.front();
            open.pop_front();
            closed.push_back( n );
            if( fillOpen( n ) ) return true;
        }
        return false;
    }
 
    int path( std::list<point<T>>& path ) {
        
        path.push_front( end );
        int cost = 1 + closed.back().cost;
        
        path.push_front( closed.back().pos );
        point<T> parent = closed.back().parent;
 
        for(typename std::list<node<T>>::reverse_iterator i = closed.rbegin(); i != closed.rend(); i++ ) {
            if( ( *i ).pos == parent && !( ( *i ).pos == start ) ) {
                path.push_front( ( *i ).pos );
                parent = ( *i ).parent;
            }
        }
        path.push_front( start );
        return cost;
    }
 
    map<T> m; point<T> end, start;
    std::vector< point<T> > neighbours;
    std::list< node<T> > open;
    std::list< node<T> > closed;
};
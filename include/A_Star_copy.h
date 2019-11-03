/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   A_Star_copy.h
 * Author: hamilton
 *
 * Created on November 1, 2019, 12:14 PM
 */

#ifndef A_STAR_COPY_H
#define A_STAR_COPY_H

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   A_Start.cpp
 * Author: hamilton
 * 
 * Created on October 24, 2019, 10:17 AM
 */

#include <list>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>

#include "Graph.hpp"

// initializing class point
template <typename T>
class point {
public:
    /* initializing x and y  */
    //point();
    point(T a=0, T b=0);
    /*Overloading operators for point class objects*/
    bool operator ==( const point<T>& o ) { return o.x == x && o.y == y; }
    point operator +( const point<T>& o ) { return point<T>( o.x + x, o.y + y ); }
    T x, y;
    virtual ~point();
};

template <typename T>
point<T>::point( T a, T b ) 
{ x = a; y = b; }

template <class Graph, typename T>
class map {
public:
    /*Initializing map. This map can be changed.*/
    map();
    map(Graph matrix);
    int operator() ( int x, int y ) { return m[x][y]; }
    std::vector<T> map_graph;
    char m[8][8];
    int w, h;
    
    virtual ~map();
};

template <class Graph, typename T>
map<Graph, T>::map() {
    char t[8][8] = {
        {0, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 1, 1, 0},
        {0, 0, 1, 0, 0, 0, 1, 0},
        {0, 0, 1, 0, 0, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0}
    };
    w = h = 8;
    for( int r = 0; r < h; r++ )
        for( int s = 0; s < w; s++ )
            m[s][r] = t[r][s];
}
 

template <class Graph, typename T>
map<Graph, T>::map(Graph matrix) {
    
    w = matrix.getnumVertices();
    h = matrix.getnumVertices();
    map_graph.resize(w*h);
    map_graph = matrix.getAdjMat();
    
}

template <typename Tp, typename Tn>
class node {
public:
    /*Overloading operators for node class objects*/
    node();
    bool operator == (const node<Tp, Tn>& o ) { return pos == o.pos; }
    bool operator == (const point<Tp>& o ) { return pos == o; }
    bool operator < (const node<Tp, Tn>& o ) { return dist + cost < o.dist + o.cost; }
    point<Tp> pos, parent;
    Tn dist, cost;
    virtual ~node();
};
 
/*aStar class contains all function for the search algorithm*/
template <typename Tp, typename Tn>
class aStar {
public:

    aStar();
    aStar(point<Tp>& p);
    Tn calcDist( point<Tp>& p );
    bool isValid( point<Tp>& p );
    bool existPoint( point<Tp>& p, Tn cost );
    bool fillOpen(node<Tp, Tn>& n);
    bool search(point<Tp>& s, point<Tp>& e, map<Graph, Tp>& mp);
    Tn path(std::list< point<Tp> > &path);
    
    map<Graph, Tp> m; point<Tp> end, start;
    std::vector< point<Tp> > neighbours(8);
    std::list<node<Tp, Tn>> open;
    std::list<node<Tp, Tn>> closed;
    
    virtual ~aStar();    
    
}; 


// initializing aStar class
template <typename Tp, typename Tn>
aStar<Tp, Tn>::aStar( point<Tp>& p ) {
    /*neighbours with respect to cell (0,0)*/
    
    /*  0   5   2
        4   C   6 
        1   7   3 */
    
    neighbours[0] = p( -1, -1 ); neighbours[1] = p(  1, -1 );
    neighbours[2] = p( -1,  1 ); neighbours[3] = p(  1,  1 );
    neighbours[4] = p(  0, -1 ); neighbours[5] = p( -1,  0 );
    neighbours[6] = p(  0,  1 ); neighbours[7] = p(  1,  0 );
}

// Ecluidian heuristic distant between points or vertex 
template <typename Tp, typename Tn>
Tn aStar<Tp, Tn>::calcDist( point<Tp>& p ){
    // need a better heuristic
    int x = end.x - p.x, y = end.y - p.y;
    return ( x * x + y * y );
}

// Checking if is a valid point within the boundaries
template <typename Tp, typename Tn>
bool aStar<Tp, Tn>::isValid( point<Tp>& p ) {
    return ( p.x >-1 && p.y > -1 && p.x < m.w && p.y < m.h );
}

template <typename Tp, typename Tn>
bool aStar<Tp, Tn>::existPoint( point<Tp>& p, Tn cost ) {
    typename std::list< node<Tp, Tn> >::iterator i;
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

template <typename Tp, typename Tn>
bool aStar<Tp, Tn>::fillOpen(node<Tp,Tn>& n) {
    int stepCost, nc, dist;
    point<Tp> neighbour;

    for (int x = 0; x < 8; x++) {
        // one can make diagonals have different cost
        stepCost = x < 4 ? 1 : 1;
        neighbour = n.pos + neighbours[x];
        if (neighbour == end) return true;

        if (isValid(neighbour) && m(neighbour.x, neighbour.y) != 1) {
            nc = stepCost + n.cost;
            dist = calcDist(neighbour);
            if (!existPoint(neighbour, nc + dist)) {
                node<Tp, Tn> m;
                m.cost = nc;
                m.dist = dist;
                m.pos = neighbour;
                m.parent = n.pos;
                open.push_back(m);
            }
        }
    }
    return false;
}

template <typename Tp, typename Tn>
bool aStar<Tp, Tn>::search( point<Tp>& s, point<Tp>& e, map<Graph, Tp>& mp) {
    node<Tp, Tn> n;
    end = e;
    start = s;
    m = mp;
    n.cost = 0;
    n.pos = s;
    n.parent = 0;
    n.dist = calcDist(s);
    open.push_back(n);
    while (!open.empty()) {
        //open.sort();
        node<Tp, Tn> n = open.front();
        open.pop_front();
        closed.push_back(n);
        if (fillOpen(n)) return true;
    }
    return false;
}

template <typename Tp, typename Tn>
Tn aStar<Tp, Tn>::path(std::list<point<Tp>>& path) {
    path.push_front(end);
    int cost = 1 + closed.back().cost;
    path.push_front(closed.back().pos);
    point<Tp> parent = closed.back().parent;
    
    typename std::list< node<Tp, Tn> >::reverse_iterator i;
    for (i = closed.rbegin(); i != closed.rend(); i++) {
        if ((*i).pos == parent && !((*i).pos == start)) {
            path.push_front((*i).pos);
            parent = (*i).parent;
        }
    }
    path.push_front(start);
    return cost;
}



#endif /* A_STAR_COPY_H */


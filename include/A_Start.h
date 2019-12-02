#include <list>
#include <vector>
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

void mpi_loop_splitter(int *Nloop, int *loopmin, int *loopmax){

    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    div_t divider;
    int quot,rem;
    
    divider = div(*Nloop,size);
    quot = divider.quot;
    rem = divider.rem;
    
    *loopmin = rank*quot;
    *loopmax = *loopmin + quot - 1;
    
    if(rank == size-1)
        *loopmax+=rem;
    
}

// initializing class point
template <typename T>
class point {
public:

    
    point(){
    };
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
        isBarrier.resize(w*h);
        
        for( int r = 0; r < h; r++ ){
            for( int s = 0; s < w; s++ ){
                isBarrier[r*w+s] = 0;        
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
    int operator() ( int inode ) { return isBarrier[inode]; }
    char m[8][8];
    int w, h;
    std::vector<int> isBarrier;
    
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
    T dist, cost;
    int isnode;
};

template <typename T>
class aStar {
private:     
    Graph<T> graph;
    std::vector< std::pair<T ,T> > vertex_x_y;
    std::vector< std::vector< std::pair<T ,int> > > Vertex_cost_and_NN;
    int numVertex;
public:
    aStar(Graph<T> &G){
    
        graph = G;
        Vertex_cost_and_NN = graph.getVertex_cost_and_NN();
        numVertex = graph.getNumVertex();
        //vertex_x_y  =  graph.getVertex_X_Y();
        
        //std::cout << " x = " << vertex_x_y[0].first << " y = " << vertex_x_y[0].second << std::endl;
        //exit(0);
        
    
    }
    aStar() {
        
        neighbours.resize(8);
        neighbours[0] = point<T>( -1, -1 ); neighbours[1] = point<T>(  1, -1 );
        neighbours[2] = point<T>( -1,  1 ); neighbours[3] = point<T>(  1,  1 );
        neighbours[4] = point<T>(  0, -1 ); neighbours[5] = point<T>( -1,  0 );
        neighbours[6] = point<T>(  0,  1 ); neighbours[7] = point<T>(  1,  0 );
        
    }
    
    struct Node_Pos{
    
        double x;
        double y;
        
    };
    
    struct Node{

        int isnode;        
        double cost; 
        double dist;
        double pos_x;
        double pos_y;
        double parent_x;
        double parent_y;
 
    };                
    
    std::pair<T ,T> getX_Y_at_point( int inode ) { 
        
        std::pair<T ,T> pair_x_y = graph.getX_Y_at_Vertex(inode); 
        return pair_x_y;
        
    }    
    
    T calcDist( point<T>& p ){
        // need a better heuristic
        T x = end.x - p.x, y = end.y - p.y;
        return( x * x + y * y );
    }
 
    bool isValid( point<T>& p ) {
        
        //if( p.x >-1 && p.y > -1 && p.x < m.w && p.y < m.h )
        //    std::cout<< " isValid False " << std::endl; 
        
        return ( p.x >-1 && p.y > -1 && p.x < m.w && p.y < m.h );
    }
 
    bool existPoint( point<T>& p, T cost ) {
        typename std::list<node<T>>::iterator i;
//        typename std::vector<node<T>>::iterator j;
        i = std::find( closed.begin(), closed.end(), p );
//        j = std::find( Closed_vect_local.begin(), Closed_vect_local.end(), p );
        
        //std::cout<< "  cost " << ( *i ).cost + ( *i ).dist << "  " << cost <<std::endl;
        
        if( i != closed.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost && ( *i ).cost + ( *i ).dist > 1e-6 ){ 
                //std::cout<< "  closed " << ( *i ).isnode << "  " << ( *i ).cost + ( *i ).dist << "  " << cost <<std::endl;
                return true;}
            else { 
                closed.erase( i ); 
//                Closed_vect_local.erase( j ); 
                return false; 
            }
        }
        
        i = std::find( open.begin(), open.end(), p );
//        j = std::find( Open_vect_local.begin(), Open_vect_local.end(), p );
        if( i != open.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost && ( *i ).cost + ( *i ).dist > 1e-6 ) {
                //std::cout<< "  open " << ( *i ).isnode  << "  " << ( *i ).cost + ( *i ).dist << "  " << cost <<std::endl;
                return true;}
            else { 
                open.erase( i ); 
//                Open_vect_local.erase(j); 
                return false; }
        }
        return false;
    }
    
    bool existLocalPoint( point<T>& p, T cost ) {
        typename std::list<node<T>>::iterator i;
        typename std::vector<node<T>>::iterator j;
        i = std::find( closed_loc.begin(), closed_loc.end(), p );
        j = std::find( Closed_vect_local.begin(), Closed_vect_local.end(), p );
        
        //std::cout<< "  cost " << ( *i ).cost + ( *i ).dist << "  " << cost <<std::endl;
        
        if( i != closed_loc.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost && ( *i ).cost + ( *i ).dist > 1e-6 ){ 
                //std::cout<< "  closed " << ( *i ).isnode << "  " << ( *i ).cost + ( *i ).dist << "  " << cost <<std::endl;
                return true;}
            else { 
                closed_loc.erase( i ); 
                Closed_vect_local.erase( j ); 
                return false; 
            }
        }
        
        i = std::find( open_loc.begin(), open_loc.end(), p );
        j = std::find( Open_vect_local.begin(), Open_vect_local.end(), p );
        if( i != open_loc.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost && ( *i ).cost + ( *i ).dist > 1e-6 ) {
                //std::cout<< "  open " << ( *i ).isnode  << "  " << ( *i ).cost + ( *i ).dist << "  " << cost <<std::endl;
                return true;}
            else { 
                open_loc.erase( i ); 
                Open_vect_local.erase(j); 
                return false; }
        }
        return false;
    }    
 
    bool fillOpen( node<T>& n , const int &icount) {
        
        int size, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);        
        
        T stepCost, nc, dist;
        point<T> neighbour;
        int inode = n.isnode;
        int Nmax = Vertex_cost_and_NN[inode].size();
        int imin, imax;
        int exitRank;
        int send_note_loc = 0;
        int send_note_glob = 0;
        
        mpi_loop_splitter(&Nmax, &imin, &imax);

        
        //open_loc.push_back( n );
        //open_loc.pop_front();
        //Closed_vect_local.push_back(n);
        //closed_loc.push_back(n);

        //std::cout << "  Nmax " << Nmax << std::endl;
        for( int in = imin; in <imax ; in++ ) {
            
            std::vector< std::pair<T ,int> > p = Vertex_cost_and_NN[inode];
            stepCost = p[in].first;

            int nn = p[in].second;
            std::pair<double, double> w = getX_Y_at_point(nn);
           
            point<T> NN(w.first, w.second);
            
            neighbour = NN;
            //np[0].x = neighbour.x;
            //np[0].y = neighbour.y;

            //np[1].x = n.pos.x;
            //np[1].y = n.pos.y;               
            
            //std::cout<< " node " << nn << "  " << in << std::endl;
            
            send_note_loc = 0;
            if( std::sqrt(calcDist(neighbour)) < 1e-3 ) { 
                
                std::cout << "  icount " << icount << " End reached. "<< rank << std::endl;
                
               send_note_loc = 1;
               return true;
/*
                int i;
                for (i = 0; i < size; i++) {
                  if (i != rank) {
                    MPI_Send(&send_note, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                  }
                }
*/                
                
                //return true;
                
            }
            //MPI_Recv(&send_note, 1, MPI_INT, rank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //MPI_Barrier(MPI_COMM_WORLD);
            //std::cout <<  " imin "  << imin <<  " imax "  << imax << std::endl;
            //if(send_note == 1) 
            //    std::cout <<  " rank "  << rank <<  " send "  << send_note << std::endl;
            //else
            //    return true;
            
            //int *bf_note = (int *) calloc(size, sizeof(int));
            //MPI_Barrier(MPI_COMM_WORLD);
            //MPI_Allgather(&send_note, 1, MPI_INT, bf_note, 1, MPI_INT, MPI_COMM_WORLD);
            
            //for(int i=0; i < size; i++)
            //    if(bf_note[i]>=1) return true;
/*
            if (rank == exitRank) {
                // If we are the root process, send our data to everyone
                int i;
                for (i = 0; i < size; i++) {
                    if (i != rank) {
                       MPI_Send(&send_note, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                    }
                }
            } else {
                // If we are a receiver process, receive the data from the root
                MPI_Recv(&send_note, 1, MPI_INT, exitRank, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }   
  */          
            
            if( isValid( neighbour ) && m(nn) != 1) {
                
                nc = stepCost + n.cost;
                dist = calcDist( neighbour );
                dist = std::sqrt(dist);
                
                if( !existPoint( neighbour, nc + dist ) ) {
                    
                    //exitRank = rank;
                    //std::cout << "  Rank " << rank << std::endl;
                    //std::cout<< " check 1 " << exitRank <<  "  " << in << std::endl;
                    
                    //MPI_Bcast( &np[0], 2, Position, exitRank, MPI_COMM_WORLD );
                    
                    //std::cout<< " check 2 " << rank <<  "  " << in << std::endl;
                    //MPI_Bcast( &nc, 1, MPI_DOUBLE, exitRank, MPI_COMM_WORLD );
                    //std::cout<< " check 3 " << rank <<  "  " << in << std::endl;
                    //MPI_Bcast( &dist, 1, MPI_DOUBLE, exitRank, MPI_COMM_WORLD );
                    //std::cout<< " check 4 " << rank <<  "  " << in << std::endl;
                    //MPI_Bcast( &nn, 1, MPI_INT, exitRank, MPI_COMM_WORLD );
                    //MPI_Bcast( &exitRank, 1, MPI_INT, exitRank, MPI_COMM_WORLD );
                    //std::cout<< " check 5 " << rank <<  "  " << in << std::endl;
                    
                    //neighbour.x = np[0].x;
                    //neighbour.y = np[0].y;
                    
                    //n.pos.x = np[1].x;
                    //n.pos.y = np[1].y;
                   
                    node<T> m;
                    m.cost = nc; m.dist = dist;
                    m.pos = neighbour; 
                    m.parent = n.pos;
                    m.isnode = nn;
                    open.push_back( m );
                    
//                    Open_vect_local.push_back(m);
                    
//                    std::cout<< "rank "<< rank << " lists " << Closed_vect_local[Closed_vect_local.size()-1].isnode << "  " <<  Open_vect_local[Open_vect_local.size()-1].isnode  <<std::endl;
                    //std::cout <<  nc + dist << " node " << nn << "  " << neighbour << "  " << neighbour.y << std::endl;
                    //std::cout <<  stepCost << " node " << nn << "  " << n.pos.x << "  " << n.pos.y << std::endl;
                    
                }

                //MPI_Barrier(MPI_COMM_WORLD);
/*                
                for(int i = 0; i<size; i++){
                
                    if(i != exitRank){
                    
                        neighbour.x = np[0].x;
                        neighbour.y = np[0].y;
                        
                        //bool set = existPoint( neighbour, nc + dist );

                        n.pos.x = np[1].x;
                        n.pos.y = np[1].y;                        
                        
                        node<T> m;
                        m.cost = nc; m.dist = dist;
                        m.pos = neighbour; 
                        m.parent = n.pos;
                        m.isnode = nn;
                        open.push_back( m );                    
                    
                    
                    }
                    
                
                }
*/ 
                
            }

           
            
        }
        
        ///MPI_Allreduce(&send_note_loc, &send_note_glob, 1, MPI_INT,MPI_SUM,MPI_COMM_WORLD);
           
        //if(send_note_glob >=1)
        //    return true;
            //std::cout <<  " rank " << rank << " " << send_note_glob  << std::endl;        
        
        //if(rank == 1){
        int hode=-1;;
        int last;
        /*
        node<T> closed_NN; // = Closed_vect_local[Closed_vect_local.size()-1];
        std::list<node<double>>::iterator ie = closed_loc.begin();
        std::advance(ie,closed_loc.size()-1);
        closed_NN =( *ie );
        */ 
        
  
        for(typename std::list<node<double>>::iterator i = open.begin(); i != open.end(); i++ ) {
            hode++;
            
            last = ( *i ).isnode;
            //if(icount == 1)
            //    std::cout <<  " rank new = "  << open.size() <<  " node = "  << last << std::endl;
            
        }
        //if(icount == 1) exit(0);
        
        //}
       
        
        std::list<node<double>>::iterator ie = closed.begin();
        std::advance(ie,closed.size()-1);
        
        if(hode > 0){
        //    std::cout<< " Last rank new = "  << rank <<  " node = "  << ( *ie ).isnode << std::endl;
            //MPI_Barrier(MPI_COMM_WORLD);           
            
            //exit(0);
        } 

        
        if(icount > 0 && size > 1 ){
            
            //std::cout<<" check 1 " << icount << std::endl;
            //MPI_Barrier(MPI_COMM_WORLD);
        
        node<T> closed_NN;
        Node new_NN;
        
        std::list<node<double>>::iterator ir = closed.begin();
        std::advance(ir,closed.size()-1);
        closed_NN =( *ir );

        new_NN.isnode = closed_NN.isnode;
        new_NN.cost = closed_NN.cost;
        new_NN.dist = closed_NN.dist;
        new_NN.pos_x = closed_NN.pos.x;
        new_NN.pos_y = closed_NN.pos.y;
        new_NN.parent_x = closed_NN.parent.x;
        new_NN.parent_y = closed_NN.parent.y;        

        //std::cout<<" check 2 " << icount << std::endl;

        //int size_nn = closed.size();        
        Node *vector_NN;
        Node *open_loc_list;
        Node *closed_loc_list;
        Node *open_glob_list;
        Node *closed_glob_list;        
        
        int open_local_size, open_global_size;
        open_local_size = open.size();
        int closed_local_size, closed_global_size;
        closed_local_size = closed.size();
        
        vector_NN = (Node *) calloc(size, sizeof(Node));
                        
        //std::cout<<" check 3 " << icount << std::endl;

        int *buf_sizes_open = (int *) calloc(size, sizeof(int));
        int *buf_sizes_closed = (int *) calloc(size, sizeof(int));
        
        MPI_Allgather(&open_local_size, 1, MPI_INT, buf_sizes_open, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&closed_local_size, 1, MPI_INT, buf_sizes_closed, 1, MPI_INT, MPI_COMM_WORLD);
        
        //std::cout<<" check 4 " << icount << std::endl;
        
        MPI_Status status;
        MPI_Datatype mynodes;
        MPI_Datatype type[7] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
        int blocklen[7] = {1, 1, 1, 1, 1, 1, 1};
        MPI_Aint disp[7];

        disp[0] = offsetof(Node, isnode);
        disp[1] = offsetof(Node, cost);
        disp[2] = offsetof(Node, dist);        
        disp[3] = offsetof(Node, pos_x);
        disp[4] = offsetof(Node, pos_y);
        disp[5] = offsetof(Node, parent_x);
        disp[6] = offsetof(Node, parent_y);
        
        MPI_Type_create_struct(7, blocklen, disp, type, &mynodes);
        MPI_Type_commit(&mynodes);         
        
        //std::cout<<" check 5 " << icount << std::endl;
        
        MPI_Allgather(&new_NN, 1, mynodes, vector_NN, 1, mynodes, MPI_COMM_WORLD);
        //std::cout<<" check 6.1 " << icount << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        
        //std::cout<<" check 7 " << icount << std::endl;

        //int *rbuf = (int *) calloc(size, sizeof(int));//  malloc(size*sizeof(int));
        //int *displs = (int *) calloc(size, sizeof(int));
        
        //for(int i=0; i<size; i++){
            
        //    rbuf[i] = buf_sizes[i];
            //displs[i] = i*open_local_size;

        //}        
            
        //displs[0] = 0;
        //for (int i=1; i<size; i++)
        //    displs[i] = displs[i-1] + rbuf[i-1];
        
        //MPI_Allgatherv(&open_loc_list[0], open_local_size, mynodes, open_glob_list, rbuf, displs, mynodes, MPI_COMM_WORLD);
        
        //MPI_Allgather(&open_loc_list[0], open_local_size, mynodes, open_glob_list, open_local_size, mynodes, MPI_COMM_WORLD);
        //std::cout<<" check 7.1 " << icount << std::endl;
        //MPI_Barrier(MPI_COMM_WORLD);
        //MPI_Allgather(&closed_loc_list[0], closed_local_size, mynodes, closed_glob_list, closed_global_size, mynodes, MPI_COMM_WORLD);        

        //MPI_Barrier(MPI_COMM_WORLD);
        //std::cout<<" check 7 " << icount << std::endl;
        //exit(0);
              
        std::vector< std::pair<T,int> > costs;
        
        for(int i=0; i<size; i++) {
            
            T c = vector_NN[i].cost;
            std::pair<T,int> p = std::make_pair(c, i);
            
            costs.push_back(p);
            
        }
        //std::cout<<" check 6 " << icount << std::endl;
         
        std::sort(costs.begin(),costs.end());        
        
        int low_cost_proc = costs[0].second;   
        int low_cost_proc_size_open = buf_sizes_open[low_cost_proc];
        int low_cost_proc_size_closed = buf_sizes_closed[low_cost_proc];
        
        open_loc_list = (Node *) calloc(low_cost_proc_size_open, sizeof(Node));
        closed_loc_list = (Node *) calloc(low_cost_proc_size_closed, sizeof(Node));        
        
        int is;
        //std::cout<<" check 7 " << icount << std::endl;
        if(rank == low_cost_proc){
            is = 0;
            for(typename std::list<node<double>>::iterator i = closed.begin(); i != closed.end(); i++ ){

                closed_loc_list[is].isnode = ( *i ).isnode;
                closed_loc_list[is].cost = ( *i ).cost;
                closed_loc_list[is].dist = ( *i ).dist;
                closed_loc_list[is].pos_x = ( *i ).pos.x;
                closed_loc_list[is].pos_y = ( *i ).pos.y;
                closed_loc_list[is].parent_x = ( *i ).parent.x;
                closed_loc_list[is].parent_y = ( *i ).parent.y;
                is++;

            }
            
            is=0;
            for(typename std::list<node<double>>::iterator i = open.begin(); i != open.end(); i++ ){

                open_loc_list[is].isnode = ( *i ).isnode;
                open_loc_list[is].cost = ( *i ).cost;
                open_loc_list[is].dist = ( *i ).dist;
                open_loc_list[is].pos_x = ( *i ).pos.x;
                open_loc_list[is].pos_y = ( *i ).pos.y;
                open_loc_list[is].parent_x = ( *i ).parent.x;
                open_loc_list[is].parent_y = ( *i ).parent.y;
                is++;

            }                    
            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        //std::cout<<" check 8 " << icount << std::endl;
        
        MPI_Bcast(&open_loc_list[0], low_cost_proc_size_open, mynodes, low_cost_proc, MPI_COMM_WORLD);
        //std::cout<<" check 9 " << icount << std::endl;
        MPI_Bcast(&closed_loc_list[0], low_cost_proc_size_closed, mynodes, low_cost_proc, MPI_COMM_WORLD);
        //std::cout<<" check 10 " << icount << std::endl;
        
        
        open.clear();
        
        for(int it=0; it < low_cost_proc_size_open ; it++){
        
            node<T> m;
            
            m.isnode = open_loc_list[it].isnode;
            m.dist   = open_loc_list[it].dist;
            m.cost   = open_loc_list[it].cost;
            m.parent.x   = open_loc_list[it].parent_x;
            m.parent.y   = open_loc_list[it].parent_y;
            m.pos.x   = open_loc_list[it].pos_x;
            m.pos.y   = open_loc_list[it].pos_y;            
            
            open.push_back(m);
        
        }
        
        closed.clear();
        for(int it=0; it < low_cost_proc_size_closed ; it++){
        
            node<T> m;
            
            m.isnode = closed_loc_list[it].isnode;
            m.dist   = closed_loc_list[it].dist;
            m.cost   = closed_loc_list[it].cost;
            m.parent.x   = closed_loc_list[it].parent_x;
            m.parent.y   = closed_loc_list[it].parent_y;
            m.pos.x   = closed_loc_list[it].pos_x;
            m.pos.y   = closed_loc_list[it].pos_y;            
            
            closed.push_back(m);
        
        }        
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(rank == 1){
            
        for(typename std::list<node<double>>::iterator i = closed.begin(); i != closed.end(); i++ ) {
            hode++;
            
            last = ( *i ).isnode;
            //std::cout <<  " rank = "  << rank <<  " node = "  << last << std::endl;
            
        }        
               
        }
        //MPI_Barrier(MPI_COMM_WORLD);
        //        exit(0);
                
        }
        
//        MPI_Allgather(&new_NN, 1, mynodes, vector_NN, 1, mynodes, MPI_COMM_WORLD);
        
        //for(int i=0; i<size_nn; i++)
        //    std::cout<< size_nn << " " << vector_NN[i].isnode << std::endl;
        //            exit(0);  
        
        /*
        
        int irr = *(&vector_NN + 1) - vector_NN; 
        if(irr == 2 && (vector_NN[0].isnode == vector_NN[1].isnode) ){
        
            size_nn--;
            
        }
       
        
        for( int in = 0; in < size_nn ; in++ ) {
            
            std::vector< std::pair<T ,int> > p = Vertex_cost_and_NN[inode];
            int inext = vector_NN[in].isnode;
            
            stepCost = p[inext].first;

            int nn = p[inext].second;
            std::pair<double, double> w = getX_Y_at_point(nn);
           
            point<T> NN(w.first, w.second);
            neighbour = NN;
           
            
            if( std::sqrt(calcDist(neighbour)) < 1e-3 ) { 
                
                std::cout << "  Rank " << rank << " End reached." << std::endl;
                std::cout << "  Rank " << rank << " x = " << neighbour.x << " y = " << neighbour.y << std::endl;
                return true; 
            }
        
            if( isValid( neighbour ) && m(nn) != 1 ) {
                                 
                nc = stepCost + n.cost;
                dist = calcDist( neighbour );
                dist = std::sqrt(dist);
                
                if( !existPoint( neighbour, nc + dist ) ) {
                    //std::cout <<  in << " node " << std::endl;
                   
                    node<T> m;
                    m.cost = nc; m.dist = dist;
                    m.pos = neighbour; 
                    m.parent = n.pos;
                    m.isnode = nn;
                    open.push_back( m );
                    //std::cout <<  in << " node " << nn << "  " << n.pos.x << "  " << n.pos.y << std::endl;
                    
                }

                
            }
            
            
        } 
   */                
                    
/* 
        for( int x = 0; x < 8; x++ ) {
            stepCost = 1; // x < 4 ? 1 : 1; // here cost
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
 
        
        exit(0);
        
        open_loc.clear();
        closed_loc.clear();
        Closed_vect_local.clear();
        Open_vect_local.clear();
 */   
        //exit(0);
        MPI_Barrier(MPI_COMM_WORLD);


        return false;

    }
 
    bool search( int start_p, int end_p, map<T>& mp ) {
        
        int size, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);          
        
        int icount;
        std::pair<T ,T> pair_x_y;
      
        pair_x_y = graph.getX_Y_at_Vertex(start_p); 
        point<T> s(pair_x_y.first, pair_x_y.second);
         
        pair_x_y = graph.getX_Y_at_Vertex(end_p);         
        point<T> e(pair_x_y.first, pair_x_y.second);
        short_path.push_back(end_p);
       
        point<T> Zero(0.0,0.0);
        node<T> n; end = e; start = s; m = mp;

        n.cost = T(0) ; n.pos = s; n.parent = Zero; n.dist = std::sqrt(calcDist( s )); 
        n.isnode = start_p;

        //std::cout << " x = " << n.dist << std::endl;
        icount = 0;
        open.push_back( n );
        while( !open.empty() ) {
            //open.sort();
            node<T> n = open.front();
            //std::cout << " node = " << n.isnode << " dist = " << n.dist << std::endl;
            open.pop_front();
            closed.push_back( n );
            
            if( fillOpen( n, icount) ){
                //MPI_Barrier(MPI_COMM_WORLD);
                return true;}
            icount++;
            
        }
        return false;
    }
 
    int path( std::list<point<T>>& path ) {
        
        path.push_front( end );
        int cost = 1 + closed.back().cost;
        
        path.push_front( closed.back().pos );
        point<T> parent = closed.back().parent;
        //short_path.push_back( closed.back().isnode );

        for(typename std::list<node<T>>::reverse_iterator i = closed.rbegin(); i != closed.rend(); i++ ) {
            
            if( ( *i ).pos == parent && !( ( *i ).pos == start ) ) {
                
                path.push_front( ( *i ).pos );
                parent = ( *i ).parent;
                
                
            }
            
            short_path.push_back(  ( *i ).isnode );
           
        }
        
        //std::sort(short_path.begin(), short_path.end());
        path.push_front( start );
        return cost;
    }
    
    void clean_lists(){
    
        open.clear();
        closed.clear();
    
    }
 
    map<T> m; point<T> end, start;
    std::vector< point<T> > neighbours;
    std::list< node<T> > open;
    std::list< node<T> > closed;
    std::vector<int> short_path;
    
    std::list< node<T> > open_loc;
    std::list< node<T> > closed_loc;    
    std::vector< node<T> > Open_vect_local;
    std::vector< node<T> > Closed_vect_local;
    
};
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
    
    int vert_min, vert_max;
    int Nvert_loc;
    T *AdjMat_loc;
    T *d_loc;
    int *pred_loc;
    int *visited_nodes;
    
    int size, rank;
    
public:
    //default constructor
    Dijkstra();

    //constructor
    Dijkstra(Graph<T> &G){
        graph = G;        
        Nvertices = graph.getNumVertex();
        AdjMat.resize(Nvertices);
        AdjMat = graph.getAdjMat();
        
        Vertex_cost_and_NN.resize(Nvertices);
        for(int i=0; i<Nvertices; i++)
            Vertex_cost_and_NN[i].resize(graph.getVertex_cost_and_NN()[i].size());
        
        Vertex_cost_and_NN = graph.getVertex_cost_and_NN();
        
        d.resize(Nvertices);
        pred.resize(Nvertices);
        
        for(int u=0; u<Nvertices; u++){
            d[u] = infinity;
            pred[u] = -1;
        }

    }
    
    void My_MPI_World(int mysize, int myrank){
    
        size = mysize; rank = myrank;
        
    
    }    

    void mpi_loop_splitter(const int Nloop, int &loopmin, int &loopmax){

        div_t divider;
        int quot,rem;

        divider = div(Nloop,size);

        quot = divider.quot;
        rem = divider.rem;

        if(rank < divider.rem){

            loopmin = rank * divider.quot  + rank;
            loopmax = loopmin + divider.quot + 1 ;        

        }else{

            loopmin = rank * divider.quot  + divider.rem;
            loopmax = loopmin + divider.quot;          

        }

    }

    
    void Find_path_Dijkstra(const int &source, const int &target, std::vector<int> &path){
    
    
        d[source] = 0;
        pred[source] = -1;
        std::set< std::pair<T, int> > Q_queue;
        Q_queue.insert ( std::make_pair(d[source], source) );
        
        while(!Q_queue.empty()){
        
            typename std::set<std::pair<T, int> >::iterator iter = Q_queue.begin();
            int u = ( *iter ).second;
            Q_queue.erase( Q_queue.begin());
            
            
            int ns = Nvertices;
            //for(int ic = 0; ic<Nvertices; ic++){
                
            //    if(AdjMat[u * Nvertices + ic] < infinity)
            //        ns++;
            
            //}
            //std::cout << ns << std::endl;
            //ns = Vertex_cost_and_NN[u].size();
            for(int v=0; v< ns; v++ ){
                
                T weight = AdjMat[u * Nvertices + v]; //Vertex_cost_and_NN[u][v].first; //
                int ver = v;  //Vertex_cost_and_NN[u][v].second;
                
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
    void printPath(std::vector<int> parent, int j) 
    { 

        // Base Case : If j is source 
        if (parent[j] == -1)
            return; 
        
        int i = parent[j];
        //std::cout << i << " " << parent[i] << std::endl;
        printPath(parent, i); 

        printf("%3d", j); 
    } 

    // A utility function to print  
    // the constructed distance 
    // array 
    void printSolution(std::vector<T> & dist, int src, std::vector<int>  parent) 
    { 
        //int src = 1; 
        printf("Vertex\t Distance\tPath"); 
        for (int i = 0; i < Nvertices; i++) 
        { 
            if(i != src){
                //std::cout << dist[i] << " " << parent[i] << std::endl;
                printf("\n%d -> %d \t\t %3g\t\t%d ", 
                              src, i, dist[i], src); 
                printPath(parent, i); 
            }
        } 
    }
    
void Print_paths(std::vector<int> pred ) {
    int v, w, *path, count, i;
    int n = Nvertices;

      path =  (int *) malloc(Nvertices*sizeof(int));
//      dist =  (T *) malloc(Nvertices*sizeof(T));

      printf("The shortest path from 0 to each vertex is:\n");
      printf("  v   distance   Path 0->v\n");
      printf("----    ---------     ---------\n");
      for (v = 1; v < n; v++) {
         printf("%d:    ", v);
         count = 0;
         w = v;
         while (w != -1) {
            path[count] = w;
            count++;
            w = pred[w];
           
         }
         printf("%3f    ", d[v]); 
         printf("0 ");
         for (i = count-1; i >= 0; i--){
            printf("%d ", path[i]);
         }   
         printf("\n");
      }

      free(path);


} 
    void Allocate_local_AdjMat_d_pred(){
    
        mpi_loop_splitter(Nvertices, vert_min, vert_max);
                
        Nvert_loc =  vert_max-vert_min;
        AdjMat_loc = (T *) malloc(Nvertices*Nvert_loc*sizeof(T));
        d_loc = (T *) malloc(Nvert_loc*sizeof(T));
        pred_loc = (int *) malloc(Nvert_loc*sizeof(int));
        visited_nodes = (int *) malloc(Nvert_loc*sizeof(int));
        
    }
    
    void Deallocate_local_AdjMat_d_pred(){
                
        free(AdjMat_loc);
        free(d_loc);
        free(pred_loc);
        free(visited_nodes);
        
    }    
    
    int Pos_Vertices_in_Global(int pos_loc, int n_pos_loc, int rank) {
        return pos_loc + vert_min; //rank*n_pos_loc;
    }     

    void Split_AdjMat(int source){
 
        
        for (int i=0; i<Nvertices; i++){
            for (int j_loc=0; j_loc<Nvert_loc; j_loc++){
                int j = Pos_Vertices_in_Global(j_loc,Nvert_loc,rank);
                AdjMat_loc[i*Nvert_loc + j_loc] = AdjMat[i*Nvertices + j];
                
            }
        }    
        
        for (int j_loc=0; j_loc<Nvert_loc; j_loc++){
            //int j = Pos_Vertices_in_Global(j_loc,Nvert_loc,rank);
            d_loc[j_loc] = AdjMat_loc[source*Nvert_loc + j_loc];
            pred_loc[j_loc] = -1;
            visited_nodes[j_loc] = 0;
            
        }
        
        if(source >= vert_min && source < vert_max){
            visited_nodes[source - vert_min] = 1;
        }
    }
    
    void Dijkstra_search(const int &source, const int &target){
    
        Allocate_local_AdjMat_d_pred();
        Split_AdjMat(source);
        Find_path_Dijkstra_in_parallel(source, target);
        
    }
    
    int find_closest_locally_vertex(int source){
        
        std::vector<std::pair<T, int> > p;
        for(int i=0; i<Nvert_loc; i++){
            if(!visited_nodes[i])
                p.push_back( std::make_pair( d_loc[i] , i) );
        }    
        std::sort(p.begin(),p.end());
        
        return p[0].second;
       
    
    } 
    
    int Find_min_dist(double loc_dist[], int loc_known[], int loc_n, int my_rank, 
        MPI_Comm comm) {
        int loc_v, loc_u;
        T loc_min_dist = T(1E6);

        loc_u = 1;
        for (loc_v = 0; loc_v < loc_n; loc_v++)
            if (!loc_known[loc_v])
                if (loc_dist[loc_v] < loc_min_dist) {
                    loc_u = loc_v;
                    loc_min_dist = loc_dist[loc_v];
                }

        return loc_u;
    }    

    void Find_path_Dijkstra_in_parallel(const int &source, const int &target){
       
        /*
        struct double_int{

            double dist;
            int vertex;
    
        } ver_dist_loc, ver_dist_global;        
        */
        
        int *loc_min_pos = (int *) malloc(size*sizeof(int));//  new int[size];
        double  *loc_min_dist_for_pos = (double *) malloc(size*sizeof(double)); // new double[size];
        int *glob_min_pos = (int *) malloc(size*sizeof(int)); // new int[size];
        double  *glob_min_dist_for_pos = (double *) malloc(size*sizeof(double)); // new double[size];        
        
        for (int iu=1; iu<Nvertices; iu++){
            
            //if(rank == 0){
            //    std::cout<< std::endl;
            //    std::cout << " LOOP: " << iu << std::endl;
            //}    
       
            //int u_low_dist_loc = find_closest_locally_vertex(source);
            int u_low_dist_loc = Find_min_dist(d_loc, visited_nodes, Nvert_loc, rank, MPI_COMM_WORLD);
            
            
            //std::cout << u_low_dist_loc << " loc_min_dist_for_pos  " << d_loc[u_low_dist_loc] << std::endl;
            double dist_u_loc = d_loc[u_low_dist_loc];
            int j = Pos_Vertices_in_Global(u_low_dist_loc,Nvert_loc, rank);
            //ver_dist_loc.dist = dist_u_loc; ver_dist_loc.vertex = j;
            
            for(int is=0; is<size; is++){
                glob_min_pos[is]=0; 
                glob_min_dist_for_pos[is]=0.0;
                if(is==rank){
                
                    loc_min_pos[is] = j;
                    loc_min_dist_for_pos[is] = dist_u_loc;
                    
                    //std::cout << "  loc_min_pos " << loc_min_pos[is] << "   loc_min_dist_for_pos  " << dist_u_loc << std::endl;
                
                }else{
                
                    loc_min_pos[is] = 0;
                    loc_min_dist_for_pos[is] = 0.0;                
                
                }
            }
            
            //MPI_Status  status;
            //MPI_Request request_pos;
            //MPI_Request request_dist;
            
            int ret_val ;
            MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
            ret_val = MPI_Allreduce( &loc_min_pos[0], &glob_min_pos[0], size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
           
            
            if (ret_val == MPI_ERR_COMM){
                std::cout << "MPI_Gather : MPI_ERR_COMM"<< std::endl;
                exit(0);
            }else if (ret_val == MPI_ERR_COUNT){
                std::cout << "MPI_Gather : MPI_ERR_COUNT"<< std::endl;
                exit(0);
            }else if (ret_val == MPI_ERR_TYPE){
                std::cout << "MPI_Gather : MPI_ERR_TYPE"<< std::endl;
                exit(0);
            }else if (ret_val == MPI_ERR_BUFFER){
                std::cout << "MPI_Gather : MPI_ERR_BUFFER"<< std::endl;
                exit(0);
            }            
            
            //MPI_Iallreduce( &loc_min_pos[0], &glob_min_pos[0], size, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &request_pos);
            //MPI_Wait(&request_pos, &status);
            
            MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
            ret_val = MPI_Allreduce( &loc_min_dist_for_pos[0], &glob_min_dist_for_pos[0], size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
           
            /*
            if (ret_val == MPI_ERR_COMM){
                std::cout << "MPI_Gather : MPI_ERR_COMM"<< std::endl;
                exit(0);
            }else if (ret_val == MPI_ERR_COUNT){
                std::cout << "MPI_Gather : MPI_ERR_COUNT"<< std::endl;
                exit(0);
            }else if (ret_val == MPI_ERR_TYPE){
                std::cout << "MPI_Gather : MPI_ERR_TYPE"<< std::endl;
                exit(0);
            }else if (ret_val == MPI_ERR_BUFFER){
                std::cout << "MPI_Gather : MPI_ERR_BUFFER"<< std::endl;
                exit(0);
            }
            */                
            //MPI_Iallreduce( &loc_min_dist_for_pos[0], &glob_min_dist_for_pos[0], size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request_dist);
            //MPI_Wait(&request_dist, &status);
            
//            std::cout<<" check 1 " << count1 << " " << count2 << std::endl;
            //MPI_Wait(&request_dist, MPI_STATUS_IGNORE);
            
            //if(check != 0) exit(0);
            
            std::vector<std::pair<T,int> > p;
            for(int i=0; i<size; i++){
                
                    p.push_back( std::make_pair(glob_min_dist_for_pos[i] , glob_min_pos[i]) );
            }    
            std::sort(p.begin(),p.end());
            
 
            //if(rank == 0){
            
                //for(int is=0; is<size; is++)
                //    std::cout << "glob_min_pos  " << glob_min_pos[is] << "   glob_min_dist_for_pos  " << glob_min_dist_for_pos[is] << std::endl;          
            
           // }
            
            //if(rank == 0){
            //    std::cout<< std::endl;
            //    std::cout << " END LOOP: " << iu << std::endl;            
            //}
  
            //if(rank == 0)
                //std::cout << " check 1 " << rank<<  " " <<  j << "  " <<ver_dist_loc.dist << "  " << u_low_dist_loc << "  " << iu<< std::endl;            
            //MPI_Barrier(MPI_COMM_WORLD);    
            //MPI_Allreduce( &ver_dist_loc, &ver_dist_global, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
             //if(rank == 0)
               // std::cout << " check 2 " << rank<<  " " << ver_dist_global.vertex << "  " << iu<< std::endl;   
            //std::cout << " check 2 "<< rank << "  " << ver_dist_global[0].vertex<<  " " << ver_dist_global[0].dist<< std::endl;  
            //MPI_Barrier(MPI_COMM_WORLD);
             //if(rank == 0)
                //std::cout << " check 3 " << rank<<  " " << iu<< std::endl;   
            //MPI_Barrier(MPI_COMM_WORLD);
            
            int u = p[0].second;//ver_dist_global.vertex;
            T glob_min_u_dist = p[0].first;//ver_dist_global.dist;
            
            if(u >= vert_min && u < vert_max){
                 //std::cout << " check 1 " << rank <<  " " << u - vert_min << "   " << iu<< std::endl;
                 if(u - vert_min < 0) exit(0);
                visited_nodes[u - vert_min] = 1;
                //std::cout << " check 2 " << rank <<  " " << visited_nodes[u - vert_min] << "   " << iu<< std::endl;
            }    
                    
            
            int ns = Nvert_loc;
            for(int v=0; v< ns; v++ ){
                
                if(!visited_nodes[v] ){
                
                    
                    T weight = AdjMat_loc[u * Nvert_loc + v]; 
                    int ver = v; 
                    
                    if( glob_min_u_dist + weight < d_loc[v] ){

                        
                        d_loc[ver] = glob_min_u_dist  + weight;
                        pred_loc[ver] = u;
                        //std::cout << " check 1 " << pred_loc[ver] << " " << u<< std::endl;

                    }
                }
                
                
                
            }
        }

        //MPI_Barrier(MPI_COMM_WORLD);
        
        free(loc_min_pos);
        free(loc_min_dist_for_pos);
        free(glob_min_pos);
        free(glob_min_dist_for_pos);
        
        //MPI_Barrier(MPI_COMM_WORLD);
        //std::cout<< " going out " <<std::endl;
        //exit(0);
        
        Gather_dist_and_pred_global();
        Deallocate_local_AdjMat_d_pred();

        //exit(0);
        //if(rank == 0)
        //    Print_paths(pred);
            //printSolution(d, source, pred);

    }
    
    void Gather_dist_and_pred_global(){
        
        T *d_glob = (T *) calloc(Nvertices, sizeof(T)); // malloc(Nvertices*Nvertices*sizeof(T));
        int *pred_glob = (int *) calloc(Nvertices, sizeof(int)); // malloc(Nvertices*Nvertices*sizeof(int));  
                
        int *buf_sizes = (int *) calloc(size, sizeof(int));
        
        MPI_Allgather(&Nvert_loc, 1, MPI_INT, buf_sizes, 1, MPI_INT, MPI_COMM_WORLD);        
        
        int *rbuf = (int *) calloc(size, sizeof(int));
        int *displs = (int *) calloc(size, sizeof(int));
        
        for(int i=0; i<size; i++){
            
            rbuf[i] = buf_sizes[i];
            displs[i] = i*Nvert_loc;

        }        
        
        //MPI_Gather( &d_loc[0], Nvert_loc, MPI_DOUBLE, d_glob, Nvert_loc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Gather( &pred_loc[0], Nvert_loc, MPI_INT, pred_glob, Nvert_loc, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Allgatherv( &pred_loc[0], Nvert_loc, MPI_INT, pred_glob, rbuf, displs, MPI_INT, MPI_COMM_WORLD);        
        MPI_Allgatherv( &d_loc[0], Nvert_loc, MPI_DOUBLE, d_glob, rbuf, displs, MPI_DOUBLE, MPI_COMM_WORLD);
                
        d.clear();
        //pred.clear();
        //if(rank == 0)
            for(int i= 0; i<Nvertices; i++){
            
                d.push_back(d_glob[i]);
                pred[i] = pred_glob[i];
                //std::cout << i << " " << visited_nodes_glob[i]  << std::endl;
            
            }
    
        free(d_glob);
        free(pred_glob);
        free(buf_sizes);
        free(rbuf);
        free(displs);
    }

    
};



#endif
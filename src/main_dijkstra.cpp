/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main_A_Start.cpp
 * Author: hamilton
 *
 * Created on October 24, 2019, 10:18 AM
 */

#include <mpi.h>
#include <chrono>  // for high_resolution_clock
#include <dijkstra.hpp>

//void mpi_loop_splitter(int *size, int *loopmin, int *loopmax);

int main( int argc, char* argv[] ) {
   
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    
    
    readInputFile<double> inFile;
   
    char **f =  argv;
    std::vector< std::vector<double> > MAP = inFile.readFileToVector(f[1]);
    
    bool read_file = true;
    Graph<double> graph(0, read_file, MAP);
    
    Dijkstra<double> dijkstra(graph);
    
    int is = 0;
    int ie = 10;

    double time_avg = 0.0;
    int icount = 0;
        
    do {
    
         std::vector<int> path;
        auto start = std::chrono::high_resolution_clock::now();
        
        dijkstra.Find_path_Dijkstra(is, ie, path);
        
        auto finish = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = finish - start;
        time_avg = time_avg + elapsed.count();
        icount++;

    } while (icount < 1);
    time_avg = time_avg / (double) (icount);
    
    //if(rank == 1)
    std::cout << std::endl << "Time wasted = " << time_avg << " Seconds"  <<std::endl;
     
    MPI_Finalize();
    return 0;
}





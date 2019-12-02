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
#include <A_Start.h>

//void mpi_loop_splitter(int *size, int *loopmin, int *loopmax);

int main( int argc, char* argv[] ) {
   
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    
    
    readInputFile<double> inFile;
   
    char **f =  argv;
    std::vector< std::vector<double> > MAP = inFile.readFileToVector(f[1]);
    
    //int index = 2;
    //inFile.displayVector(MAP,index);
    
    bool read_file = true;
    Graph<double> graph(0, read_file, MAP);
    
    
    //std::pair<double, double> v = graph.getX_Y_at_Vertex(1);
    //std::cout<< v.first << "  " << v.second << std::endl;
 
   /* 
    std::vector< std::vector< std::pair<double ,int> > > Vertex_cost_and_NN;
    Vertex_cost_and_NN = graph.getVertex_cost_and_NN();
    
    for (int i=0; i< Vertex_cost_and_NN[1].size(); i++)
        std::cout << "(" << 1 <<  "," << Vertex_cost_and_NN[1][i].second << ")" << Vertex_cost_and_NN[1][i].first << std::endl;
    

    std::vector< std::pair<double ,double> > Vertex_X_Y;
    Vertex_X_Y = graph.getVertex_X_Y();    
    for (int i=0; i< Vertex_X_Y.size(); i++)
        std::cout << "(" << Vertex_X_Y[i].first <<  "," << Vertex_X_Y[i].second << ")" << std::endl;    
    
   std::vector<double> AdjMat = graph.getAdjMat();
    int nver = graph.getNumVertex();
    
    for (int i=0; i<  nver; i++){
        for (int j=0; j<  nver; j++){
            
            std::cout << "(" << i << "," << j << ") :" ;
            std::cout << AdjMat[i*nver + j] << " ";
            std::cout << std::endl;
            
        }
        
    } 
    exit(0);
    */ 

    
    //aStar<double> A(graph);
    //std::pair<double, double> w = A.getX_Y_at_point(1);
    //std::cout<< w.first << "  " << w.second << std::endl;
    
    map<double> m(graph);
    //point<double> s(0, 0), e( 7, 7 );
    
    aStar<double> as(graph);
    
    int is = 0;
    int ie = 2;

    double time_avg = 0.0;
    int icount = 0;
        
    do {
    
        auto start = std::chrono::high_resolution_clock::now();

//        for (int is = imin; is < imax; is++ ){

//            for (int ie = 0; ie < NNodes; ie++ ){

                if( as.search( is, ie, m ) ) {

//                    std::pair<double, double> v = graph.getX_Y_at_Vertex(ie);

                    //std::cout << "Check: " << std::endl;
                    //std::cout << " x = " << v.first << " y = " << v.second << std::endl;

                    //std::list<point<double>> path;
                    //int c = as.path( path );        
                    //std::cout<< "rank =" << rank <<std::endl;
                    //for(int i=0;  i < as.short_path.size(); i++ ) {
                    //    std::cout<< "(" << as.short_path[i]  <<  ") ";
                   //}
                    as.clean_lists();
                    
                }else{
                    
                    std::cout<< "  Rank " << rank << " Path from is = "<< is << " to " << " ie = "<< ie << " not found"<< std::endl;
                    std::cout << "\n\n";            

//                }else{

//                    std::cout << "Path from is = "<< is << " to " << " ie = "<< ie << " not found"<< std::endl;

                }

//            }

//        }

        auto finish = std::chrono::high_resolution_clock::now();
        
//        MPI_Barrier(MPI_COMM_WORLD);
        std::chrono::duration<double> elapsed = finish - start;
        time_avg = time_avg + elapsed.count();
        icount++;

    } while (icount < 1);
    time_avg = time_avg / (double) (icount);
    
    //if(rank == 1)
    std::cout << "Time wasted = " << time_avg << " Seconds"  <<std::endl;
    MPI_Abort(MPI_COMM_WORLD, 911);
    
    
 /*
    if( as.search( s, e, m ) ) {
        std::list<point<double>> path;
        int c = as.path( path );
        for( int y = -1; y < 9; y++ ) {
            for( int x = -1; x < 9; x++ ) {
                
                if( x < 0 || y < 0 || x > 7 || y > 7 || m( x, y ) == 1 )
                    std::cout << char(0xdb);
                else {
                    if( std::find( path.begin(), path.end(), point<double>( x, y ) ) != path.end() )
                        std::cout << "x";
                    else std::cout << ".";
                }
            }
            std::cout << "\n";
        }
 
        std::cout << "\nPath cost " << c << ": ";
        for(typename std::list<point<double>>::iterator i = path.begin(); i != path.end(); i++ ) {
            std::cout<< "(" << ( *i ).x << ", " << ( *i ).y << ") ";
        }
    }
    std::cout << "\n\n";
*/    
    MPI_Finalize();
    return 0;
}



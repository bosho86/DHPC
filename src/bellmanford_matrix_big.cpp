/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bellmanford_matrix_big.cpp
 * Author: Paulina
 *
 * Created on October 15, 2019, 10:18 AM
 */


#include <iostream>
#include <stdio.h>
#include <bits/stdc++.h>
#include <cmath>
#include <vector>
#include <sys/time.h>
#include <mpi.h>
#include "read_input_file.h"
#include "Graph.hpp"


using namespace std;

void graphmatrix(int V, int matrix[5][5]){
    int *mat;
	mat = (int *) malloc(V * V * sizeof(int));
	        for (int i = 0; i < V; i++){
	            for (int j = 0; j < V; j++) {
	                mat[V*i+j] = matrix[i][j];
	            }
	            cout << mat[i] << endl;
	        }
    cout << mat;

}

int convert_dimension_2D_1D(int x, int y, int n) {
       return x * n + y;
   }



void bellmanford(double *mat, double *distance, int V, int src, int my_rank, int p, bool *has_negative_cycle){
int local_V;
int local_start, local_finish;
double *local_matrix;
double *local_distance;

//MPI_Bcast(variable, size, MPI_INIT, 0, comunicator)
if (my_rank == 0) {
       local_V = V;
   }

MPI_Bcast(&local_V, 1, MPI_INT, 0, MPI_COMM_WORLD);
//Split
local_start = (local_V / p) * my_rank ; //my rank=0, empieza en 0
local_finish = (local_V / p) * (my_rank + 1);

if(my_rank == p - 1){
   local_finish = local_V;
}

//allocate the memory
local_matrix = (double *) malloc(local_V * local_V * sizeof(double));
local_distance = (double *) malloc(local_V * sizeof(double));

// broadcast matrix
if (my_rank == 0){
    memcpy(local_matrix, mat, sizeof(int) * local_V * local_V);
}

MPI_Bcast(local_matrix, local_V * local_V, MPI_INT, 0, MPI_COMM_WORLD);


///WHY I HAVE A SEGMENTATION FAUL WITH THE LOCAL DISTANCE ???
for(int ii = 0; ii < local_V; ii++){
    local_distance[ii] = INT_MAX;
}

local_distance[0] = 0;




bool local_change;
int local_iteration = 0;

for (int iter = 0; iter < local_V-1; iter++){
	local_change = false;
	local_iteration ++;

	for (int u = local_start; u < local_finish; u++){
		for(int v = 0; v < local_V; v++){
			int weight = local_matrix[convert_dimension_2D_1D(u, v, local_V)];
			
			if(local_distance[u] + weight < local_distance[v]){
				local_distance[v] = local_distance[u] + weight;
			
				local_change = true;

			}

		}

	}
//MPI reduction
//int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
//Collective communication can occur ``in place'' for intracommunicators, with the output buffer being identical to
//the input buffer. This is specified by providing a special argument value, MPI_IN_PLACE, instead of the send buffer or the receive buffer argument.
//MPI_LOR..return the logical
MPI_Allreduce(MPI_IN_PLACE, &local_change, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
if(!local_change)
	break;
MPI_Allreduce(MPI_IN_PLACE, local_distance, local_V, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
}

//checking for the negative cycle
if(local_iteration = local_V -1){
	local_change = false;

	for (int u = local_start; u < local_finish; u++){
			for(int v = 0; v < local_V; v++){
				int weight = local_matrix[convert_dimension_2D_1D(u, v, local_V)];
				
				if(local_distance[u] + weight < local_distance[v]){
					local_distance[v] = local_distance[u] + weight;

					local_change = true;
					break;

				}

			}

		}
MPI_Allreduce(&local_change, has_negative_cycle, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
}
//check if the results go back
if (my_rank = 0){
	memcpy(distance, local_distance, local_V * sizeof(int));
}


free(local_matrix);
free(local_distance);

}




int main(int argc, char **argv){
	int V;
    // definir graph
	//std::vector<double> graph;
	std::vector<double> AdjMat;
	double *mat;
	double *distance;
	bool has_negative_cycle = false;

	readInputFile<double> inFile;

	    char **f =  argv;
	    std::vector< std::vector<double> > MAP = inFile.readFileToVector(f[1]);

	    bool read_file = true;
	    Graph<double> graph(0, read_file, MAP);
   // la matriz es doble.
         V = graph.getNumVertex();
         AdjMat.resize(V * V);
	     AdjMat = graph.getAdjMat();


    mat = (double *) malloc(V * V * sizeof(double));
    for (int i = 0; i < V; i++){
        for (int j = 0; j < V; j++) {
             mat[V*i+j] = AdjMat[V*i+j];
        }
     }


    distance = (double *) malloc( V * sizeof(double));


    // Initialize MPI

    MPI_Init(&argc, &argv);

    int p; // number of processes
    int my_rank; // ranks

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    double t1, t2;
   
    t1 = MPI_Wtime();

    
    bellmanford(mat, distance, V, 0, my_rank, p, &has_negative_cycle);

    t2 = MPI_Wtime();

    if(my_rank == 0){
      std::cerr.setf(std::ios::fixed);
      std::cerr << std::setprecision(6) << "Time(s): " << (t2 - t1) << endl;
    }


   MPI_Finalize();


   return 0;

}

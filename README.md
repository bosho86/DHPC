# DHPC
Authors: H. Carrillo Nunez-Meder, D.P. Aguirre-Hirt

The aim of this project is to develop and implement a parallelization scheme, based on the Message-Passing-Interface
(MPI), of two different algorithms to find the shortest path, given a source node, to all other nodes in a graph. The well-
known Dijkstra’s [1] and Bellman-Ford’s [2, 3] algorithms are considered. Our findings suggest that the former requires more communication, leading the latter to a better scaling when the number of processors increases.

Code Organization
1. Workstation.ipynb -> Map or Graph generator.(ie. map.3000.dat -> Adjancent matrix of 3000 vertices)
2. Include:
A_Star.h -> Try to paralellize the A*algo.
disktra.hpp- >spliting fuctions of the adjacent matrix.
read_input_file.h -> helpter to read the map file.
3. Src:
1. bellmanford_matrix_big.cpp -> parallelized version of bellman Ford
2. main_dijkstra.cpp -> parallelized version of the Dijkstra Algortihm
3. main_A_Start.cpp -> parallelized attempt of the A* Algorithm

Please go to Report.pdf, to see the performance and detailed explanation of this code.



References

[1] E. W. Dijkstra, “A note on two problems in connexion with graphs,” Numerische Mathematik, vol. 1, pp. 269–271, 1959.
[2] R. Bellman, “On a routing problem,” Quarterly of Applied Mathematics, vol. 16, pp. 87–90, 1958.
[3] L. R. Ford, “Network flow theory,” Rand Corporation,p. 923, 1956.

# DHPC
Authors: H. Carrillo Nunez-Meder, D.P. Aguirre-Hirt

The aim of this project is to develop and implement a parallelization scheme, based on the Message-Passing-Interface
(MPI), of two different algorithms to find the shortest path, given a source node, to all other nodes in a graph. The well-
known Dijkstra’s [1] and Bellman-Ford’s [2, 3] algorithms are considered. Our findings suggest that the former requires more communication, leading the latter to a better scaling when the number of processors increases.


References

[1] E. W. Dijkstra, “A note on two problems in connexion with graphs,” Numerische Mathematik, vol. 1, pp. 269–271, 1959.
[2] R. Bellman, “On a routing problem,” Quarterly of Applied Mathematics, vol. 16, pp. 87–90, 1958.
[3] L. R. Ford, “Network flow theory,” Rand Corporation,p. 923, 1956.
[4] N. Edmonds, A. Breuer, D. Gregor, and A. Lumsdaine,“Single-source shortest paths with the parallel boost graph library,” .
[5] X. Han, Q. Sun, and J. Fan, “Parallel dijkstra’s algorithm based on multi-core and mpi,” Applied Mechanics and Materials, vol. 441, pp. 750–753, 2013.
[6] V. T. Chakaravarthy, F. Checconi, F. Petrini, and Y. Sabharwal, “Scalable single source shortest path algo-
rithms for massively parallel systems,” in 2014 IEEE 28th International Parallel and Distributed Processing
Symposium, May 2014, pp. 889–901.
[7] K. Nikas, N. Anastopoulos, G. Goumas, and N. Koziris, “Employing transactional memory and helper threads
to speedup dijkstra’s algorithm,” in In ICPP, 2009, pp.
388–395.
[8] A Grama, A. Gupta, G. Karypis, and V. Kumar, Introduction to Parallel Computing, chapter 10, p. 429, Pearson, 2 edition, 2003.
[9] A. A. Hagberg, D. A. Schult, and P. J. Swart, “Exploring network structure, dynamics, and function us-
ing networkx,” in Proceedings of the 7th Python in Science Conference (SciPy2008), August 2008, pp. 11–15.

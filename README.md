# Description
A graph $G$ is quasiplanar if it admits a drawing in the plane where no 3 edges are pairwise intersecting.
This code enumerates all simple quasiplane drawings of a graph.

Based on the code from "On maximal 3-planar graphs" by Michael Hoffmann, Meghana M. Reddy, and Shengzhe Wang (EuroCG 2024), which has been adapted to enumerate quasiplanar graphs 
and $k$-planar graphs (for arbitrary $k$).

 - The file `hds_kplanar.h` contains the code to enumerate $k$-planar graphs for any choice of $k\geq 0$ (generalization of code from paper, which works for $k \leq 3$).
 - The file `hds_quasiplanar.h` contains the code to enumerate quasiplanar graphs with local crossing number $k\geq0$ (parameter).

The local crossing number (`klim`) controls the depth of the recurrence: the higher `klim` the more branching the algorithm do.

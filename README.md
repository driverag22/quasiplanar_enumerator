# Description
A graph $G$ is quasiplanar if it admits a drawing in the plane where no 3 edges are pairwise intersecting.
This code enumerates all simple quasiplane drawings of a graph up to strong isomorphism.

Based on the code from "On maximal 3-planar graphs" by Michael Hoffmann, Meghana M. Reddy, and Shengzhe Wang (EuroCG 2024), which has been adapted to enumerate quasiplanar graphs 
and $k$-planar graphs (for arbitrary $k$).

 - The file `hds_kplanar.h` contains the code to enumerate $k$-planar graphs for any choice of $k\geq 0$ (generalization of code from paper, which works for $k \leq 3$).
 - The file `hds_quasiplanar.h` contains the code to enumerate quasiplanar graphs with local crossing number $k\geq0$ (parameter).

The local crossing number (`klim`) controls the depth of the recurrence: the higher `klim` the more branching the algorithm do.

The current implementation has a hard-coded limit of 96 edges (via `inline static constexpr std::size_t MAX_EDGES = 96`).

# JSON drawing representation

Created drawings can be saved into a JSON format that describes the drawing. 
From it, the drawing can be recreated by hand, or reloaded in the code as a `Drawing` struct.

The function `hds_quasiplanar/serialize_to_json()` generates the JSON object using the `nlohmann` library, which can then be saved as a 
JSON file using `nlohmann::json/dump()`.

When constructing a new `Drawing` struct, one can pass a `nlohmann::json` object, and the code will then follow the recipe specified by the json file to make the drawing.

## JSON format

First the abstract graph is described, giving the list of edges, described by the two endpoints:

```json
"abstract_graph": [ [0,1], [0,2], ...]
```

The edges are ordered lexicographically.

Then the drawing recipe is given. It describes how to draw the edges in order, meaning that one can just follow the list from top to bottom 
to draw the graph by hand.

For the first edge it contains the edge label and the two endpoints `u,v`. The first edge is identified by having label `0` and being the first in the recipe.

For subsequent edges it contains a unique label, the labels of edges it crosses, the endpoints `u` and `v`, and label of the edge it "starts after" in the counter-clockwise rotation around `u`.
Note that therefore the start after edge is also incident on `u`.

This drawing_recipe follows the same order as the code itself when constructing a new drawing.

The `kplane` and `num_vertices` of the drawing is also given. 
Note that `kplane` is not necessarily the local crossing number of the drawing, but the `klim` that the code was subject to when it constructed the drawing (so the local crossing number is at most `kplane`, but could be less).

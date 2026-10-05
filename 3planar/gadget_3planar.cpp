#include "hds_kplanar.h"
#include "iso.h"

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const Edges edges = {
    // K8
    {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},
    {1,2},{1,3},{1,4},{1,5},{1,6},{1,7},
    {2,3},{2,4},{2,5},{2,6},{2,7},
    {3,4},{3,5},{3,6},{3,7},
    {4,5},{4,6},{4,7},
    {5,6},{5,7},
    {6,7},
    // degree 5 vertex
    {0,8},{1,8},{2,8},{3,8},{4,8},
    // rest of vertices
    {0,9},{2,9},{5,9},
    {1,10},{3,10},{6,10},

    {2,11},{5,11},{9,11},
    {3,12},{6,12},{10,12},

    // vertex 13 used to verify that (0,1) is isolated from the rest of the vertices
    // by running assert(!d.is_drawing_extensible())
    {0,13},{1,13}, 
};

const std::size_t n = 14;
const std::size_t klim = 3;
const std::string src = "extension_8_rel";
const std::string dir = "gadget_drawings_full";

int main() {
    std::vector< Drawing<klim> > solutions;
    std::size_t cnt = 0;
    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);

    const auto start_edge = edges.begin() + 1;
    for (auto e = start_edge;;) {
        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];

        HdsPath p = d.first_path(u, v);
        if (p.empty()) {
BACKUP:
            // no way to add uv -> do previous edges differently
            do {
                if (e == start_edge) {
                    goto END;
                }
                --e;

                u = (*e)[0];
                assert(u == d.edges.back().u);
                v = (*e)[1];
                assert(v == d.edges.back().v);
                p = d.edges.back().built;
                d.remove_edge();
            } while (!d.next_path(p, v));
        }
        d.add_edge(p, v);

        if (++e == edges.end()) {
            bool newSol = true;
            for (auto it = solutions.begin(); it != solutions.end(); it++) {
                if(are_isomorphic((*it),d)) {
                    newSol = false;
                    break;
                }
            }
            if (newSol) {
                assert(!d.is_drawing_extensible());
                std::cout << cnt++ << std::endl;
                solutions.push_back(d);
                for (const auto& e : d.edges) {
                    if (e.u == 0 && e.v == 1) assert(e.ncr == 0);
                }
            }
            goto BACKUP;
        }
    }
END:

    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    // because vertex 13 can be in either of the two faces incident to (0,1)
    // we expect 2 * 1444 = 2888 solutions
    assert(solutions.size() == 2888);
    return 0;
}

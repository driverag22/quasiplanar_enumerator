#include "hds_kplanar.h"
#include "iso.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12;
const std::size_t klim = 3;
const std::string split = "3";

const Edges edges = {
    // base K6
    {0,1},{0,2},{0,3},{0,4},{0,5},
    {1,2},{1,3},{1,4},{1,5},
    {2,3},{2,4},{2,5},
    {3,4},{3,5},
    {4,5},

    // connection to outer 6-vertex-graph
    {0,6},{0,7},{0,11},
    {1,6},{1,7},{1,8},{1,9},{1,11},
    {2,7},{2,8},{2,9},
    {3,8},{3,9},{3,10},{3,11},
    {4,9},{4,10},{4,11},
    {5,6},{5,10},{5,11},

    // outer 6-vertex-graph
    {6,7},{6,8},{6,10},{6,11},
    {7,8},{7,9},{7,11},
    {8,9},{8,10},
    {9,10},{9,11},
    {10,11},
};

std::vector<std::pair<std::size_t, std::size_t>> get_missing_edges(const Drawing<klim>& d, std::size_t num_vertices) {
    std::vector<std::vector<bool>> adj(num_vertices, std::vector<bool>(num_vertices, false));
    for (const auto& edge : d.edges) {
        adj[edge.u][edge.v] = true;
        adj[edge.v][edge.u] = true;
    }

    std::vector<std::pair<std::size_t, std::size_t>> missing;
    for (std::size_t u = 0; u < num_vertices; ++u)
        for (std::size_t v = u + 1; v < num_vertices; ++v)
            if (!adj[u][v]) missing.push_back({u, v});
    return missing;
}

bool is_drawing_extendable(const Drawing<klim>& d, std::size_t num_vertices) {
    auto missing_edges = get_missing_edges(d, num_vertices);

    for (const auto& [u, v] : missing_edges) {
        Drawing<klim> d_search(d);
        HdsPath p = d_search.first_path(u, v);

        if (!p.empty()) {
            std::cout << "\n  [!] Edge (" << u << ", " << v << ") can be legally added!";
            return true;
        }
    }
    return false;
}

int main() {
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(10000,0); // assume no more than 10000 unique drawings up to iso

    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);
    std::size_t num_fixed_edges = 5;
    for (std::size_t i = 2; i <= num_fixed_edges; ++i) {
        HdsPath p = d.first_path(0, i);
        if (p.empty()) {
            throw std::runtime_error("Failed to build the initial star!");
        }
        d.add_edge(p, i);
    }
    std::cout << "3crossing: split " << split << " | Star built" << std::endl;

    auto start_edge = edges.begin() + num_fixed_edges;
    std::size_t counter = 0;

    for (auto e = start_edge;;) {
        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];

        HdsPath p = d.first_path(u, v);
        if (p.empty()) {
BACKUP:
            // no way to add uv -> do previous edges differently
            do {
                if (e == start_edge) {
                    goto NEXT_PERM;
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
            std::size_t d_ind = 0;
            for (auto it = solutions.begin(); it != solutions.end(); it++) {
                if(are_isomorphic((*it),d)) {
                    newSol = false;
                    d_cnt[d_ind]++;
                    break;
                }
                d_ind++;
            }
            if (newSol) {
                solutions.push_back(d);
                d_cnt[d_ind] = 1;
                std::cout << "Solution found: " << counter++ << "\n\n";
                if (is_drawing_extendable(d,n)) {
                    std::cout << "Extendable drawing!!!: " << d_ind << std::endl;
                }
            }
            goto BACKUP;
        }
    }
NEXT_PERM:
std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
if(solutions.size() == 0) return 0;

std::size_t idx = 0;
for (auto it = solutions.begin();it!=solutions.end();it++) {
    std::string filename = "../quasiDrawings/K12_3planar/3cross/split" + split + "/" + std::to_string(idx) + ".json";
    std::ofstream of_json(filename);
    nlohmann::ordered_json output_json = (*it).serialize_to_json();
    of_json << output_json.dump(4);
    of_json.close();

    std::string filename2 = "../quasiDrawings/K12_3planar/3cross/split" + split + "/" + std::to_string(idx) + ".graphml";
    std::ofstream of_graphml(filename2);
    (*it).graphml_output(of_graphml);
    of_graphml.close();
    idx++;
}

std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;
for (std::size_t i = 0; i < solutions.size(); i++) {
    std::cout << "Drawing-" << i << " appeared " << d_cnt[i] << " times.\n\n";
}
return 0;
}

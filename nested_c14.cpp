// #include "hds_quasiplanar.h"
#include "hds_kplanar.h"
#include "iso.h"
//#include <sstream>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

Edges generate3C14Graph() {
    Edges edges;

    // 1. inner C14 cycle (Vertices 0..13) - FIXED FIRST EDGES
    for (std::size_t i = 0; i < 14; ++i)
        edges.push_back({i, (i + 1) % 14});

    // 1b. inner C14 chords at distance 2
    // for (std::size_t i = 0; i < 14; ++i)
    //     edges.push_back({i, (i + 2) % 14});

    // 2. inner to middle connections at distance 3
    for (std::size_t i = 0; i < 14; ++i) {
        if (i % 2 == 0) edges.push_back({i, 14 + ((i + 3) % 14)});
        else            edges.push_back({i, 14 + ((i + 11) % 14)});
    }

    // 3. middle C14 cycle (Vertices 14..27)
    for (std::size_t i = 0; i < 14; ++i)
        edges.push_back({14 + i, 14 + ((i + 1) % 14)});

    // 4. middle to outer connections at distance 3
    for (std::size_t i = 0; i < 14; ++i) {
        if (i % 2 == 0) edges.push_back({i, 14 + ((i + 3) % 14)});
        else if (i % 2 == 1) edges.push_back({i, 14 + ((i + 11) % 14)});
    }

    // 5. outer C14 cycle (Vertices 28..41)
    for (std::size_t i = 0; i < 14; ++i)
        edges.push_back({28 + i, 28 + ((i + 1) % 14)});

    // 5b. Outer C14 chords at distance 2
    // for (std::size_t i = 0; i < 14; ++i)
    //     edges.push_back({28 + i, 28 + ((i + 2) % 14)});

    return edges;
}

const std::size_t n = 28;
const std::size_t klim = 3;

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
        if (!p.empty()) return true;
    }

    return false;
}


int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n << std::endl;
    const Edges edges = generate3C14Graph();
    // const std::size_t minimal_cr = 0;
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(10000,0); // assume no more than 10000 unique drawings up to iso

    Drawing<klim> d(n);
    // Edge 0: (0,1)
    d.add_first_edge(edges[0][0], edges[0][1]);
    // Edges 1..12: (1,2), (2,3), ..., (12,13)
    for (std::size_t i = 1; i < 13; ++i) {
        std::size_t u = edges[i][0];
        std::size_t v = edges[i][1];
        HdsPath p = d.first_path(u, v);
        if (p.empty()) {
            throw std::runtime_error("Failed to build the initial C14 cycle path!");
        }
        d.add_edge(p, v);
    }
    // Edge 13: Closing the first C14 cycle (13, 0)
    {
        std::size_t u = edges[13][0];
        std::size_t v = edges[13][1];
        HdsPath p = d.first_path(u, v);
        if (p.empty()) {
            throw std::runtime_error("Failed to close the initial C14 cycle!");
        }
        d.add_edge(p, v);
    }

    auto start_edge = edges.begin() + 14;

    int counter = 0;
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
                // if (is_drawing_extendable(d,n)) std::cout << "Not maximal!" << std::endl;
                // minimal_cr = std::min(minimal_cr, d.crossings.size());
                if (++counter % 25 == 0) std::cout << counter << std::endl;
            }
            goto BACKUP;
        }
    }
END:
    std::cout << "Found " << solutions.size() << " min crossing drawings in total." << std::endl;
    if(solutions.size() == 0) {
        return 0;
    }

    // std::size_t idx = 0;
    // for (auto it = solutions.begin();it!=solutions.end();it++) {
    //     std::string filename = "../drawingsQuasi/K9/minGlobalCr/jsons/" + std::to_string(idx) + ".json";
    //     std::ofstream of_json(filename);
    //     nlohmann::ordered_json output_json = (*it).serialize_to_json();
    //     of_json << output_json.dump(4);
    //     of_json.close();

    //     // std::string filename2 = "drawings/K7_prop_test/minCr_" + std::to_string(idx) + ".graphml";
    //     // std::ofstream of_graphml(filename2);
    //     // (*it).graphml_output(of_graphml);
    //     // of_graphml.close();
    //     idx++;
    // }

    std::cout << "Found " << counter << " drawings in total." << std::endl;
    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;
    // std::cout << "Found " << idx << " min crossing drawings in total." << std::endl;
    // std::cout << "Minimal Crossing Number is "<<minimal_cr<<std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;
}

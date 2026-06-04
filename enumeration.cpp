#include "hds_quasiplanar.h"
//#include "hds_kplanar.h"
#include "iso.h"
//#include <sstream>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

Edges generateCompleteGraph(std::size_t n) {
    Edges edges;
    edges.reserve(n * (n - 1) / 2);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            edges.push_back({i, j});
        }
    }
    return edges;
}

const std::size_t n = 7; // note hard-coded limit of 64 edges for quasiplanar...
//const std::string split = "3t1i";
// const Edges edges =
// {
//     {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},{0,8},{0,9},{0,10},
//     {1,2}, {1,3},{1,4},{1,5},{1,6},{1,7},{1,8},{1,9},
//     {2,3},{2,4},{2,5},{2,6},{2,7},{2,8},{2,9},
//     {3,4},{3,5},{3,6},{3,7},{3,8},{3,9},
//     {4,5},{4,6},{4,7},{4,8},{4,9},
//     {5,6},{5,7},{5,8},{5,9},
//     {6,7},{6,8},{6,9},
//     {7,8},{7,9},
//     {8,9},
// };

const std::size_t klim = 9;

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n << std::endl;
    const Edges edges = generateCompleteGraph(n);
    std::size_t minimal_cr = 0x3f3f3f3f;
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(10000,0); // assume no more than 10000 unique drawings up to iso

    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);
    std::size_t num_fixed_edges = n;
    for (std::size_t i = 2; i < num_fixed_edges; ++i) {
        HdsPath p = d.first_path(0, i);
        if (p.empty()) {
            throw std::runtime_error("Failed to build the initial star!");
        }
        d.add_edge(p, i);
    }
    std::cout << "Star built" << std::endl;

    auto start_edge = edges.begin() + (num_fixed_edges-1);
    // auto start_edge = edges.begin() + 1;

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
                minimal_cr = std::min(minimal_cr,d.crossings.size());

                std::cout << ++counter << std::endl;
            }
            goto BACKUP;
        }
    }
END:
    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) {
        return 0;
    }

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "drawings/K7_k9/jsons/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        // std::string filename2 = "drawings/K7_prop_test/minCr_" + std::to_string(idx) + ".graphml";
        // std::ofstream of_graphml(filename2);
        // (*it).graphml_output(of_graphml);
        // of_graphml.close();
        idx++;
    }

    std::cout << "Found " << counter << " drawings in total." << std::endl;
    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;
    std::cout << "Found " << idx << " min crossing drawings in total." << std::endl;
    std::cout << "Minimal Crossing Number is "<<minimal_cr<<std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;
}

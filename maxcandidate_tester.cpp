#include "hds_quasiplanar.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12;      // Total vertices
const std::size_t klim = 17;   // 2n - 7 = 17
const std::string split = "subgraph_c_0_to_5_block";

struct EdgeSetWithMissing {
    Edges edges;
    std::vector<std::pair<std::size_t, std::size_t>> missing_edges;
};

std::vector<EdgeSetWithMissing> generate_edge_sets() {
    std::vector<EdgeSetWithMissing> all_edge_sets;
    for (std::size_t c = 0; c < 6; ++c) {
        std::vector<std::size_t> S; // S = {0..5} \ {c}
        for (std::size_t u = 0; u < 6; ++u)
            if (u != c) S.push_back(u);

        // iterate over non-neighbor (nn) of c in S
        for (std::size_t nn : S) {
            EdgeSetWithMissing item;

            // edges inside {6..9}
            for (std::size_t u = 6; u < 10; ++u) {
                for (std::size_t v = u + 1; v < 10; ++v) {
                    item.edges.push_back({u, v});
                }
            }
            // connections from {6..9} to new 10 (old c)
            for (std::size_t u = 6; u < 10; ++u)
                item.edges.push_back({u, 10});

            // connections from {6..9} to S
            for (std::size_t u = 6; u < 10; ++u)
                for (std::size_t s : S)
                    item.edges.push_back({s, u});

            // connections from new 10 (old c) to S \ {nn}
            for (std::size_t s : S) if (s != nn) 
                    item.edges.push_back({s, 10});

            // connection between c (old 10) and new 10
            item.edges.push_back({c, 10});

            // connections for 11: connects to all of S and new 10
            for (std::size_t s : S) item.edges.push_back({s, 11});
            item.edges.push_back({10, 11});

            std::sort(item.edges.begin(), item.edges.end());

            // edge (nn, 10) omitted between new 10 and non-neighbor nn
            item.missing_edges.push_back({nn, 10});

            // edge (c, 11) omitted between old 10 (now c) and 11
            item.missing_edges.push_back({c, 11});

            // edges between {6..9} and, c (old 10) and 11
            for (std::size_t u = 6; u < 10; ++u) {
                item.missing_edges.push_back({c, u}); item.missing_edges.push_back({u, 11});
            }

            all_edge_sets.push_back(item);
        }
    }

    return all_edge_sets;
}

bool is_drawing_extendable(const Drawing<klim>& d, 
        const std::vector<std::pair<std::size_t, std::size_t>>& missingEdges) {
    for (const auto& [u, v] : missingEdges) {
        Drawing<klim> d_search(d);
        HdsPath p = d_search.first_path(u, v);

        while (!p.empty()) {
            Drawing<klim> d_test(d);
            d_test.add_edge(p, v);
            if (d_test.verify_quasiplanarity()) {
                std::cout << "\n  [!] Edge (" << u << ", " << v << ") can be legally added!";
                return true;
            } else {
                std::cout << "not quasi!\n";
            }
            if (!d_search.next_path(p, v)) break;
        }
    }
    return false;
}

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "max candidate, k = " << klim << ", n = " << n << ", (proper) split = " << split << std::endl;

    std::vector<EdgeSetWithMissing> all_edge_sets = generate_edge_sets();
    std::size_t total_edge_sets = all_edge_sets.size();
    std::cout << "Generated " << total_edge_sets << " edge set configurations." << std::endl;

    // Iterate through all drawings of K_6
    int idx = 0;
    for (int i = 0; i < 63; i++) {
        std::cout << "Drawing " << std::to_string(i) << std::endl;
        for (const auto& item : all_edge_sets) {
            idx++;
            std::cout << "\rPermutation [" << idx << " / " << total_edge_sets << "]" << std::flush;
            std::ifstream input_file("../quasiDrawings/K6/jsons/" + std::to_string(i) + ".json");
            if (!input_file.is_open()) {
                std::cerr << "Could not open drawing file " << i << std::endl;
                continue;
            }

            nlohmann::json import_data; 
            input_file >> import_data; 
            input_file.close();

            Drawing<klim> base_d(import_data, n);


            Drawing<klim> d = base_d;
            const Edges& current_edges = item.edges;

            auto start_edge = current_edges.begin();

            for (auto e = start_edge;;) {
                std::size_t u = (*e)[0];
                std::size_t v = (*e)[1];
                HdsPath p = d.first_path(u, v);

                if (p.empty()) {
                    // Backtrack
                    do {
                        if (e == start_edge) {
                            goto NEXT_ITEM; // Fail this permutation, try the next edge set
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

                if (++e == current_edges.end()) {
                    std::cout << "\nFound solution for drawing " << i << "!" << std::endl;
                    assert(d.verify_quasiplanarity());
                    if (is_drawing_extendable(d, item.missing_edges)) {
                        std::cout << "extendable\n";
                        return 0;
                    }
                    goto NEXT_ITEM;
                }
            }
NEXT_ITEM:;
        }
        std::cout << "\nPerm " << std::to_string(idx) << " finished." << std::endl;
    }

    std::cout << "Found no extensible solutions with k = " << klim << " for split = " << split << std::endl;
    return 0;
}

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
const std::string split = "1";

struct EdgeSetWithMissing {
    Edges edges;
    std::vector<std::pair<std::size_t, std::size_t>> missing_edges;
};

std::vector<EdgeSetWithMissing> generate_edge_sets() {
    std::vector<EdgeSetWithMissing> all_edge_sets;

    // base internal K_6 edges on V2 = {6, 7, 8, 9, 10, 11}
    Edges base_v2_edges;
    for (std::size_t u = 6; u < 12; ++u) for (std::size_t v = u + 1; v < 12; ++v)
            base_v2_edges.push_back({u, v});

    // select s1 (source vertex connecting to target t1)
    for (std::size_t s1 = 0; s1 < 6; ++s1) {
        // 2. select s2,s3 (sources connecting to target t2)
        for (std::size_t s2 = 0; s2 < 6; ++s2) {
            if (s2 == s1) continue;
            for (std::size_t s3 = s2 + 1; s3 < 6; ++s3) {
                if (s3 == s1) continue;

                // select target t1 (1 non-neigh)
                for (std::size_t t1 = 6; t1 < 12; ++t1) {
                    if (t1 == 6 + s1) continue; // matching edge

                    // target t2 (2 non-neighs)
                    for (std::size_t t2 = 6; t2 < 12; ++t2) {
                        if (t2 == t1) continue; // distinct targets
                        if (t2 == 6 + s2 || t2 == 6 + s3) continue;

                        EdgeSetWithMissing item;
                        item.edges = base_v2_edges;

                        // add bipartite edges with boolean checks
                        for (std::size_t u = 0; u < 6; ++u) {
                            for (std::size_t v = 6; v < 12; ++v) {
                                // skip matching edge
                                if (v == 6 + u) continue;
                                // skip 3 extra missing edges
                                if (u == s1 && v == t1) continue;
                                if (u == s2 && v == t2) continue;
                                if (u == s3 && v == t2) continue;

                                item.edges.push_back({u, v});
                            }
                        }

                        // Populate the 9 missing edges directly
                        for (std::size_t u = 0; u < 6; ++u) {
                            item.missing_edges.push_back({u, 6 + u});
                        }
                        item.missing_edges.push_back({s1, t1});
                        item.missing_edges.push_back({s2, t2});
                        item.missing_edges.push_back({s3, t2});

                        std::sort(item.edges.begin(), item.edges.end());
                        all_edge_sets.push_back(item);
                    }
                }
            }
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
    std::cout << "k = " << klim << ", n = " << n << ", split = " << split << std::endl;

    std::vector<EdgeSetWithMissing> all_edge_sets = generate_edge_sets();
    std::size_t total_edge_sets = all_edge_sets.size();
    std::cout << "Generated " << total_edge_sets << " edge set configurations." << std::endl;

    // Iterate through all 63 drawings of K_6
    for (int i = 0; i < 63; i++) {
        std::cout << "Drawing " << std::to_string(i) << std::endl;

        std::ifstream input_file("../quasiDrawings/K6/jsons/" + std::to_string(i) + ".json");
        if (!input_file.is_open()) {
            std::cerr << "Could not open drawing file " << i << std::endl;
            continue;
        }

        nlohmann::json import_data; 
        input_file >> import_data; 
        input_file.close();

        Drawing<klim> base_d(import_data, n);

        int idx = 0;
        for (const auto& item : all_edge_sets) {
            idx++;
            std::cout << "\rPermutation [" << idx << " / " << total_edge_sets << "]" << std::flush;

            Drawing<klim> d = base_d;
            const Edges& current_edges = item.edges;
            const auto& missing_edges = item.missing_edges;

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
                    if (!d.verify_quasiplanarity()) std::cout << "not quasi?\n";
                    if (is_drawing_extendable(d, missing_edges)) {
                        std::cout << "extendable\n";
                        return 0;
                    }
                    goto NEXT_ITEM;
                }
            }
NEXT_ITEM:;
        }
        std::cout << "\nDrawing " << std::to_string(i) << " finished." << std::endl;
    }

    std::cout << "Found no extendable solutions with k = " << klim << " for split = " << split << std::endl;
    return 0;
}

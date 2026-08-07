#include "hds_quasiplanar.h"
#include <fstream>
#include <vector>
#include <set>
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

    // 1. Base internal K_6 edges on V2 = {6, 7, 8, 9, 10, 11}
    Edges base_v2_edges;
    for (std::size_t u = 6; u < 12; ++u) for (std::size_t v = u + 1; v < 12; ++v)
            base_v2_edges.push_back({u, v});

    // 2. Select 3 distinct source vertices S from V1 = {0..5}
    std::vector<int> mask_s(6, 0);
    std::fill(mask_s.end() - 3, mask_s.end(), 1);

    do {
        std::vector<std::size_t> S;
        for (std::size_t i = 0; i < 6; ++i) {
            if (mask_s[i]) S.push_back(i);
        }

        // Select 3 distinct target vertices T from V2 = {6..11}
        std::vector<int> mask_t(6, 0);
        std::fill(mask_t.end() - 3, mask_t.end(), 1);

        do {
            std::vector<std::size_t> T_base;
            for (std::size_t j = 0; j < 6; ++j) {
                if (mask_t[j]) T_base.push_back(6 + j);
            }

            // Permute T to test all pairings (S[k] -> T[k])
            std::vector<std::size_t> T = T_base;
            do {
                // Ensure no extra missing edge overlaps with the matching edge (S[k], 6 + S[k])
                bool valid = true;
                for (std::size_t k = 0; k < 3; ++k) {
                    if (T[k] == 6 + S[k]) {
                        valid = false;
                        break;
                    }
                }
                if (!valid) continue;

                EdgeSetWithMissing item;
                item.edges = base_v2_edges;

                std::set<std::pair<std::size_t, std::size_t>> missing_set;

                // Add 6 matching missing edges across V1 and V2
                for (std::size_t i = 0; i < 6; ++i)
                    missing_set.insert({i, 6 + i});

                // Add 3 extra missing edges
                for (std::size_t k = 0; k < 3; ++k)
                    missing_set.insert({S[k], T[k]});

                // Add all non-missing bipartite edges between V1 and V2
                for (std::size_t u = 0; u < 6; ++u) for (std::size_t v = 6; v < 12; ++v)
                        if (missing_set.find({u, v}) == missing_set.end()) 
                            item.edges.push_back({u, v});

                for (const auto& edge : missing_set) {
                    item.missing_edges.push_back(edge);
                }

                std::sort(item.edges.begin(), item.edges.end());
                all_edge_sets.push_back(item);

            } while (std::next_permutation(T.begin(), T.end()));
        } while (std::next_permutation(mask_t.begin(), mask_t.end()));
    } while (std::next_permutation(mask_s.begin(), mask_s.end()));

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

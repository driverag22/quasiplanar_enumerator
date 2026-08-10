#include "hds_quasiplanar.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12; // note hard-coded limit of 64 edges for quasiplanar...
const std::size_t klim = 17; // 2n-7 = 17
const std::string split = "1b";

struct EdgeSetWithMissing {
    Edges edges;
    std::vector<std::pair<std::size_t, std::size_t>> missing_edges;
};

std::vector<EdgeSetWithMissing> generate_edge_sets() {
    const bool r_masks[4][5] = {
        {1, 0, 1, 0, 0}, // ra: {7, 9}
        {1, 1, 1, 1, 0}, // rb: {7, 8, 9, 10}
        {1, 1, 0, 1, 1}, // rc: {7, 8, 10, 11}
        {0, 1, 0, 0, 1}  // rd: {8, 11}
    };
    std::vector<EdgeSetWithMissing> all_edge_sets;

    // 7 choose 3 mask to select the 3 vertices connected to all {7,8,9,10,11}
    std::vector<int> mask(7, 0);
    std::fill(mask.end() - 3, mask.end(), 1);

    do {
        Edges base_edges = {
        {7,8},{7,9},{7,10},{7,11},
        {8,9},{8,10},{8,11},
        {9,10},{9,11},
        {10,11}
        };
        std::vector<std::size_t> C; // 3 vertices connected to 7,8,9,10,11
        std::vector<std::size_t> R; // 4 remaining vertices
        for (std::size_t i = 0; i < 7; ++i) {
            if (mask[i]) for (std::size_t j=7;j<12;++j) base_edges.push_back({i,j});
            else R.push_back(i);
        }

        do {
            EdgeSetWithMissing item;
            item.edges = base_edges;

            for (std::size_t j = 0; j < 5; ++j) {
                std::size_t target = 7 + j;
                for (std::size_t k = 0; k < 4; ++k) {
                    if (r_masks[k][j]) item.edges.push_back({R[k], target});
                    else item.missing_edges.push_back({R[k], target});
                }
            }

            std::sort(item.edges.begin(), item.edges.end());
            all_edge_sets.push_back(item);

        } while (std::next_permutation(R.begin(), R.end()));
    } while (std::next_permutation(mask.begin(), mask.end()));
    return all_edge_sets;
}

bool is_drawing_extendable(const Drawing<klim>& d, 
        const std::vector<std::pair<std::size_t, std::size_t>>& missingEdges) {
    for (const auto& [u,v]: missingEdges) {
        Drawing<klim> d_search(d);
        HdsPath p = d_search.first_path(u,v);

        while (!p.empty()) {
            Drawing<klim> d_test(d);
            d_test.add_edge(p, v);
            if (d_test.verify_quasiplanarity()) {
                std::cout << "\n  [!] Edge (" << u << ", " << v << ") can be legally added!";
                return true;
            } else std::cout << "not quasi!\n";
            if (!d_search.next_path(p, v)) break;
        }
    }
    return false;
}

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n  << ", split = " << split << std::endl;

    std::vector<EdgeSetWithMissing> all_edge_sets = generate_edge_sets();
    std::size_t total_edge_sets = all_edge_sets.size();
    for (int i = 0; i < 1607; i++) {
        std::cout << "Drawing " << std::to_string(i) << std::endl;

        std::ifstream input_file("../quasiDrawings/K7/klim9_jsons/" + std::to_string(i) + ".json");
        nlohmann::json import_data; input_file >> import_data; input_file.close();
        // loading drawing
        Drawing<klim> d(import_data, n);

        int idx = 0;
        for (const auto& item : all_edge_sets) {
            std::cout << "\rPermutation [" << idx << " / " << total_edge_sets << "]" << std::flush;
            const Edges& current_edges = item.edges;
            const auto& missing_edges = item.missing_edges;

            auto start_edge = current_edges.begin();

            for (auto e = start_edge;;) {
                std::size_t u = (*e)[0];
                std::size_t v = (*e)[1];
                HdsPath p = d.first_path(u, v);
                if (p.empty()) {
                    // no way to add uv -> do previous edges differently
                    do {
                        if (e == start_edge) {
                            goto NEXT_JSON;
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
                    std::cout << "Found sol for: " << split << std::endl;
                    if (!d.verify_quasiplanarity()) std::cout << "not quasi?\n";
                    if (is_drawing_extendable(d, missing_edges)) {
                        std::cout << "extendable\n";
                        return 0;
                    }
                }
            }
        }
        std::cout << "Drawing "<< std::to_string(i) << " finished." << std::endl;
NEXT_JSON:;
    }
    std::cout << "Found no extendable solutions with k = " << klim << " for split = " << split << std::endl;
    return 0;
}

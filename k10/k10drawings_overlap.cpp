#include "hds_quasiplanar.h"
#include "iso.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <nlohmann/json.hpp>

constexpr int KPLANE_LIMIT = 15;
constexpr int NUM_CLASSES = 9;
constexpr int TOTAL_EDGES = 45;

// Get recipe according to mask in_S
nlohmann::json filter_recipe_by_subset(const nlohmann::json& full_json, const std::vector<bool>& in_S) {
    nlohmann::json filtered = full_json;
    nlohmann::json new_recipe = nlohmann::json::array();

    // full spatial rotation system at each vertex (0..9)
    std::vector<std::vector<int>> full_rot(10);
    // edges currently present in subdrawing S at each vertex
    std::vector<std::vector<int>> rot_in_S(10);

    if (full_json.contains("drawing_recipe")) {
        for (const auto& step : full_json["drawing_recipe"]) {
            int e = step["edge_label"].get<int>();
            int u = step["u"].get<int>();
            int v = step["v"].get<int>();

            // --- 1. Update Full Spatial Rotation System at u ---
            int insert_pos = full_rot[u].size();
            if (step.contains("start_after_edge") && !step["start_after_edge"].is_null()) {
                int anchor = step["start_after_edge"].get<int>();
                auto it = std::find(full_rot[u].begin(), full_rot[u].end(), anchor);
                if (it != full_rot[u].end()) {
                    insert_pos = std::distance(full_rot[u].begin(), it) + 1;
                }
            }
            full_rot[u].insert(full_rot[u].begin() + insert_pos, e);
            full_rot[v].push_back(e);

            // --- 2. Process Step if e is in Subset S ---
            if (e < static_cast<int>(in_S.size()) && in_S[e]) {
                nlohmann::json step_copy = step;

                // Filter crossed_edges
                nlohmann::json filtered_cr = nlohmann::json::array();
                if (step.contains("crossed_edges")) {
                    for (int cr : step["crossed_edges"]) {
                        if (cr < static_cast<int>(in_S.size()) && in_S[cr]) {
                            filtered_cr.push_back(cr);
                        }
                    }
                }
                step_copy["crossed_edges"] = filtered_cr;

                // --- 3. Remap Anchor based on Spatial Rotation System ---
                if (step_copy.contains("start_after_edge") && !step_copy["start_after_edge"].is_null()) {
                    int anchor = step_copy["start_after_edge"].get<int>();
                    const auto& u_s_list = rot_in_S[u];
                    bool anchor_in_S = (std::find(u_s_list.begin(), u_s_list.end(), anchor) != u_s_list.end());

                    if (!anchor_in_S) {
                        // Walk backward in full spatial rotation system from anchor
                        auto anchor_it = std::find(full_rot[u].begin(), full_rot[u].end(), anchor);
                        int remapped_anchor = -1;

                        if (anchor_it != full_rot[u].end()) {
                            // Search backwards cyclically for nearest edge in S
                            int idx = std::distance(full_rot[u].begin(), anchor_it);
                            int n = full_rot[u].size();
                            for (int step_back = 1; step_back < n; ++step_back) {
                                int prev_idx = (idx - step_back + n) % n;
                                int prev_e = full_rot[u][prev_idx];
                                if (std::find(u_s_list.begin(), u_s_list.end(), prev_e) != u_s_list.end()) {
                                    remapped_anchor = prev_e;
                                    break;
                                }
                            }
                        }

                        if (remapped_anchor != -1) {
                            step_copy["start_after_edge"] = remapped_anchor;
                        } else {
                            step_copy.erase("start_after_edge");
                        }
                    }
                }

                // Record insertion in subdrawing S
                rot_in_S[u].push_back(e);
                rot_in_S[v].push_back(e);

                new_recipe.push_back(step_copy);
            }
        }
        filtered["drawing_recipe"] = new_recipe;
        filtered["num_edges"] = new_recipe.size();
    }
    return filtered;
}

// Check strong isomorphism across all classes for subset S
bool check_global_strong_isomorphism(
    const std::vector<nlohmann::json>& recipes,
    const std::vector<bool>& in_S
) {
    try { // in case we try to draw "impossible" drawing
        // since transitive, no need to try all pairs
        nlohmann::json f0 = filter_recipe_by_subset(recipes[0], in_S);
        Drawing<KPLANE_LIMIT> d0(f0);

        for (std::size_t i = 1; i < recipes.size(); ++i) {
            nlohmann::json fi = filter_recipe_by_subset(recipes[i], in_S);
            Drawing<KPLANE_LIMIT> di(fi);

            if (!are_isomorphic(d0, di)) return false;
        }
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

// compute prefix divergence pepth via LIFO edge popping
int compute_prefix_divergence_depth(const std::vector<nlohmann::json>& recipes) {
    int min_depth = TOTAL_EDGES;

    for (std::size_t i = 0; i < recipes.size(); ++i) {
        for (std::size_t j = i + 1; j < recipes.size(); ++j) {
            Drawing<KPLANE_LIMIT> d_i(recipes[i]);
            Drawing<KPLANE_LIMIT> d_j(recipes[j]);

            int depth = TOTAL_EDGES;
            for (std::size_t k = TOTAL_EDGES; k > 0; --k) {
                // remove edges until tey are isomorphic
                if (are_isomorphic(d_i, d_j)) {
                    depth = static_cast<int>(k);
                    break;
                }
                d_i.remove_edge(); d_j.remove_edge();
            }

            if (depth < min_depth) min_depth = depth;
        }
    }
    return min_depth;
}

// MSICS backtracking solver (Branch-and-Bound)
struct MSICSSolver {
    const std::vector<nlohmann::json>& recipes;
    int prefix_depth;
    std::vector<int> candidates; // Late edges that individually pass
    std::vector<bool> current_S;
    std::vector<bool> best_S;
    int max_size = 0;
    uint64_t nodes_explored = 0;

    MSICSSolver(const std::vector<nlohmann::json>& r, int pd)
        : recipes(r), prefix_depth(pd), current_S(TOTAL_EDGES, false) {
        
        for (int e = 0; e < prefix_depth; ++e) current_S[e] = true;

        // filter candidates: only consider "late-stage" edges that individually work with prefix
        for (int e = prefix_depth; e < TOTAL_EDGES; ++e) {
            current_S[e] = true;
            if (check_global_strong_isomorphism(recipes, current_S)) candidates.push_back(e);
            current_S[e] = false;
        }

        max_size = prefix_depth;
        best_S = current_S;
    }

    void solve_dfs(int idx, int current_size) {
        nodes_explored++;

        // Early pruning: cannot beat max_size even if all remaining candidates are taken
        if (current_size + (static_cast<int>(candidates.size()) - idx) <= max_size) return;

        // Base case: evaluated all candidates
        if (idx == static_cast<int>(candidates.size())) {
            if (current_size > max_size) {
                max_size = current_size;
                best_S = current_S;
                std::cout << "  [+] Found new larger MSICS subset of size " << max_size << " edges.\n";
            }
            return;
        }

        int e = candidates[idx];

        // 1. include edge e
        current_S[e] = true;
        if (check_global_strong_isomorphism(recipes, current_S)) solve_dfs(idx + 1, current_size + 1);

        // 2: exclude edge e
        current_S[e] = false;
        solve_dfs(idx + 1, current_size);
    }
};

int main() {
    std::string base_path = "../quasiDrawings/K10_all_quasi/";
    std::vector<nlohmann::json> recipes(NUM_CLASSES);

    for (int i = 0; i < NUM_CLASSES; ++i) {
        std::string filename = base_path + std::to_string(i) + ".json";
        std::ifstream in(filename);
        in >> recipes[i];
    }

    int prefix_depth = compute_prefix_divergence_depth(recipes);
    std::cout << "Prefix Core Depth: " << prefix_depth << " edges\n";

    MSICSSolver solver(recipes, prefix_depth);
    std::cout << "Viable Late-Stage Candidates: " << solver.candidates.size() << " / " << (TOTAL_EDGES - prefix_depth) << " edges\n";

    solver.solve_dfs(0, prefix_depth);

    std::vector<int> msics_subset, prefix_core, late_core;

    // extract sol from solver
    for (int e = 0; e < TOTAL_EDGES; ++e) {
        if (solver.best_S[e]) {
            msics_subset.push_back(e);
            if (e < prefix_depth) prefix_core.push_back(e);
            else late_core.push_back(e);
        }
    }

    std::cout << "Max Strongly Isomorphic Subset Size (|S|): " << msics_subset.size() 
              << " / " << TOTAL_EDGES << " edges (" 
              << std::fixed << std::setprecision(1)
              << (static_cast<double>(msics_subset.size()) / TOTAL_EDGES * 100.0) << "%)\n\n";

    std::cout << "[1] Prefix Core (Shared prior to first divergence, e < " << prefix_depth << "): "
              << prefix_core.size() << " edges\n    ";
    for (std::size_t idx = 0; idx < prefix_core.size(); ++idx) {
        std::cout << "e" << prefix_core[idx] << (idx + 1 < prefix_core.size() ? ", " : "");
    }
    std::cout << "\n\n";

    std::cout << "[2] Late-Stage Invariant Core (Edges added AFTER divergence, e >= " << prefix_depth << "): "
              << late_core.size() << " edges\n    ";
    if (late_core.empty()) {
        std::cout << "None (Proven that no combination of late edges can extend the prefix core without breaking strong isomorphism)\n";
    } else {
        for (std::size_t idx = 0; idx < late_core.size(); ++idx) {
            std::cout << "e" << late_core[idx] << (idx + 1 < late_core.size() ? ", " : "");
        }
        std::cout << "\n";
    }

    // Final Verification
    bool verified = check_global_strong_isomorphism(recipes, solver.best_S);

    std::cout << "\n-----------------------------------------------------\n";
    std::cout << "[Verification] Subset S strongly isomorphic across all 9 classes: "
              << (verified ? "TRUE" : "FALSE") << "\n";
    std::cout << "-----------------------------------------------------\n\n";

    // std::cout << "Exported MSICS output to 'msics_result.json'\n";
    // std::cout << "Constructing HDS Drawing from filtered Class 0 recipe...\n";
    // nlohmann::json msics_filtered_json = filter_recipe_by_subset(recipes[0], solver.best_S);

    // Instantiate Drawing from filtered subdrawing recipe
    // Drawing<KPLANE_LIMIT> msics_drawing(msics_filtered_json);

    // // Export using native graphml_output method in hds_quasiplanar
    // std::ofstream graphml_out("dataAnalysis/msics_subdrawing.graphml");
    // msics_drawing.graphml_output(graphml_out);
    // std::cout << "Exported native GraphML to 'msics_subdrawing.graphml' (open in yEd)\n";

    // // Export JSON summary
    // nlohmann::json out_json;
    // out_json["msics_subset"] = msics_subset;
    // out_json["nodes_explored"] = solver.nodes_explored;

    // std::ofstream out_file("dataAnalysis/msics_result.json");
    // out_file << out_json.dump(4);

    return 0;
}

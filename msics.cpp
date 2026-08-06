// Maximum strongly isomorphic common subdrawing (msics)
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

// Data structure storing full topological drawing invariants
struct DrawingData {
    // rot[u] = clockwise cyclic sequence of edge labels incident to vertex u
    std::vector<std::vector<std::size_t>> rot;
    // cross_seq[e] = sequence of edge labels crossed along edge e from min(u,v) to max(u,v)
    std::vector<std::vector<std::size_t>> cross_seq;
};

// Helper: Extract rotation systems and crossing sequences cleanly from HDS
template <int kplane>
DrawingData extract_drawing_data(const Drawing<kplane>& d) {
    DrawingData data;
    std::size_t num_v = d.vertices.size();
    std::size_t num_e = d.edges.size();
    data.rot.resize(num_v);
    data.cross_seq.resize(num_e);

    // 1. Rotation system around each vertex u (clockwise order using HDS pointers)
    for (std::size_t u = 0; u < num_v; ++u) {
        const HdsHalfedge* start_he = d.vertices[u].halfedge;
        if (start_he != nullptr) {
            const HdsHalfedge* curr = start_he;
            do {
                data.rot[u].push_back(curr->edge->label);
                curr = curr->twin->prev; // Step clockwise around vertex u
            } while (curr != start_he);
        }
    }

    // 2. Crossing sequences along each edge e from min(u,v) -> max(u,v)
    for (auto it = d.edges.begin(); it != d.edges.end(); ++it) {
        std::size_t e_label = it->label;
        std::size_t u_min = std::min(it->u, it->v);
        std::size_t v_max = std::max(it->u, it->v);

        // Find halfedge pointing TO u_min along edge e_label
        const HdsHalfedge* e_he = nullptr;
        const HdsHalfedge* start_he = d.vertices[u_min].halfedge;
        if (start_he != nullptr) {
            const HdsHalfedge* curr = start_he;
            do {
                if (curr->edge->label == e_label) {
                    e_he = curr;
                    break;
                }
                curr = curr->twin->prev;
            } while (curr != start_he);
        }

        if (e_he != nullptr) {
            // e_he points to u_min; e_he->twin leaves u_min towards v_max
            const HdsHalfedge* f = e_he->twin;
            while (f->vertex->label != v_max) {
                // f->vertex is a crossing vertex; f->next belongs to the crossed edge
                std::size_t crossed_label = f->next->edge->label;
                data.cross_seq[e_label].push_back(crossed_label);
                f = f->next->twin->next; // Advance along edge e past crossing
            }
        }
    }

    return data;
}

// check if two edge sequences are cyclically equivalent
bool are_cyclically_equivalent(const std::vector<std::size_t>& A, const std::vector<std::size_t>& B) {
    if (A.size() != B.size()) return false;
    if (A.empty()) return true;
    std::size_t n = A.size();
    for (std::size_t shift = 0; shift < n; ++shift) {
        bool match = true;
        for (std::size_t i = 0; i < n; ++i) {
            if (A[i] != B[(i + shift) % n]) {
                match = false;
                break;
            }
        }
        if (match) return true;
    }
    return false;
}

// check strong isomorphism for subset S across all 9 classes using INDUCED rotations and crossings
bool check_global_strong_isomorphism(
    const std::vector<DrawingData>& all_data,
    const std::vector<bool>& in_S
) {
    std::size_t num_v = all_data[0].rot.size();
    std::size_t num_e = all_data[0].cross_seq.size();

    // since reflexive, we just compare class 0 to the rest

    // verify INDUCED rotation systems at each vertex u
    for (std::size_t u = 0; u < num_v; ++u) {
        std::vector<std::size_t> r0;
        for (std::size_t e : all_data[0].rot[u])
            if (e < in_S.size() && in_S[e]) r0.push_back(e);

        for (std::size_t k = 1; k < all_data.size(); ++k) {
            std::vector<std::size_t> rk;
            for (std::size_t e : all_data[k].rot[u])
                if (e < in_S.size() && in_S[e]) rk.push_back(e);

            if (!are_cyclically_equivalent(r0, rk)) return false;
        }
    }

    // verify INDUCED crossing sequences along each edge e
    for (std::size_t e = 0; e < num_e; ++e) {
        if (e >= in_S.size() || !in_S[e]) continue;

        std::vector<std::size_t> c0;
        for (std::size_t cr : all_data[0].cross_seq[e])
            if (cr < in_S.size() && in_S[cr]) c0.push_back(cr);

        for (std::size_t k = 1; k < all_data.size(); ++k) {
            std::vector<std::size_t> ck;
            for (std::size_t cr : all_data[k].cross_seq[e])
                if (cr < in_S.size() && in_S[cr]) ck.push_back(cr);
            if (c0 != ck) return false;
        }
    }

    return true;
}

// Compute prefix divergence depth via LIFO edge popping
int compute_prefix_divergence_depth(const std::vector<nlohmann::json>& recipes) {
    int min_depth = TOTAL_EDGES;

    for (std::size_t i = 0; i < recipes.size(); ++i) {
        for (std::size_t j = i + 1; j < recipes.size(); ++j) {
            Drawing<KPLANE_LIMIT> d_i(recipes[i]);
            Drawing<KPLANE_LIMIT> d_j(recipes[j]);

            int depth = TOTAL_EDGES;
            for (std::size_t k = TOTAL_EDGES; k > 0; --k) {
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

// MSICS Backtracking Solver (Branch-and-Bound)
struct MSICSSolver {
    const std::vector<DrawingData>& all_data;
    int prefix_depth;
    std::vector<int> candidates;
    std::vector<bool> current_S;
    std::vector<bool> best_S;
    int max_size = 0;
    uint64_t nodes_explored = 0;

    MSICSSolver(const std::vector<DrawingData>& data, int pd)
        : all_data(data), prefix_depth(pd), current_S(TOTAL_EDGES, false) {
        
        for (int e = 0; e < prefix_depth; ++e) current_S[e] = true;

        // filter candidates: late-stage edges that individually work with the prefix core
        for (int e = prefix_depth; e < TOTAL_EDGES; ++e) {
            current_S[e] = true;
            if (check_global_strong_isomorphism(all_data, current_S)) candidates.push_back(e);
            current_S[e] = false;
        }

        max_size = prefix_depth;
        best_S = current_S;
    }

    void solve_dfs(int idx, int current_size) {
        nodes_explored++;

        // Pruning
        if (current_size + (static_cast<int>(candidates.size()) - idx) <= max_size) return;

        // Base case
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
        if (check_global_strong_isomorphism(all_data, current_S)) solve_dfs(idx + 1, current_size + 1);

        // 2. exclude edge e
        current_S[e] = false;
        solve_dfs(idx + 1, current_size);
    }
};

nlohmann::json create_valid_subrecipe(
    const nlohmann::json& original_recipe,
    const DrawingData& full_data,
    const std::vector<bool>& in_S
) {
    nlohmann::json sub_recipe = nlohmann::json::array();
    std::vector<bool> inserted_in_sub(in_S.size(), false);

    for (const auto& step : original_recipe["drawing_recipe"]) {
        std::size_t e = step["edge_label"].get<std::size_t>();
        
        // skip edges not in the MSICS subset
        if (e >= in_S.size() || !in_S[e]) continue;

        nlohmann::json sub_step;
        sub_step["edge_label"] = e;
        std::size_t u = step["u"].get<std::size_t>();
        std::size_t v = step["v"].get<std::size_t>();
        sub_step["u"] = u;
        sub_step["v"] = v;

        // Repair start_after_edge relative to vertex u
        const auto& rot_u = full_data.rot[u];
        auto it = std::find(rot_u.begin(), rot_u.end(), e);
        
        if (it != rot_u.end()) {
            std::size_t n = rot_u.size();
            std::size_t idx = std::distance(rot_u.begin(), it);
            
            // search counter-clockwise for the nearest already-inserted edge in S
            for (std::size_t step_back = 1; step_back < n; ++step_back) {
                std::size_t prev_idx = (idx + n - step_back) % n;
                std::size_t candidate_edge = rot_u[prev_idx];
                
                if (candidate_edge < in_S.size() && in_S[candidate_edge] && inserted_in_sub[candidate_edge]) {
                    sub_step["start_after_edge"] = candidate_edge;
                    break;
                }
            }
            // if no anchor edge exists yet at vertex u, no "start_after_edge" is set
        }

        // Filter crossed_edges to keep only those present in S
        if (step.contains("crossed_edges")) {
            auto filtered_crossed = nlohmann::json::array();
            for (const auto& cr : step["crossed_edges"]) {
                std::size_t cr_e = cr.get<std::size_t>();
                if (cr_e < in_S.size() && in_S[cr_e]) filtered_crossed.push_back(cr_e);
            }
            sub_step["crossed_edges"] = filtered_crossed;
        }

        sub_recipe.push_back(sub_step);
        inserted_in_sub[e] = true; // Mark edge as active in subdrawing
    }

    nlohmann::json result_json = original_recipe;
    result_json["drawing_recipe"] = sub_recipe;
    return result_json;
}

int main() {
    std::string base_path = "../quasiDrawings/K10_all_quasi/";
    std::vector<nlohmann::json> recipes(NUM_CLASSES);

    for (int i = 0; i < NUM_CLASSES; ++i) {
        std::string filename = base_path + std::to_string(i) + ".json";
        std::ifstream in(filename);
        if (!in.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return 1;
        }
        in >> recipes[i];
    }

    std::cout << "=====================================================\n";
    std::cout << " MAXIMUM STRONGLY ISOMORPHIC COMMON SUBDRAWING (MSICS)\n";
    std::cout << "              (EXACT combinatorial SEARCH)           \n";
    std::cout << "=====================================================\n\n";

    int prefix_depth = compute_prefix_divergence_depth(recipes);
    std::cout << "Prefix Core Depth: " << prefix_depth << " edges\n";

    // Build full HDS drawings once and extract drawing data
    std::vector<DrawingData> all_data;
    all_data.reserve(NUM_CLASSES);
    for (int i = 0; i < NUM_CLASSES; ++i) {
        Drawing<KPLANE_LIMIT> full_d(recipes[i]);
        all_data.push_back(extract_drawing_data(full_d));
    }

    MSICSSolver solver(all_data, prefix_depth);
    // std::cout << "Viable Late-Stage Candidates: " << solver.candidates.size() << " / " << (TOTAL_EDGES - prefix_depth) << " edges\n";

    solver.solve_dfs(0, prefix_depth);

    std::vector<int> msics_subset, prefix_core, late_core;

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
    bool verified = check_global_strong_isomorphism(all_data, solver.best_S);

    std::cout << "\n-----------------------------------------------------\n";
    std::cout << "[Verification] Subset S strongly isomorphic across all 9 classes: "
              << (verified ? "TRUE" : "FALSE") << "\n";
    std::cout << "-----------------------------------------------------\n\n";

    // 1. Generate repaired recipe for MSICS
    nlohmann::json valid_h24_json = create_valid_subrecipe(recipes[0], all_data[0], solver.best_S);


    // 2. Construct the Drawing object safely
    Drawing<KPLANE_LIMIT> h24_drawing(valid_h24_json);

    // // 3. Export to GraphML
    // std::ofstream graphml_out("dataAnalysis/msics_subdrawing.graphml");
    // if (graphml_out.is_open()) {
    //     h24_drawing.graphml_output(graphml_out);
    //     std::cout << "Successfully exported valid GraphML subdrawing!\n";
    // }
    
    std::ofstream of_json("dataAnalysis/msics_subdrawing.json");
    nlohmann::ordered_json output_json = h24_drawing.serialize_to_json();
    of_json << output_json.dump(4);
    of_json.close();

    // // Export JSON summary
    // nlohmann::json out_json;
    // out_json["msics_subset"] = msics_subset;
    // out_json["nodes_explored"] = solver.nodes_explored;

    // std::ofstream out_file("dataAnalysis/msics_result.json");
    // if (out_file.is_open()) {
    //     out_file << out_json.dump(4);
    //     std::cout << "Exported MSICS result to 'dataAnalysis/msics_result.json'\n";
    // }

    return 0;
}

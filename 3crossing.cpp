#include "hds_kplanar.h"
#include "iso.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12;
const std::size_t klim = 3;

const Edges outer_v2_edges = {
    {6,7},{6,8},{6,9},{6,10},{6,11},
    {7,8},{7,9},{7,10},{7,11},
    {8,9},{8,10},
    {9,10},{9,11},
    {10,11}
};

const std::vector<std::vector<std::size_t>> v1_targets_template = {
    {6,7,11}, // Slot 0 connections
    {6,7,8,9,11},  // Slot 1 connections
    {7,8,9},  // Slot 2 connections
    {8,9,10}, // Slot 3 connections
    {9,10,11},// Slot 4 connections
    {6,10,11} // Slot 5 connections
};

std::vector<Edges> generate_all_edge_sets() {
    std::vector<Edges> all_edge_sets;
    std::vector<std::size_t> p = {0, 1, 2, 3, 4, 5};
    do {
        Edges current_edges = outer_v2_edges;
        for (std::size_t slot = 0; slot < 6; ++slot) {
            std::size_t u = p[slot];
            for (std::size_t target : v1_targets_template[slot]) {
                current_edges.push_back({u, target});
            }
        }
        std::sort(current_edges.begin(), current_edges.end());
        all_edge_sets.push_back(current_edges);
    } while (std::next_permutation(p.begin()+1, p.end()));

    return all_edge_sets;
}

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

    std::vector<Edges> all_edge_sets = generate_all_edge_sets();
    std::size_t total_edge_sets = all_edge_sets.size();
    std::cout << "Generated " << total_edge_sets << " edge set permutations." << std::endl;

    for (std::size_t drawing_id = 0; drawing_id < 39; ++drawing_id) {
        std::ifstream input_file("../quasiDrawings/K6_3planar/" + std::to_string(drawing_id) + ".json");
        nlohmann::json import_data; input_file >> import_data; input_file.close();
        Drawing<klim> base_d(import_data, n);
        std::cout << "Loaded drawing " << drawing_id << std::endl;

        std::size_t perm_idx = 0;
        for (const auto& current_edges : all_edge_sets) {
            std::cout << "\rDrawing " << drawing_id << " | Permutation [" << perm_idx << " / " << total_edge_sets << "]" << std::flush;
            Drawing<klim> d = base_d;

            auto start_edge = current_edges.begin();
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

                if (++e == current_edges.end()) {
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
            perm_idx++;
        }
        std::cout << "\n Drawing " << drawing_id++ << " done\n\n";
    }
    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/K12_3planar/3cross/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/K12_3planar/3cross/" + std::to_string(idx) + ".graphml";
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

#include "hds_quasiplanar.h"
#include "iso.h"
#include <cwchar>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12; // 58 edges for optimal quasi
const std::size_t klim = 17; // leq 2n-7

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

        while (!p.empty()) {
            Drawing<klim> d_test(d);
            d_test.add_edge(p, v);
            if (d_test.verify_quasiplanarity()) {
                std::cout << "\n  [!] Edge (" << u << ", " << v << ") can be legally added!";
                return true;
            }
            if (!d_search.next_path(p, v)) {
                break;
            }
        }
    }
    return false;
}

// function to generate the 10choose5 combinations of edges
std::vector<Edges> generate_edge_sets() {
    std::vector<Edges> all_edge_sets;

    // 10 choose 6
    std::vector<int> mask_a(10, 0);
    std::fill(mask_a.end() - 6, mask_a.end(), 1);

    do {
        std::vector<std::size_t> Ca, R; // Ca = 6 neighbors of a, R = 4 remaining vertices
        for (std::size_t i = 0; i < 10; ++i) {
            if (mask_a[i]) Ca.push_back(i);
            else R.push_back(i);
        }

        // 6 choose 2
        std::vector<int> mask_b(6, 0);
        std::fill(mask_b.end() - 2, mask_b.end(), 1);

        do {
            Edges current_edges;

            // connect a (10) to its 6 chosen vertices in K10
            for (std::size_t v : Ca) current_edges.push_back({v, 10});

            // connect b (11) to all 4 remaining vertices (R)
            for (std::size_t v : R) current_edges.push_back({v, 11});

            // connect b (11) to 2 overlapping vertices
            for (std::size_t i = 0; i < 6; ++i)
                if (mask_b[i]) current_edges.push_back({Ca[i], 11});

            // current_edges.push_back({10, 11});
            all_edge_sets.push_back(current_edges);
        } while (std::next_permutation(mask_b.begin(), mask_b.end()));
    } while (std::next_permutation(mask_a.begin(), mask_a.end()));

    return all_edge_sets;
}

int main() {
    std::vector< Drawing<klim> > solutions;
    std::cout << "\n\n ===================================================== \n\n";
    std::cout << "k = " << klim << ", n = " << n  << std::endl;

    std::vector<Edges> all_edge_sets = generate_edge_sets();
    std::size_t total_edge_sets = all_edge_sets.size();

    for (int i = 0; i < 9; i++) {
        // std::cout << "K10 drawing " << i << std::endl;
        std::string name = "../quasiDrawings/K10_all_quasi/" + std::to_string(i);
        std::string jsonFile = name + ".json";
        std::ifstream input_file(jsonFile);
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        // loading drawing
        Drawing<klim> base_drawing(import_data, n);
        std::size_t verified_sets = 0;

        for (const auto& current_edges : all_edge_sets) {
            verified_sets++;

            // 1. Dynamic status line updated in-place using \r and std::flush
            std::cout << "\r[Drawing " << (i+1) << "/9] Solutions found: " << solutions.size()
                << " | Edge sets verified: " << verified_sets << "/" << total_edge_sets
                << "   " << std::flush;

            auto start_edge = current_edges.begin();
            // 3. Utilize the implemented copy constructor to reset state[cite: 1]
            Drawing<klim> d(base_drawing);
            for (auto e = start_edge;;) {
                std::size_t u = (*e)[0];
                std::size_t v = (*e)[1];
                HdsPath p = d.first_path(u, v);
                if (p.empty()) {
BACKUP:
                    // no way to add uv -> do previous edges differently
                    do {
                        if (e == start_edge) {
                            goto NEXT_EDGE_SET;
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
                    if (!d.verify_quasiplanarity()) std::cout << "not quasi?\n";

                    bool newSol = true;
                    for (auto it = solutions.begin(); it != solutions.end(); it++) {
                        if(are_isomorphic((*it),d)) {
                            newSol = false;
                            break;
                        }
                    }
                    if (newSol) {
                        if (is_drawing_extendable(d, n)) return 0;
                        solutions.push_back(d); 
                    }
                    goto BACKUP;
                }
            }
NEXT_EDGE_SET:;
        }
        std::cout << std::endl;
    }

    if (solutions.size() == 0) {
        std::cout << "NO SOLUTIONS!"  << std::endl;
        return 0;
    }
    std::cout << "GRAPH IS MAXIMAL QUASIPLANAR, TOTAL NUMBER OF SOLUTIONS: " << solutions.size() << std::endl;
    // if (solutions.size() > 50) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/maxQuasi/13_with_K3/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/maxQuasi/13_with_K3/" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }
    return 0;
}

#include "hds_quasiplanar.h"
#include "iso.h"
#include <cwchar>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 14; // 58 edges for optimal quasi
                          // const Edges edges = // 57 + 1 = 58
                          // {
                          //     {0,10},{1,10},{2,10},{3,10},{4,10},{5,10},
                          //     {4,11},{5,11},{6,11},{7,11},{8,11},{9,11},   //bistar_sep
                          // };

const std::size_t klim = 21; // leq 2n-7=21

// Helper function to generate the 10choose5 combinations of edges
std::vector<Edges> generate_edge_sets() {
    std::vector<Edges> all_edge_sets;

    // Mask for 10 choose 5 (5 zeros, 5 ones ensures sorted order)
    std::vector<int> mask_a(10, 0);
    std::fill(mask_a.end() - 5, mask_a.end(), 1);

    do {
        Edges current_edges;
        for (std::size_t i = 0; i < 10; ++i) {
            if (mask_a[i]) {
                current_edges.push_back({i, 10});
                current_edges.push_back({i, 12});
            } else {
                current_edges.push_back({i, 11});
                current_edges.push_back({i, 13});
            }
        }

        // edges between 10,11,12,13
        current_edges.push_back({10, 11}); // a-b
        current_edges.push_back({10, 12}); // a-c
        current_edges.push_back({10, 13}); // a-d
        current_edges.push_back({11, 12}); // b-c
        current_edges.push_back({11, 13}); // b-d
        // current_edges.push_back({12, 13}); // c-d

        all_edge_sets.push_back(current_edges);
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
                    if (newSol) solutions.push_back(d); 
                    goto BACKUP;
                }
            }
NEXT_EDGE_SET:;
        }
        std::cout << std::endl;
    }

    if (solutions.size() == 0) {
        std::cout << "No solutions!"  << std::endl;
        return 0;
    }
    std::cout << "TOTAL NUMBER OF SOLUTIONS: " << solutions.size() << std::endl;
    if (solutions.size() > 50) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/maxQuasi/14_v2/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/maxQuasi/14_v2/" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }
    return 0;
}

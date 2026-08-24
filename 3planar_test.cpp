#include "hds_kplanar.h"
// #include "hds_quasiplanar.h"
#include "iso.h"
//#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>

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

const Edges edges = {
    // {0,1},{0,2},{0,3},
    // {1,2},{1,3},
    // {2,3},

    {0,4},
    {4,5},{4,6},{4,7},
    {5,6},{5,7},
    {6,7},

    {0,8},
    {8,9},{8,10},{8,11},
    {9,10},{9,11},
    {10,11},

    {4,8},
    {1,5},{1,9},{5,9},
    {2,6},{2,10},{6,10},
    {3,7},{3,11},{7,11},

    {0,5},{0,11},
    {1,6},{1,8},
    {2,7},{2,9},
    {3,4},{3,10},
};

const std::size_t n = 12;
const std::size_t klim = 3;

int main() {
    // const Edges edges = generateCompleteGraph(n);
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(10000,0); // assume no more than 10000 unique drawings up to iso

    // Drawing<klim> d(n);
    // d.add_first_edge(edges[0][0], edges[0][1]);
    std::size_t num_fixed_edges = 0;
    // for (std::size_t i = 2; i < num_fixed_edges; ++i) {
    //     HdsPath p = d.first_path(0, i);
    //     if (p.empty()) {
    //         throw std::runtime_error("Failed to build the initial star!");
    //     }
    //     d.add_edge(p, i);
    // }
    // std::cout << "Star built" << std::endl;

    for (int i = 0; i < 2; i++) {
        std::cout << "Drawing " << std::to_string(i) << std::endl;

        std::ifstream input_file("../quasiDrawings/K4_3planar/" + std::to_string(i) + ".json");
        nlohmann::json import_data; input_file >> import_data; input_file.close();
        // loading drawing
        Drawing<klim> d(import_data, n);

        auto start_edge = edges.begin() + num_fixed_edges;

        std::vector<uint64_t> fail_counts(edges.size(), 0);
        std::size_t max_depth = 0;
        uint64_t step_counter = 0;
        uint64_t total_fails = 0;
        auto start_time = std::chrono::steady_clock::now();

        for (auto e = start_edge;;) {
            std::size_t e_idx = std::distance(edges.begin(), e);
            std::size_t u = (*e)[0];
            std::size_t v = (*e)[1];


            if (e_idx > max_depth) max_depth = e_idx;
            // Periodic progress report (every 1K steps)
            if (++step_counter % 10000 == 0) {
                auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
                std::cout << "\033[2J\033[1;1H";
                std::cout << "================== SEARCH PROGRESS ==================\n";
                std::cout 
                    << " | Drawing: " << i 
                    << " | Steps: " << (step_counter / 1000) << "K"
                    << " | Time: " << static_cast<uint64_t>(elapsed) << " s"
                    << " | Speed: " << static_cast<uint64_t>(step_counter / elapsed) << " st/s"
                    << " | Current Depth: " << e_idx << "/" << edges.size() - 1
                    << " | Max Depth: " << max_depth
                    << " | Solutions: " << solutions.size()
                    << " | Total fails: " << total_fails
                    << " | (n,k): (" << n << "," << klim << ")\n";
                std::cout << "-----------------------------------------------------\n";
                std::cout << "Edge Failures Breakdown:\n";

                // Print failure grid (4 edges per line)
                int cols = 0;
                for (std::size_t i = num_fixed_edges; i < edges.size(); ++i) {
                    std::string edge_label = "(" + std::to_string(edges[i][0]) + "->" + std::to_string(edges[i][1]) + ")";

                    std::cout << "E" << std::right << std::setw(2) << i << "(" 
                        << std::left << std::setw(8) << edge_label << ": " 
                        << std::right << std::setw(8) << fail_counts[i] << " | ";
                    if (++cols % 4 == 0) std::cout << "\n";
                }
                std::cout << "\n" << std::flush;
            }

            HdsPath p = d.first_path(u, v);
            if (p.empty()) {
BACKUP:
                fail_counts[std::distance(edges.begin(), e)]++; // Count edge failure
                total_fails++;
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
                if (newSol) solutions.push_back(d);
                goto BACKUP;
            }
        }
END:
    }
    // std::cout << "\n\n=== Edge Bottleneck Analysis ===\n";
    // for (std::size_t i = num_fixed_edges; i < edges.size(); ++i) {
    //     std::cout << "Edge " << std::setw(2) << i << " (" << edges[i][0] << "->" << edges[i][1] << "): " 
    //         << fail_counts[i] << " failures/backtracks\n";
    // }

    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/3planar_test/quad_triangle/" + std::to_string(idx) + "_quasi.json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/3planar_test/quad_triangle/" + std::to_string(idx) + "_quasi.graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }

    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;
}

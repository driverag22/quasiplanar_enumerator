#include "hds_kplanar.h"
// #include "hds_3_planar.h"
// #include "hds_quasiplanar.h"
#include "iso.h"
//#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

// Edges generateCompleteGraph(std::size_t n) {
//     Edges edges;
//     edges.reserve(n * (n - 1) / 2);
//     for (std::size_t i = 0; i < n; ++i) {
//         for (std::size_t j = i + 1; j < n; ++j) {
//             edges.push_back({i, j});
//         }
//     }
//     return edges;
// }

const Edges edges = {
    {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},{0,8},{0,9},{0,10},{0,11},
    {1,2},{1,3},{1,4},{1,5},{1,6},{1,7},{1,8},{1,9},{1,10},{1,11},

    {2,3},{3,4},{4,5},{5,6},
    {2,7},{7,8},{8,9},{9,10},{10,11},

    {2,4},{3,5},{7,9},{8,10},

    {3,8},{3,6},
    {4,6},{4,11},
    {5,10},
    {6,9},
    {6,11},

    // {2,3},{2,4},{2,7},
    // {3,4},{3,5},{3,8},
    // {4,5},{4,11},
    // {5,6},{5,10},
    // {6,9},{6,11},
    // {7,8},{7,9},
    // {8,9},{8,10},
    // {9,10},
    // {10,11},
};

const std::size_t n = 12;
const std::size_t klim = 3;
const std::string split = "1c";

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
    // const Edges edges = generateCompleteGraph(n);
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(10000,0); // assume no more than 10000 unique drawings up to iso

    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);
    std::size_t num_fixed_edges = 1;
    // for (std::size_t i = 2; i <= num_fixed_edges; ++i) {
    //     HdsPath p = d.first_path(0, i);
    //     if (p.empty()) {
    //         throw std::runtime_error("Failed to build the initial star!");
    //     }
    //     d.add_edge(p, i);
    // }
    // std::cout << "Star built" << std::endl;

    auto start_edge = edges.begin() + num_fixed_edges;
    // auto start_edge = edges.begin() + 1;

    int counter = 0;
    std::vector<uint64_t> fail_counts(edges.size(), 0);
    std::size_t max_depth = 0;
    uint64_t step_counter = 0;

    auto start_time = std::chrono::steady_clock::now();
    auto last_log_time = start_time;
const auto log_interval = std::chrono::seconds(60);

    for (auto e = start_edge;;) {
        std::size_t e_idx = std::distance(edges.begin(), e);
        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];
        auto now = std::chrono::steady_clock::now();
        step_counter++;

        if (e_idx > max_depth) max_depth = e_idx;

        // Periodic progress report (every 1M steps)
        if (now - last_log_time >= log_interval) {
            last_log_time = now;
            auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
            uint64_t speed = static_cast<uint64_t>(step_counter / (elapsed > 0 ? elapsed : 1));
            std::cout << "\n============================== SPLIT " << split << ", SEARCH PROGRESS (t = " << elapsed << "s) ==============================\n";
            std::cout << "Steps: " << (step_counter / 1000000) << "M | Speed: " << speed << " st/s"
                << " | Time: " << static_cast<uint64_t>(elapsed) << " s"
                << " | Current Edge: " << e_idx << "/" << edges.size() - 1
                << " | Max Depth: " << max_depth
                << " | Sols: " << solutions.size()
                << " | (n,k,|E|): (" << n << "," << klim << "," << edges.size() << ")\n";
            std::cout << "-----------------------------------------------------------------------------------------\n";
            std::cout << "Edge Failures Breakdown:\n";

            // Print failure grid (4 edges per line)
            int cols = 0;
            for (std::size_t i = num_fixed_edges; i < edges.size(); ++i) {
                std::cout << "E" << std::setw(2) << i << "(" << edges[i][0] << "->" << edges[i][1] << "): " 
                    << std::setw(9) << fail_counts[i] << " | ";
                if (++cols % 4 == 0) std::cout << "\n";
            }
            if (cols % 4 != 0) std::cout << "\n";
            std::cout << "=========================================================================================\n" << std::endl;
        }

        HdsPath p = d.first_path(u, v);
        if (p.empty()) {
BACKUP:
            fail_counts[std::distance(edges.begin(), e)]++; // Count edge failure
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
            if (newSol) {
                if (is_drawing_extendable(d,n)) {
                    std::cout << "EXTENDABLE\n\n";
                }
                solutions.push_back(d);
                std::cout << ++counter << std::endl;
            }
            goto BACKUP;
        }
    }
END:
    std::cout << "\n\n=== Edge Bottleneck Analysis ===\n";
    for (std::size_t i = num_fixed_edges; i < edges.size(); ++i) {
        std::cout << "Edge " << std::setw(2) << i << " (" << edges[i][0] << "->" << edges[i][1] << "): " 
            << fail_counts[i] << " failures/backtracks\n";
    }

    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/3planar_test/split1/" + std::to_string(idx) + "_quasi.json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();
    
        std::string filename2 = "../quasiDrawings/3planar_test/split1/" + std::to_string(idx) + "_quasi.graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }

    std::cout << "Found " << counter << " drawings in total." << std::endl;
    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;
}

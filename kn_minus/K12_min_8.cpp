#include "hds_quasiplanar.h"
#include <fstream>
#include <chrono>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12; // note hard-coded limit of 64 edges for quasiplanar...
const std::string split = "path";
const Edges edges =
{
    {0,7},{1,7},{2,7},{3,7},{4,7},{5,7},{6,7},{7,8},{7,9},{7,10},{7,11},
    {0,8},{1,8},{4,8},{5,8},{6,8},{8,9},{8,10},{8,11},
    {0,9},{1,9},{2,9},{5,9},{6,9},{9,10},{9,11},
    {0,10},{1,10},{2,10},{3,10},{6,10},{10,11},
    {0,11},{1,11},{2,11},{3,11},{4,11},
};

const std::size_t klim = 12;

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n  << ", split = " << split << std::endl;
    std::cout << "Total edges (exp. 58): 21 + " << edges.size() << " = " << (21+edges.size()) << std::endl;

    for (int i = 0; i < 1607; i++) {
        std::cout << i << std::endl;
        std::string name = "../quasiDrawings/K7/klim9_jsons/" + std::to_string(i);
        std::string jsonFile = name + ".json";
        std::ifstream input_file(jsonFile);
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        // loading drawing
        Drawing<klim> d(import_data, n);
        auto start_edge = edges.begin();

        std::vector<uint64_t> fail_counts(edges.size(), 0);
        std::size_t max_depth = 0;
        uint64_t step_counter = 0;
        uint64_t total_fails = 0;
        auto start_time = std::chrono::steady_clock::now();

        for (auto e = start_edge;;) {
            std::size_t e_idx = std::distance(start_edge, e);
            std::size_t u = (*e)[0];
            std::size_t v = (*e)[1];

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
                    << " | Total fails: " << total_fails
                    << " | (n,k): (" << n << "," << klim << ")\n";
                std::cout << "-----------------------------------------------------\n";
                std::cout << "Edge Failures Breakdown:\n";

                // Print failure grid (4 edges per line)
                int cols = 0;
                for (std::size_t j = 0; j < edges.size(); ++j) {
                    std::string edge_label = "(" + std::to_string(edges[j][0]) + "->" + std::to_string(edges[j][1]) + ")";

                    std::cout << "E" << std::right << std::setw(2) << j << "(" 
                        << std::left << std::setw(8) << edge_label << ": " 
                        << std::right << std::setw(8) << fail_counts[j] << " | ";
                    if (++cols % 4 == 0) std::cout << "\n";
                }
                std::cout << "\n" << std::flush;
            }

            HdsPath p = d.first_path(u, v);
            if (p.empty()) {
                fail_counts[std::distance(edges.begin(), e)]++; // Count edge failure
                total_fails++;
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

            if (++e == edges.end()) {
                std::cout << "Found sol for: " << split << std::endl;
                if (!d.verify_quasiplanarity()) std::cout << "not quasi?\n";

                // output json
                std::string name = "../quasiDrawings/K12_minus_8/" + split + "/klim" + std::to_string(klim) + "_" + std::to_string(i);
                std::string jsonOut = name + ".json";
                std::ofstream of_json(jsonOut);
                nlohmann::ordered_json output_json = d.serialize_to_json();
                of_json << output_json.dump(4);
                of_json.close();

                // output graphml
                std::string graphmlOut = name + ".graphml";
                std::ofstream of_graphml(graphmlOut);
                d.graphml_output(of_graphml);
                of_graphml.close();
                return 0;
            }
        }
NEXT_JSON:;
    }
    std::cout << "Found no solutions with k = " << klim << " for split = " << split << std::endl;
    return 0;
}

#include "hds_quasiplanar.h"
#include "iso.h"
//#include <sstream>
#include <fstream>

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

const std::size_t n = 8; // note hard-coded limit of 64 edges for quasiplanar...
const std::size_t klim = 12;

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n << std::endl;
    const Edges edges = generateCompleteGraph(n);
    std::size_t minimal_cr = 0x3f3f3f3f;

    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);
    std::size_t num_fixed_edges = n;
    for (std::size_t i = 2; i < num_fixed_edges; ++i) {
        HdsPath p = d.first_path(0, i);
        if (p.empty()) {
            throw std::runtime_error("Failed to build the initial star!");
        }
        d.add_edge(p, i);
    }
    std::cout << "Star built" << std::endl;
    std::size_t cnt = 0;

    auto start_edge = edges.begin() + (num_fixed_edges-1);
    int counter = 0;
    for (auto e = start_edge;;) {
        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];
        HdsPath p = d.first_path(u, v);
        if (p.empty()) {
BACKUP:
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
            if (d.verify_quasiplanarity()) minimal_cr = std::min(minimal_cr, d.crossings.size());
            else {
                std::cout << "not quasiplanar" << std::endl;

                std::string jsonOut = "./drawings/K8/notQuasi" + std::to_string(cnt) + ".json";
                std::ofstream of_json(jsonOut);
                nlohmann::ordered_json output_json = d.serialize_to_json();
                of_json << output_json.dump(4);
                of_json.close();

                cnt++;
            }
            goto BACKUP;
        }
    }
END:
    std::cout << "Min crossing drawing for K" << n << " has " << minimal_cr << " crossings." << std::endl;
    // std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    // if(solutions.size() == 0) {
    //     return 0;
    // }

    // std::size_t idx = 0;
    // for (auto it = solutions.begin();it!=solutions.end();it++) {
    //     std::string filename = "drawings/K7_k9/jsons/" + std::to_string(idx) + ".json";
    //     std::ofstream of_json(filename);
    //     nlohmann::ordered_json output_json = (*it).serialize_to_json();
    //     of_json << output_json.dump(4);
    //     of_json.close();

    //     // std::string filename2 = "drawings/K7_prop_test/minCr_" + std::to_string(idx) + ".graphml";
    //     // std::ofstream of_graphml(filename2);
    //     // (*it).graphml_output(of_graphml);
    //     // of_graphml.close();
    //     idx++;
    // }

    // std::cout << "Found " << counter << " drawings in total." << std::endl;
    // std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;
    // std::cout << "Found " << idx << " min crossing drawings in total." << std::endl;
    // std::cout << "Minimal Crossing Number is "<<minimal_cr<<std::endl;

    // for (std::size_t i = 0; i < solutions.size(); i++) {
    //     std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    // }
    return 0;
}

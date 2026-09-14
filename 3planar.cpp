#include "hds_kplanar.h"
#include "iso.h"
#include <fstream>
// #include <set>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const Edges edges = {
    {2,9}, 
    {3,9},
    {5,9},
    {0,10},
    {1,10},
    {6,10},
    {1,11},
    {6,11},
    {3,12},
    {5,12},
    {10,11},
    {9,12},
};

const std::size_t n = 13;
const std::size_t klim = 3;
const std::string src = "extension_8_rel";
const std::string dir = "k8_8_9_10";

int main() {
    std::vector< Drawing<klim> > solutions;
    std::size_t cnt = 0;
    for (std::size_t d_n = 0; d_n < 5; d_n++) {
        std::ifstream input_file("../quasiDrawings/K8_3planar/extension_deg5_vertex/" + std::to_string(d_n) + ".json");
        nlohmann::json import_data; 
        input_file >> import_data; 
        input_file.close();
        Drawing<klim> d(import_data, n);

        const auto start_edge = edges.begin();
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
                if (d_n != 4) std::cout << "Different solution!!!\n";
                bool newSol = true;
                for (auto it = solutions.begin(); it != solutions.end(); it++) {
                    if(are_isomorphic((*it),d)) {
                        newSol = false;
                        break;
                    }
                }
                if (newSol) {
                    solutions.push_back(d);
                    for (const auto& e : d.edges) {
                        if (e.u == 0 && e.v == 2) {
                            if (e.ncr != 0) {
                                std::cout << "CROSSED\n";
                            }
                        } else if (e.u == 1 && e.v == 7) {
                            if (e.ncr != 3) {
                                std::cout << "(1,7) NOT CROSSED ENOUGH\n";
                            }
                        } else if (e.u == 3 && e.v == 7) {
                            if (e.ncr != 3) {
                                std::cout << "(3,7) NOT CROSSED ENOUGH\n";
                            }
                        }
                    }
                }
                goto BACKUP;
            }
        }
END:
        std::cout << cnt++ << std::endl;
    }

    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::cout << "Drawing " << idx << "\n";
        if ((*it).is_drawing_extendable()) {
            std::cout << "EXTENDABLE" << std::endl;
        }
        // std::string filename = "../quasiDrawings/K8_3planar/extension_deg5_vertex/" + dir + "/" + std::to_string(idx) + ".json";
        // std::ofstream of_json(filename);
        // nlohmann::ordered_json output_json = (*it).serialize_to_json();
        // of_json << output_json.dump(4);
        // of_json.close();

        // std::string filename2 = "../quasiDrawings/K8_3planar/extension_deg5_vertex/" + dir + "/" + std::to_string(idx) + ".graphml";
        // std::ofstream of_graphml(filename2);
        // (*it).graphml_output(of_graphml);
        // of_graphml.close();
        idx++;
    }

    return 0;
}

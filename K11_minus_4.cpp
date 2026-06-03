#include "hds_quasiplanar.h"
#include "iso.h"
#include <sstream>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 11;
const std::size_t klim = 9;

int main() {
    std::string name = "jsons/K10/K10_plus_one";
    std::string filename = name + ".json";
    std::ifstream input_file(filename);
    nlohmann::json import_data;
    input_file >> import_data;
    input_file.close();

    // loading drawing
    Drawing<klim> loaded_drawing(import_data);

    const Edges edges = {
        {1,10},{5,10},{7,10},{8,10},{9,10}
    };

    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(300,1); // assume no more than 300 unique drawings up to iso
    std::size_t counter = 1;

    auto start_edge = edges.begin();
    for (auto e = start_edge;;) {
        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];
        HdsPath p = loaded_drawing.first_path(u, v);
        if (p.empty()) {
BACKUP:
            // no way to add uv -> do previous edges differently
            do {
                if (e == start_edge) {
                    goto END;
                }
                --e;

                u = (*e)[0];
                assert(u == loaded_drawing.edges.back().u);
                v = (*e)[1];
                assert(v == loaded_drawing.edges.back().v);
                p = loaded_drawing.edges.back().built;
                loaded_drawing.remove_edge();
            } while (!loaded_drawing.next_path(p, v));
        }
        loaded_drawing.add_edge(p, v);

        if (++e == edges.end()) {
            bool newSol = true;
            std::size_t d_ind = 0;
            for (auto it = solutions.begin(); it != solutions.end(); it++) {
                if(are_isomorphic((*it),loaded_drawing)) {
                    newSol = false;
                    d_cnt[d_ind]++;
                    break;
                }
                d_ind++;
            }
            if (newSol) {
                if (!loaded_drawing.verify_quasiplanarity()) {
                } else {
                    solutions.push_back(loaded_drawing);
                    std::cout << ++counter << std::endl;
                    goto BACKUP;
                }
            }
            goto BACKUP;
        }
    }
END:
    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) {
        return 0;
    }

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "jsons/K11_minus_4/star/star_" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "jsons/K11_minus_4/star/star_" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();

        idx++;
    }

    std::cout << "Found " << counter << " drawings in total." << std::endl;
    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;
    // std::cout << "Minimal Crossing Number is "<<minimal_cr<<std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;

    return 0;
}

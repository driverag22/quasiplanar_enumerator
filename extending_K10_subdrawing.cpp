#include "hds_quasiplanar.h"
#include "iso.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 13; // note hard-coded limit of 64 edges for quasiplanar...
const Edges edges =
{
    {0,10}, {1,10}, {2,10}, {3,10}, {4,10}, {5,10}, {6,10}, {7,10}, {8,10}, {9,10},
    {0,11}, {1,11}, {2,11}, {3,11}, {4,11}, {5,11}, {6,11}, {7,11}, {8,11}, {9,11},
    {10,11},
    {0,12}, {1,12}, {2,12}, {3,12}, {4,12}, {5,12}, {6,12}, {7,12}, {8,12}, {9,12},
    {10,12},{11,12},
};
const std::size_t klim = 19;

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n << std::endl;
    std::vector< Drawing<klim> > solutions;
    std::string name = "../quasiDrawings/K10_subdrawing/subdrawing";
    std::string jsonFile = name + ".json";
    std::ifstream input_file(jsonFile);
    nlohmann::json import_data;
    input_file >> import_data;
    input_file.close();
    int counter = 0;

    // loading drawing
    Drawing<klim> d(import_data, n);
    auto start_edge = edges.begin();

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
            if (!d.verify_quasiplanarity()) std::cout << "not quasi?\n";
            bool newSol = true;
            for (auto it = solutions.begin(); it != solutions.end(); it++) {
                if(are_isomorphic((*it),d)) {
                    newSol = false;
                    break;
                }
            }
            if (newSol) {
                solutions.push_back(d);
                std::cout << ++counter << std::endl;
            }
            goto BACKUP;
        }
    }
END:
    std::cout << "Found " << solutions.size() << " solutions" << std::endl;
    // for (auto it = solutions.begin(); it != solutions.end(); it++) {
            // // output json
            // std::string name = "../quasiDrawings/K11_minus_4/" + split + "/klim" + std::to_string(klim) + "_" + std::to_string(i);
            // std::string jsonOut = name + ".json";
            // std::ofstream of_json(jsonOut);
            // nlohmann::ordered_json output_json = d.serialize_to_json();
            // of_json << output_json.dump(4);
            // of_json.close();

            // // output graphml
            // std::string graphmlOut = name + ".graphml";
            // std::ofstream of_graphml(graphmlOut);
            // d.graphml_output(of_graphml);
            // of_graphml.close();
    // }
    return 0;
}

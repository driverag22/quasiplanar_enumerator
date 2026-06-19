#include "hds_quasiplanar.h"
#include "iso.h"
#include <cwchar>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12; // 58 edges for optimal quasi
const Edges edges = // 6+6+1 (+45 = 58)
{
    {0,10},{1,10},{2,10},{3,10},{4,10},{5,10},
    {4,11},{5,11},{6,11},{7,11},{8,11},{9,11},
    {10,11},
};

const std::size_t klim = 17; // leq 2n-7=17

int main() {
    std::vector< Drawing<klim> > solutions;
    std::size_t counter = 0;
    std::cout << "\n\n ===================================================== \n star\n";
    std::cout << "k = " << klim << ", n = " << n  << std::endl;

    for (int i = 0; i < 9; i++) {
        std::cout << "K10 drawing " << i << std::endl;
        std::string name = "../quasiDrawings/K10_all_quasi/" + std::to_string(i);
        std::string jsonFile = name + ".json";
        std::ifstream input_file(jsonFile);
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

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
                    std::cout << counter++ << std::endl;
                }
                goto BACKUP;
            }
        }
NEXT_JSON:;
    }
    std::size_t idx = 0;
    if(solutions.size() == 0) {
        std::cout << "Not quasiplanar!" << std::endl; return 0;
    }
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/maxQuasi/K12_min_star/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/maxQuasi/K12_min_star/" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }
    return 0;
}

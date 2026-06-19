#include "hds_quasiplanar.h"
#include "iso.h"
#include <cwchar>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12; // 58 edges for optimal quasi
const Edges edges = // 57 + 1 = 58
{
    {0,11},{1,11},{2,11},{3,11},
    {6,10},{7,10},{8,10},{9,10},
    {10,11},

};

const std::size_t klim = 17; // leq 2n-7=17

int main() {
    std::vector< Drawing<klim> > solutions;
    std::size_t counter = 0;
    std::cout << "\n\n ===================================================== \n star\n";
    std::cout << "k = " << klim << ", n = " << n  << std::endl;

    for (int i = 0; i < 6; i++) {
        std::cout << "57 edge drawing " << i << std::endl;
        std::string name = "../quasiDrawings/maxQuasi/K12_min_star/" + std::to_string(i);
        std::string jsonFile = name + ".json";
        std::ifstream input_file(jsonFile);
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        // loading drawing
        Drawing<klim> d(import_data, n);

        std::size_t eC = 0;
        for (const auto& e : edges) {
            std::cout << eC++ << std::endl;
            std::size_t u = e[0];
            std::size_t v = e[1];
            
            HdsPath p = d.first_path(u, v);
            // check all possible paths for this single edge
            while (!p.empty()) {
                d.add_edge(p,v);
                if (d.verify_quasiplanarity()) {
                    std::cout << "  -> SUCCESS: Added edge " << u << "-" << v << ". Graph is NOT maximal." << std::endl;
                    solutions.push_back(d); 
                    std::cout << counter++ << std::endl;
                }
                d.remove_edge();

                if (!d.next_path(p,v)) {
                    break;
                }
            }
        }
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

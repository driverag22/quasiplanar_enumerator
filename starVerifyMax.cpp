#include "hds_quasiplanar.h"
#include "iso.h"
#include <cwchar>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12; // 58 edges for optimal quasi
const Edges edges = // 57 + 1 = 58
{
    // for bistar_eq
    // {6,10},{7,10},{8,10},{9,10},
    // {5,11},{6,11},{7,11},{8,11},{9,11},

    // for bistar_sep
    {6,10},{7,10},{8,10},{9,10},
    {0,11},{1,11},{2,11},{3,11},
    {10,11},
};

const std::size_t klim = 17; // leq 2n-7=17

int main() {
    std::vector< Drawing<klim> > solutions;
    std::size_t counter = 0;
    std::cout << "\n\n ===================================================== \n star\n";
    std::cout << "k = " << klim << ", n = " << n  << std::endl;

    for (int i = 0; i < 6; i++) {
        std::cout << "56 edge drawing " << i << std::endl;
        std::string name = "../quasiDrawings/maxQuasi/K12_min_bistar_sep/" + std::to_string(i);
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
     if (solutions.size() == 0) 
         std::cout << "Maximal quasiplanar"  << std::endl;
     else 
         std::cout << "Not max quasiplanar!" << std::endl; 

    return 0;
}

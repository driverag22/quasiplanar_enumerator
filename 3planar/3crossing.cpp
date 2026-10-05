#include "hds_kplanar.h"
#include "iso.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12;
const std::size_t klim = 3;
const std::string split = "3";

const Edges edges = {
    {3,11}
};

int main() {
    const std::string base_path = "../quasiDrawings/K12_3planar/3cross/split2/";
    for (int i = 0; i <= 17; i++) {
        std::string filename = base_path + std::to_string(i) + ".json";
        std::ifstream input_file(filename);
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        Drawing<klim> d(import_data);

        // Check if edge (3, 11) can be added legally
        HdsPath p = d.first_path(3, 11);

        if (!p.empty()) {
            std::cout << "[SUCCESS] Drawing " << i << ": Edge CAN be legally added.\n";
            d.add_edge(p,11);

            std::string filename = "../quasiDrawings/K12_3planar/3cross/split2/" + std::to_string(i) + "_plus_3_11.graphml";
            std::ofstream of_graphml(filename);
            d.graphml_output(of_graphml);
            of_graphml.close();
        } else {
            std::cout << "[FAILED]  Drawing " << i << ": Edge CANNOT be added.\n";
        }
    }
    return 0;
}

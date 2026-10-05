#include "hds_quasiplanar.h"
#include "iso.h"
#include <cwchar>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 14;
const std::size_t klim = 21; // leq 2n-7=21

int main() {
    for (int i = 0; i < 18; i++) {
        std::cout << "Drawing " << i << std::endl;
        std::string filename = "../quasiDrawings/maxQuasi/14_with_K4/" + std::to_string(i) + ".json";
        std::ifstream input_file(filename);
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        Drawing<klim> d(import_data, n);

        if (d.is_drawing_extensible()) {
            std::cout << " -> EXTENSIBLE (Not maximal)\n";
        } else {
            std::cout << " -> MAXIMAL\n";
        }
    }

    return 0;
}

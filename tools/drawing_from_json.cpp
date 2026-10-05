#include "hds_quasiplanar.h"
#include <fstream>

const std::size_t klim = 9;

int main() {
    std::string name = "jsons/K10/K10_plus_one";
    // reading json
    std::string filename = name + ".json";
    std::ifstream input_file(filename);
    nlohmann::json import_data;
    input_file >> import_data;
    input_file.close();

    // loading drawing
    Drawing<klim> loaded_drawing(import_data);
    std::cout << loaded_drawing << std::endl;

    // output graphml
    std::string filename2 = name + "_test.graphml";
    std::ofstream of_graphml(filename2);
    loaded_drawing.graphml_output(of_graphml);
    of_graphml.close();

    // output json
    std::ofstream of_json("K0_plus_one_test.json");
    nlohmann::ordered_json output_json = loaded_drawing.serialize_to_json();
    of_json << output_json.dump(4);
    of_json.close();

  return 0;
}

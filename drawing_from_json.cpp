#include "hds_quasiplanar.h"
#include "iso.h"
#include <sstream>
#include <fstream>

const std::size_t klim = 5;

int main() {
  // reading json
  std::string filename = "K9_original.json";
  std::ifstream input_file(filename);
  nlohmann::json import_data;
  input_file >> import_data;
  input_file.close();

  // loading drawing
  Drawing<klim> loaded_drawing(import_data);

  // output graphml
  std::ofstream of_graphml("K9_test.graphml");
  loaded_drawing.graphml_output(of_graphml);
  of_graphml.close();

  // output json
  std::ofstream of_json("K9_test.json");
  nlohmann::json output_json = loaded_drawing.serialize_to_json();
  of_json << output_json.dump(4);
  of_json.close();

  return 0;
}

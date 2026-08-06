#include "hds_quasiplanar.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 12; // note hard-coded limit of 64 edges for quasiplanar...
const std::string split = "1";
const Edges edges = // 57 edges
{
    {0,7},{1,7},{2,7},{3,7},{4,7},
    // {5,7},
    {1,8},{2,8},{3,8},{4,8},{5,8},{6,8},
    {0,9},{1,9},{2,9},{3,9},{4,9},
    {1,10},{2,10},{3,10},{4,10},{5,10},
    {2,11},{3,11},{4,11},{5,11},{6,11},
    // K5
    {7,8},{7,9},{7,10},{7,11},
    {8,9},{8,10},{8,11},
    {9,10},{9,11},
    {10,11},
};

const std::vector<std::pair<std::size_t,std::size_t>> missingEdges =
{
    {5,7},
    {6,7},
    {0,8},
    {5,9},{6,9},
    {0,10},{6,10},
    {0,11},{1,11},
};
const std::size_t klim = 17; // 2n-7 = 17

bool is_drawing_extendable(const Drawing<klim>& d, std::size_t num_vertices) {
    for (const auto& [u,v]: missingEdges) {
        Drawing<klim> d_search(d);
        HdsPath p = d_search.first_path(u,v);

        while (!p.empty()) {
            Drawing<klim> d_test(d);
            d_test.add_edge(p, v);
            if (d_test.verify_quasiplanarity()) {
                std::cout << "\n  [!] Edge (" << u << ", " << v << ") can be legally added!";
                return true;
            } else std::cout << "not quasi!\n";
            if (!d_search.next_path(p, v)) break;
        }
    }
    return false;
}

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n  << ", split = " << split << std::endl;

    for (int i = 0; i < 1607; i++) {
        std::cout << "Drawing " << std::to_string(i) << std::endl;
        std::string name = "../quasiDrawings/K7/klim9_jsons/" + std::to_string(i);
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
                std::cout << "Found sol for: " << split << std::endl;
                if (!d.verify_quasiplanarity()) std::cout << "not quasi?\n";
                if (is_drawing_extendable(d, n)) std::cout << "extendable\n";

                // output json
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
                return 0;
            }
        }
NEXT_JSON:;
    }
    std::cout << "Found no solutions with k = " << klim << " for split = " << split << std::endl;
    return 0;
}

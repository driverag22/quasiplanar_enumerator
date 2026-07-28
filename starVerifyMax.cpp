#include "hds_quasiplanar.h"
#include "iso.h"
#include <cwchar>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 14;
const std::size_t klim = 21; // leq 2n-7=21

// Helper function to compute missing edges for a drawing
std::vector<std::pair<std::size_t, std::size_t>> get_missing_edges(const Drawing<klim>& d, std::size_t num_vertices) {
    std::vector<std::vector<bool>> adj(num_vertices, std::vector<bool>(num_vertices, false));
    for (const auto& edge : d.edges) {
        adj[edge.u][edge.v] = true;
        adj[edge.v][edge.u] = true;
    }

    std::vector<std::pair<std::size_t, std::size_t>> missing;
    for (std::size_t u = 0; u < num_vertices; ++u)
        for (std::size_t v = u + 1; v < num_vertices; ++v)
            if (!adj[u][v]) missing.push_back({u, v});
    return missing;
}

bool is_drawing_extendable(const Drawing<klim>& d, std::size_t num_vertices) {
    auto missing_edges = get_missing_edges(d, num_vertices);

    for (const auto& [u, v] : missing_edges) {
        Drawing<klim> d_search(d);
        HdsPath p = d_search.first_path(u, v);

        while (!p.empty()) {
            Drawing<klim> d_test(d);
            d_test.add_edge(p, v);
            if (d_test.verify_quasiplanarity()) {
                std::cout << "\n  [!] Edge (" << u << ", " << v << ") can be legally added!";
                return true;
            }
            if (!d_search.next_path(p, v)) {
                break;
            }
        }
    }
    return false;
}

int main() {
    bool graph_is_maximal = true;

    for (int i = 0; i < 18; i++) {
        std::cout << "Drawing " << i << std::endl;
        std::string filename = "../quasiDrawings/maxQuasi/14_with_K4/" + std::to_string(i) + ".json";
        std::ifstream input_file(filename);
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        Drawing<klim> d(import_data, n);

        if (is_drawing_extendable(d, n)) {
            std::cout << " -> EXTENDABLE (Not maximal)\n";
            graph_is_maximal = false;
        } else {
            std::cout << " -> MAXIMAL\n";
        }
    }

    return 0;
}

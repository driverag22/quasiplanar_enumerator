#include "hds_kplanar.h"
#include "iso.h"
#include <fstream>
// #include <set>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const Edges edges = {
    {2,9}, 
    {3,9},
    {5,9},
    {0,10},
    {1,10},
    {6,10},
};

const std::size_t n = 11;
const std::size_t klim = 3;
const std::string src = "k8_8_9_10";
const std::string dir = "k8_8_9_10";

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

        if (!p.empty()) {
            std::cout << "Edge: " << u << " , " << v << "\n";
            return true;
        }
    }
    return false;
}

void checkUncr(std::size_t u, std::size_t v) {
    std::cout << "\n\n\n\n" << "(" << u << "," << v << ")\n\n";
    for (std::size_t d_n = 0; d_n < 5; ++d_n) {
        std::ifstream input_file("../quasiDrawings/K8_3planar/extension_deg5_vertex/" + src + "/" + std::to_string(d_n) + ".json");
        nlohmann::json import_data; input_file >> import_data; input_file.close();
        Drawing<klim> d(import_data, n);
        for (const auto& e : d.edges) {
            if (e.u == u && e.v == v) {
                if (e.ncr != 0) {
                    std::cout << "CROSSED\n";
                    return;
                }
            }
        }
    }
}
void checkExt() {
    for (std::size_t d_n = 0; d_n < 5; ++d_n) {
        std::ifstream input_file("../quasiDrawings/K8_3planar/extension_deg5_vertex/" + src + "/" + std::to_string(d_n) + ".json");
        nlohmann::json import_data; input_file >> import_data; input_file.close();
        Drawing<klim> d(import_data, n);
        if (is_drawing_extendable(d,n)) {
            std::cout << "extendable " << std::to_string(d_n) << "\n";
        }
    }
}

int main() {
    std::vector< Drawing<klim> > solutions;
    std::size_t cnt = 0;
    for (std::size_t d_n = 0; d_n < 5; d_n++) {
        std::ifstream input_file("../quasiDrawings/K8_3planar/extension_deg5_vertex/" + src + "/" + std::to_string(d_n) + ".json");
        nlohmann::json import_data; 
        input_file >> import_data; 
        input_file.close();
        Drawing<klim> d(import_data, n);

        const auto start_edge = edges.begin();
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
                bool newSol = true;
                for (auto it = solutions.begin(); it != solutions.end(); it++) {
                    if(are_isomorphic((*it),d)) {
                        newSol = false;
                        break;
                    }
                }
                if (newSol) {
                    solutions.push_back(d);
                    for (const auto& e : d.edges) {
                        if (e.u == 0 && e.v == 2) {
                            if (e.ncr != 0) {
                                std::cout << "CROSSED\n";
                            }
                        }
                    }
                }
                goto BACKUP;
            }
        }
END:
        std::cout << cnt++ << std::endl;
    }

    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::cout << "Drawing " << idx << "\n";
        if (is_drawing_extendable(*it,n)) {
            std::cout << "EXTENDABLE" << std::endl;
        }
        std::string filename = "../quasiDrawings/K8_3planar/extension_deg5_vertex/" + dir + "/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/K8_3planar/extension_deg5_vertex/" + dir + "/" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }

    return 0;
}

#include "hds_quasiplanar.h"
#include "iso.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const Edges edges = 
{
    {0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},{0,8},{0,9},
    {1,2},{1,3},{1,4},{1,5},{1,6},{1,7},{1,8},{1,9},
    {2,3},{2,4},{2,5},{2,6},{2,7},{2,8},{2,9},
    {3,4},{3,5},{3,6},{3,7},{3,8},{3,9},
    {4,5},{4,6},{4,7},{4,8},{4,9},
    {5,6},{5,7},{5,8},{5,9},
    {6,7},{6,8},{6,9},
    {7,8},{7,9},
    {8,9},
};

const std::size_t n = 10;
const std::size_t klim = 13;

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n << std::endl;
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt; // d_cnt[i]: number of further drawings found isomorphic to solutions[i]
    std::size_t minimal_cr = SIZE_MAX;


    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);
    for (std::size_t i = 1; i < n-1; ++i) {
        HdsPath p = d.first_path(edges[i][0], edges[i][1]);
        if (p.empty()) throw std::runtime_error("Failed to build the initial star!");
        d.add_edge(p, edges[i][1]);
    }
    std::cout << "Star built" << std::endl;
    auto start_edge = edges.begin() + n - 1;

    int counter = 0;
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
                v = (*e)[1];
                p = d.edges.back().built;
                d.remove_edge();
            } while (!d.next_path(p, v));
        }
        d.add_edge(p, v);

        if (++e == edges.end()) {
            bool newSol = true;
            std::size_t d_ind = 0;
            for (auto it = solutions.begin(); it != solutions.end(); it++) {
                if(are_isomorphic((*it),d)) {
                    newSol = false;
                    d_cnt[d_ind]++;
                    break;
                }
                d_ind++;
            }
            if (newSol) {
                std::cout << "Found sol with gcr: " << d.crossings.size() << "  and lcr: " << klim << std::endl;
                if (d.crossings.size() < minimal_cr) minimal_cr = d.crossings.size();
                solutions.push_back(d);
                d_cnt.push_back(0);
                std::cout << ++counter << std::endl;
            }
            goto BACKUP;
        }
    }
END:
    if(counter == 0) return 0;

    // std::size_t idx = 0;
    // for (auto it = solutions.begin();it!=solutions.end();it++) {
    //     std::string filename = "../quasiDrawings/K10_all_quasi/" + std::to_string(idx) + ".json";
    //     std::ofstream of_json(filename);
    //     nlohmann::ordered_json output_json = (*it).serialize_to_json();
    //     of_json << output_json.dump(4);
    //     of_json.close();

    //     std::string filename2 = "../quasiDrawings/K10_all_quasi/" + std::to_string(idx) + ".graphml";
    //     std::ofstream of_graphml(filename2);
    //     (*it).graphml_output(of_graphml);
    //     of_graphml.close();
    //     idx++;
    // }

    std::cout << "Found " << counter << " drawings in total." << std::endl;
    std::cout << "Minimal Crossing Number is " << minimal_cr << std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++)
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    return 0;
}

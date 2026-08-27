#include "hds_kplanar.h"
#include "iso.h"
#include <fstream>
#include <bitset>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const Edges edges = {
    {8,9},
    {9,10},{9,11},{9,12},{9,13},{9,14},
    {10,11},{10,12},{10,13},{10,14},{10,15},{10,16},{10,17},
    {11,12},{11,13},{11,14},{11,15},{11,16},{11,17},
    {12,13},{12,14},{12,15},{12,16},{12,17},
    {13,14},{13,15},{13,16},{13,17},
    {14,15},{14,16},{14,17},
    {15,16},{15,17},
    {16,17},
};

const std::size_t n = 18;
const std::size_t klim = 3;

int main() {
    std::vector< Drawing<klim> > solutions;
    std::vector<std::bitset<3>> solution_sources;
    std::vector<std::size_t> d_cnt(10000,0); // assume no more than 10000 unique drawings up to iso

    // std::vector<int> mask = {0,0,0,1,1,1,1,1};
    // int i = 0;
    // do {
    //     std::cout << "\n\n PERMUTATION: " << i++ << "\n";
        for (int i = 0; i < 5; i++) {
            std::cout << "Drawing " << std::to_string(i) << std::endl;

            std::ifstream input_file("../quasiDrawings/K8_3planar/extension_deg5_vertex/" + std::to_string(i) + ".json");
            nlohmann::json import_data; input_file >> import_data; input_file.close();
            // loading drawing
            Drawing<klim> d(import_data, n);
            // for (std::size_t j = 0; j < 8; j++) 
            //     if (mask[j]) edges.push_back({j,8});

            auto start_edge = edges.begin();

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
                    std::size_t d_ind = 0;
                    for (auto it = solutions.begin(); it != solutions.end(); it++) {
                        if(are_isomorphic((*it),d)) {
                            newSol = false;
                            d_cnt[d_ind]++;
                            solution_sources[d_ind].set(i);
                            break;
                        }
                        d_ind++;
                    }
                    if (newSol) {
                        solutions.push_back(d);
                        std::bitset<3> src;
                        src.set(i);
                        solution_sources.push_back(src);
                        d_cnt[d_ind] = 1;
                    }
                    goto BACKUP;
                }
            }
END:
            std::cout << "NEXT DRAWING\n\n\n";
        }
    // } while (std::next_permutation(mask.begin(), mask.end()));

    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) return 0;

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/K8_3planar/extension_deg5_vertex/double/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/K8_3planar/extension_deg5_vertex/double/" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }

    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;
    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] 
                  << " isomorphic drawings, sources (bits 2,1,0): " 
                  << solution_sources[i] << std::endl;
    }
    return 0;
}

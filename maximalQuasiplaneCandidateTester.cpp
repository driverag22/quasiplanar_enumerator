#include "hds_quasiplanar.h"
#include "iso.h"
#include <fstream>

const std::size_t n = 16;
const std::size_t klim = 10;

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;
const Edges edges =
{
    {0,6},{0,12},{0,13},{0,14},{0,15},
    {1,5},{1,7},{1,9},{1,13},{1,14},{1,15},
    {2,5},{2,6},{2,10},{2,12},{2,14},{2,15},
    {3,5},{3,6},{3,7},{3,13},{3,15},
    {4,10},
    {5,9},{5,11},{5,13},
    {6,9},{6,10},{6,14},
    {7,9},{7,10},{7,11},
    {8,14},
    {9,13},{9,15},
    {10,13},{10,14},
    {11,13},{11,14},{11,15},
    {12,13},{12,14},{12,15},
    {13,14},{13,15},
    {14,15}
};
// {
//     {0,6},{0,12},{0,13},{0,14},{0,15},
//     {1,5},{1,7},{1,9},{1,13},{1,14},{1,15},
//     {2,5},{2,6},{2,10},{2,12},{2,14},{2,15},
//     {3,5},{3,6},{3,7},{3,13},{3,15},
//     {4,5},{4,6},{4,7},{4,8},{4,10},
//     {5,6},{5,7},{5,8},{5,9},{5,11},{5,13},
//     {6,7},{6,8},{6,9},{6,10},{6,14},
//     {7,8},{7,9},{7,10},{7,11},
//     {8,9},{8,10},{8,11},{8,12},{8,14},
//     {9,10},{9,11},{9,12},{9,13},{9,15},
//     {10,11},{10,12},{10,13},{10,14},
//     {11,12},{11,13},{11,14},{11,15},
//     {12,13},{12,14},{12,15},
//     {13,14},{13,15},
//     {14,15}
// };

int main() {
    std::vector< Drawing<klim> > solutions;
    std::size_t counter = 0;

    // Loop over each of the 5 K_5 drawings
    for (int i = 0; i < 5; i++) {
        nlohmann::json data_i; {
            std::ifstream f("../quasiDrawings/K5/jsons/" + std::to_string(i) + ".json");
            f >> data_i;
        }

        for (int j = 0; j < 5; j++) {
            nlohmann::json data_j; {
                std::ifstream f("../quasiDrawings/K5/jsons/" + std::to_string(j) + ".json");
                f >> data_j;
            }

            for (int k = 0; k < 5; k++) {
                nlohmann::json data_k; {
                    std::ifstream f("../quasiDrawings/K5/jsons/" + std::to_string(k) + ".json");
                    f >> data_k;
                }
                
                Drawing<klim> base_d(data_i, n); // pass n to ensure vertices initialized large enough

                std::vector<HdsHalfedge*> anchors_v4;
                HdsHalfedge* start_v4 = base_d.vertices[4].halfedge;
                if (start_v4 != nullptr) {
                    HdsHalfedge* curr = start_v4;
                    do {
                        anchors_v4.push_back(curr);
                        curr = curr->twin->prev; // rotates through incident edges
                    } while (curr != start_v4);
                }

                for (HdsHalfedge* anchor4 : anchors_v4) {
                    Drawing<klim> d_with_2(base_d);

                    // get exact memory location of the anchor inside the new object
                    HdsHalfedge* corresponding_anchor4 = nullptr;
                    for (auto& he : d_with_2.halfedges) {
                        if (he.label == anchor4->label) {
                            corresponding_anchor4 = &he;
                            break;
                        }
                    }

                    // inject second K5 {4, ..., 8}
                    if (!d_with_2.inject_shifted_recipe(data_j, 4, 10, corresponding_anchor4))
                        continue;

                    // find the 4 face slot anchors pointing to vertex 8
                    std::vector<HdsHalfedge*> anchors_v8;
                    HdsHalfedge* start_v8 = d_with_2.vertices[8].halfedge;
                    if (start_v8 != nullptr) {
                        HdsHalfedge* curr = start_v8;
                        do {
                            anchors_v8.push_back(curr);
                            curr = curr->twin->prev;
                        } while (curr != start_v8);
                    }

                    for (HdsHalfedge* anchor8 : anchors_v8) {
                        Drawing<klim> d(d_with_2);

                        HdsHalfedge* corresponding_anchor8 = nullptr;
                        for (auto& he : d.halfedges) {
                            if (he.label == anchor8->label) {
                                corresponding_anchor8 = &he;
                                break;
                            }
                        }

                        // inject third K5 {8, ..., 12}
                        if (!d.inject_shifted_recipe(data_k, 8, 20, corresponding_anchor8))
                            continue;

                        std::cout << "Combined drawings: " << i << ", " << j << ", " << k << std::endl;
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
                                        goto NEXT_COMBINATION;
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
                                if (!d.verify_quasiplanarity()) std::cout << "not quasi?\n";

                                bool newSol = true;
                                for (auto it = solutions.begin(); it != solutions.end(); it++) {
                                    if(are_isomorphic((*it),d)) {
                                        newSol = false;
                                        break;
                                    }
                                }
                                if (newSol) {
                                    solutions.push_back(d); 
                                    std::cout << counter++ << std::endl;
                                }
                                goto BACKUP;
                            }
                        }
NEXT_COMBINATION:;
                 return 0;
                    }
                }
            }
        }
    }

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/maxQuasi/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/maxQuasi/" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }
    return 0;
}

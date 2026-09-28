#include "hds_quasiplanar.h"
#include "iso.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t n = 32;
const Edges edges = // expected: 32 * 5.25 = 168
{
    {0,1},{0,6},{0,7},{0,8},{0,9},{0,11},{0,15},{0,25},
    {1,2},{1,6},{1,7},{1,8},{1,9},{1,10},{1,11},{1,12},{1,24},{1,25},{1,26},{1,28},
    {2,3},{2,4},{2,5},{2,9},{2,11},{2,12},{2,13},{2,25},{2,26},{2,27},{2,28},{2,29},
    {3,4},{3,5},{3,12},{3,26},{3,28},{3,29},{3,30},
    {4,5},{4,11},{4,12},{4,13},{4,15},{4,29},
    {5,6},{5,8},{5,12},{5,13},{5,14},{5,15},{5,28},{5,29},{5,30},{5,24},
    {6,7},{6,8},{6,9},{6,13},{6,15},{6,29},{6,30},{6,31},{6,24},{6,25},
    {7,8},{7,24},{7,25},{7,26},{7,30},
    {8,9},{8,10},{8,11},{8,15},{8,16},{8,17},{8,18},{8,22},
    {9,10},{9,11},{9,17},
    {10,11},{10,16},{10,17},{10,18},{10,20},
    {11,12},{11,17},{11,18},{11,19},{11,20},{11,21},
    {12,13},{12,14},{12,15},{12,18},{12,20},{12,21},{12,22},
    {13,14},{13,15},{13,21},
    {14,15},{14,16},{14,20},{14,21},{14,22},
    {15,16},{15,17},{15,21},{15,22},{15,23},
    {16,17},{16,21},{16,22},{16,23},{16,24},{16,25},{16,26},{16,30},{16,31},
    {17,18},{17,19},{17,20},{17,25},{17,26},{17,27},{17,31},
    {18,19},{18,20},{18,26},
    {19,20},{19,25},{19,26},{19,27},{19,29},
    {20,21},{20,26},{20,27},{20,28},{20,29},{20,30},
    {21,22},{21,23},{21,27},{21,29},{21,30},{21,31},
    {22,23},{22,30},
    {23,25},{23,29},{23,30},{23,31},
    {24,25},{24,30},{24,31},
    {25,26},{25,30},{25,31},
    {26,27},{26,28},{26,29},
    {27,28},{27,29},
    {28,29},
    {29,30},
    {30,31},
};

const std::size_t klim = 10; // 2n-7 = 57...

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "k = " << klim << ", n = " << n << std::endl;
    assert(edges.size() == 168);
    // const Edges edges = generateCompleteGraph(n);
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(10000,1); // assume no more than 10000 unique drawings up to iso

    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);
    auto start_edge = edges.begin() + 1;

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
                    break;
                }
                d_ind++;
            }
            if (newSol) {
                solutions.push_back(d);
                std::cout << ++counter << std::endl;
            }
            goto BACKUP;
        }
    }
END:
    std::cout << "Found " << solutions.size() << " min crossing drawings in total." << std::endl;
    if(solutions.size() == 0) {
        return 0;
    }

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::string filename = "../quasiDrawings/hexagon/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/hexagon/" + std::to_string(idx) + ".graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }

    std::cout << "Found " << counter << " drawings in total." << std::endl;
    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++)
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    return 0;
}

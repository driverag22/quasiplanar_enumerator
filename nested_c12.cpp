#include "hds_kplanar.h"
#include "iso.h"
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

const std::size_t klim = 3;
const std::size_t C = 12;
const std::size_t n = 2 * C;

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "C" << C << ", k = " << klim << ", n = " << n << std::endl;
    
    const Edges edges = {
        // inner uncrossed cycle
        // {0,1,3},
        {1,2,  klim},
        {2,3,  klim},
        {3,4,  klim},
        {4,5,  klim},
        {5,6,  klim},
        {6,7,  klim},
        {7,8,  klim},
        {8,9,  klim},
        {9,10, klim},
        {10,11,klim},
        {11,0, klim},
        // one matching
        {0,C+9},
        // outer uncrossed cycle
        {C+9 ,C+10,klim},
        {C+10,C+11,klim},
        {C+11,C   ,klim},
        {C   ,C+1 ,klim},
        {C+1 ,C+2 ,klim},
        {C+2 ,C+3 ,klim},
        {C+3 ,C+4 ,klim},
        {C+4 ,C+5 ,klim},
        {C+5 ,C+6 ,klim},
        {C+6 ,C+7 ,klim},
        {C+7 ,C+8 ,klim},
        {C+8 ,C+9 ,klim},
        // rest of matching
        // EVEN: (i, C + ( (i+9) % 12)) == (i, C + ( (i-3) % 12))
        // ODD: (i, C + ( (i+3) % 12)) == (i, C + ( (i-9) % 12))
        {2,C+11},{4,C+1},{6,C+3},{8,C+5},{10,C+7},
        {1,C+4},{3,C+6},{5,C+8},{7,C+10},{9,C},{11,C+2},
    };
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(100,1); // assume no more than 100 unique drawings up to iso

    Drawing<klim> d(n);
    // Edge 0: (0,1)
    d.add_first_edge(0, 1, klim);

    auto start_edge = edges.begin();
    for (auto e = start_edge;;) {
        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];
        int pcr = (e->size() == 3) ? (*e)[2] : 0;
        HdsPath p = d.first_path(u, v, pcr);
        if (p.empty()) {
BACKUP:
            // no way to add uv -> do previous edges differently
            do {
                if (--e == start_edge) {
                    goto END;
                }

                u = (*e)[0];
                assert(u == d.edges.back().u);
                v = (*e)[1];
                assert(v == d.edges.back().v);
                pcr = (e->size() == 3) ? (*e)[2] : 0;

                p = d.edges.back().built;
                d.remove_edge();
            } while (!d.next_path(p, v, pcr));
        } 
        d.add_edge(p, v, pcr);
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
            if (newSol) solutions.push_back(d);
            goto BACKUP;
        }
    }
END:
    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0) {
        return 0;
    }

    std::size_t idx = 0;
    for (auto it = solutions.begin();it!=solutions.end();it++) {
        std::cout << "Drawing " << idx << std::endl;
        // if ((*it).is_drawing_extensible()) {
        //     std::cout << "is extensible!" << std::endl;
        // }
        std::string filename = "../quasiDrawings/nested_c" + std::to_string(C) + "/" + std::to_string(idx) + ".json";
        std::ofstream of_json(filename);
        nlohmann::ordered_json output_json = (*it).serialize_to_json();
        of_json << output_json.dump(4);
        of_json.close();

        std::string filename2 = "../quasiDrawings/nested_c" + std::to_string(C) + "/" + std::to_string(idx) + ".graphml";
        // std::string filename2 = "../quasiDrawings/nested_c" + std::to_string(C) + "/" + std::to_string(idx) + "_iso.graphml";
        std::ofstream of_graphml(filename2);
        (*it).graphml_output(of_graphml);
        of_graphml.close();
        idx++;
    }

    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;
}

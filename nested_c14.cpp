// #include "hds_quasiplanar.h"
#include "hds_kplanar.h"
#include "iso.h"
//#include <sstream>
// #include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

Edges generate3C14Graph() {
    Edges edges;

    // 1. inner C14 cycle (Vertices 0..13) - FIXED FIRST EDGES
    for (std::size_t i = 0; i < 14; ++i)
        edges.push_back({i, (i + 1) % 14});

    // 1b. inner C14 chords at distance 2
    // for (std::size_t i = 0; i < 14; ++i)
    //     edges.push_back({i, (i + 2) % 14});

    // 2. inner to middle connections at distance 3
    for (std::size_t i = 0; i < 14; ++i) {
        if (i % 2 == 0) edges.push_back({i, 14 + ((i + 3) % 14)});
        else            edges.push_back({i, 14 + ((i + 11) % 14)});
    }

    // 3. middle C14 cycle (Vertices 14..27)
    for (std::size_t i = 0; i < 14; ++i)
        edges.push_back({14 + i, 14 + ((i + 1) % 14)});

    // 3b. middle chords at dist 2
    // for (std::size_t i = 0; i < 14; ++i)
    //     edges.push_back({14 + i, 14 + ((i + 2) % 14)});

    // 4. middle to outer connections at distance 3
    // for (std::size_t i = 0; i < 14; ++i) {
    //     if (i % 2 == 0) edges.push_back({14 + i, 14 + ((i + 3) % 14)});
    //     else if (i % 2 == 1) edges.push_back({14 + i, 14 + ((i + 11) % 14)});
    // }

    // 5. outer C14 cycle (Vertices 28..41)
    // for (std::size_t i = 0; i < 14; ++i)
    //     edges.push_back({28 + i, 28 + ((i + 1) % 14)});

    // 5b. Outer C14 chords at distance 2
    // for (std::size_t i = 0; i < 14; ++i)
    //     edges.push_back({28 + i, 28 + ((i + 2) % 14)});

    return edges;
}

const std::size_t n = 28;
const std::size_t klim = 3;


// innermost and outermost C_{14} edges are uncrossed
inline bool is_protected_edge(std::size_t edge_idx) {
    return (edge_idx < 14) || (edge_idx >= 28 && edge_idx < 42);
}

bool path_crosses_protected(const HdsPath& p, std::size_t current_edge_idx) {
    // if the current edge being added is an inner or outer cycle edge, it must have 0 crossings
    if (is_protected_edge(current_edge_idx) && p.size() > 2) {
        return true;
    }

    // check if path 'p' crosses any protected edge along the way
    for (std::size_t i = 1; i + 1 < p.size(); ++i) {
        if (p[i] != nullptr && p[i]->edge != nullptr) {
            if (is_protected_edge(p[i]->edge->label)) {
                return true;
            }
        }
    }
    return false;
}

int main() {
    std::cout << "\n\n ===================================================== \n";
    std::cout << "c14, k = " << klim << ", n = " << n << std::endl;
    const Edges edges = generate3C14Graph();
    std::vector< Drawing<klim> > solutions;
    std::vector<std::size_t> d_cnt(30000,0); // assume no more than 10000 unique drawings up to iso

    Drawing<klim> d(n);
    // Edge 0: (0,1)
    d.add_first_edge(edges[0][0], edges[0][1]);
    // Edges 1..12: (1,2), (2,3), ..., (12,13)
    for (std::size_t i = 1; i < 14; ++i) {
        std::size_t u = edges[i][0];
        std::size_t v = edges[i][1];
        HdsPath p = d.first_path(u, v);
        while (!p.empty() && path_crosses_protected(p, i)) {
            if (!d.next_path(p, v)) {
                p.clear();
                break;
            }
        }
        d.add_edge(p, v);
    }

    auto start_edge = edges.begin() + 14;
    const size_t total_edges = edges.size() - 14;
    std::vector<int> path_choice(edges.size(), 0); 
    size_t max_depth_reached = 0;
    uint64_t total_iterations = 0;

    int counter = 0;
    for (auto e = start_edge;;) {
        ++total_iterations;
        size_t current_depth = std::distance(start_edge, e);
        size_t current_edge_idx = std::distance(edges.begin(), e);

        if (current_depth > max_depth_reached) max_depth_reached = current_depth;
        if (total_iterations % 1000 == 0) {
            std::cout << "\r\033[K Iter: " << total_iterations << " | Depth: " << current_depth << "/" << total_edges << " | Max Depth: " << max_depth_reached << " | Choices: [";
            for (size_t i = 0; i <= current_depth; ++i) std::cout << path_choice[i] << " ";
            std::cout << "]" << std::flush;
        }

        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];
        HdsPath p = d.first_path(u, v);

        // Advance until finding a path that does NOT cross protected edges
        while (!p.empty() && path_crosses_protected(p, current_edge_idx)) {
            if (!d.next_path(p, v)) {
                p.clear();
                break;
            }
        }

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

                current_edge_idx = std::distance(edges.begin(), e);

                bool found_valid = false;
                while (d.next_path(p, v)) {
                    if (!path_crosses_protected(p, current_edge_idx)) {
                        found_valid = true;
                        break;
                    }
                }

                if (found_valid) break;
            } while (true);

            size_t backtracked_depth = std::distance(start_edge, e);
            path_choice[backtracked_depth]++;
        } else {
            // reset choice count when entering a fresh edge branch
            path_choice[current_depth] = 0;
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
                // if (++counter % 25 == 0) std::cout << "\r Solutions found: " << counter << std::flush;
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
        std::cout << "Drawing " << idx++ << std::endl;
        if ((*it).is_drawing_extensible()) {
            std::cout << "is extensible!" << std::endl;
        }
    //     std::string filename = "../drawingsQuasi/K9/minGlobalCr/jsons/" + std::to_string(idx) + ".json";
    //     std::ofstream of_json(filename);
    //     nlohmann::ordered_json output_json = (*it).serialize_to_json();
    //     of_json << output_json.dump(4);
    //     of_json.close();

    //     // std::string filename2 = "drawings/K7_prop_test/minCr_" + std::to_string(idx) + ".graphml";
    //     // std::ofstream of_graphml(filename2);
    //     // (*it).graphml_output(of_graphml);
    //     // of_graphml.close();
    //     idx++;
    }

    std::cout << "Found " << counter << " drawings in total." << std::endl;
    std::cout << "Found " << solutions.size() << " unique drawings in total." << std::endl;

    for (std::size_t i = 0; i < solutions.size(); i++) {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;
}

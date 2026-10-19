#include "hds_quasiplanar.h"
#include "iso.h"
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include <cassert>
#include <set>
#include <string>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;
typedef std::pair<std::size_t, std::size_t> VPair;

// G + u: K10 on 0..9 with its vertices split into P (4), S (2), Q (4), and
//   c1 = 10 ~ Q,S,a1    a1 = 11 ~ P,S,c1    a2 = 12 ~ P,S    c2 = 13 ~ Q,S
//   u  = 14 ~ two vertices of P, one of S, one of Q, c1, c2
// Draws G + u in every possible way on top of every quasiplanar drawing of K10, and reports
//   - how many drawings of G + u each K10 drawing extends to,
//   - the drawings of G + u up to isomorphism, and which K10 drawings they come from,
//   - whether any drawing of G + u is extensible (G + u maximal).
const std::size_t n = 15;
const std::size_t klim = 23; // 2n - 7
const int num_drawings = 9;
const std::size_t c1 = 10, a1 = 11, a2 = 12, c2 = 13, u = 14;

// all k-subsets of the given vertices
std::vector<Edge> subsets(const Edge& from, std::size_t k) {
    std::vector<Edge> res;
    std::vector<bool> pick(from.size(), false);
    std::fill(pick.begin(), pick.begin() + k, true);
    do {
        Edge s;
        for (std::size_t i = 0; i < from.size(); ++i) if (pick[i]) s.push_back(from[i]);
        res.push_back(s);
    } while (std::prev_permutation(pick.begin(), pick.end()));
    return res;
}

// the vertices of from that are not in out (set difference, order kept)
Edge minus(const Edge& from, const Edge& out) {
    Edge res;
    for (std::size_t v : from) if (std::find(out.begin(), out.end(), v) == out.end()) res.push_back(v);
    return res;
}

// draws edges[i], edges[i+1], ... in all possible ways and calls done(d) for every complete drawing
void draw_edges(Drawing<klim>& d, const Edges& edges, std::size_t i, const std::function<void(Drawing<klim>&)>& done) {
    if (i == edges.size()) {
        done(d);
        return;
    }
    std::size_t x = edges[i][0];
    std::size_t y = edges[i][1];
    for (HdsPath p = d.first_path(x, y); !p.empty();) {
        d.add_edge(p, y);
        draw_edges(d, edges, i + 1, done);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, y)) break;
    }
}

// non-edges that can be inserted into d (d is unchanged afterwards)
std::vector<VPair> insertable_non_edges(Drawing<klim>& d) {
    std::set<VPair> present;
    for (const auto& e : d.edges) present.insert(std::minmax(e.u, e.v));
    std::vector<VPair> res;
    for (std::size_t x = 0; x < n; ++x) {
        for (std::size_t y = x + 1; y < n; ++y) {
            if (present.count({x, y})) continue;
            for (HdsPath p = d.first_path(x, y); !p.empty();) {
                d.add_edge(p, y);
                bool quasi = d.verify_quasiplanarity();
                p = d.edges.back().built;
                d.remove_edge();
                if (quasi) { res.push_back({x, y}); break; }
                if (!d.next_path(p, y)) break;
            }
        }
    }
    return res;
}

void save(const Drawing<klim>& d, const std::string& name) {
    std::ofstream(name + ".json") << d.serialize_to_json().dump(4);
    std::ofstream gml(name + ".graphml");
    d.graphml_output(gml);
}

int main() {
    const std::string dir = "../quasiDrawings/K10_all_quasi";
    std::cout << "K10 + c1,a1,a2,c2,u, k = " << klim << ", n = " << n << std::endl;

    std::size_t total = 0, extensible = 0;
    std::vector<Drawing<klim>> classes;            // one drawing of G + u per isomorphism class
    std::vector<std::size_t> class_size;           // number of drawings found in each class
    std::vector<std::set<int>> class_origin;       // the K10 drawings each class comes from
    const Edge k10 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    for (int i = 0; i < num_drawings; ++i) {
        std::ifstream input_file(dir + "/" + std::to_string(i) + ".json");
        if (!input_file.is_open()) {
            std::cerr << "Could not open drawing file " << i << std::endl;
            return 1;
        }
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        std::size_t found = 0;
        Drawing<klim> d(import_data, n);

        // every way of splitting the K10 into S, P, Q (the graph is symmetric in the K10 otherwise)
        for (const Edge& S : subsets(k10, 2)) {
            for (const Edge& P : subsets(minus(k10, S), 4)) {
                const Edge Q = minus(minus(k10, S), P);

                // ---- G: K10 + c1, a1, a2, c2 ----
                Edges g_edges;
                for (std::size_t v : Q) g_edges.push_back({v, c1});
                for (std::size_t v : S) g_edges.push_back({v, c1});
                for (std::size_t v : P) g_edges.push_back({v, a1});
                for (std::size_t v : S) g_edges.push_back({v, a1});
                g_edges.push_back({c1, a1});
                for (std::size_t v : P) g_edges.push_back({v, a2});
                for (std::size_t v : S) g_edges.push_back({v, a2});
                for (std::size_t v : Q) g_edges.push_back({v, c2});
                for (std::size_t v : S) g_edges.push_back({v, c2});

                draw_edges(d, g_edges, 0, [&](Drawing<klim>& dg) {
                    // ---- G + u: u ~ two vertices of P, one of S, one of Q, c1, c2 ----
                    for (const Edge& Pu : subsets(P, 2)) {
                        for (std::size_t su : S) {
                            for (std::size_t qu : Q) {
                                Edges u_edges;
                                for (std::size_t v : Pu) u_edges.push_back({v, u});
                                u_edges.push_back({su, u});
                                u_edges.push_back({qu, u});
                                u_edges.push_back({c1, u});
                                u_edges.push_back({c2, u});

                                draw_edges(dg, u_edges, 0, [&](Drawing<klim>& du) {
                                    assert(du.verify_quasiplanarity());
                                    ++total;
                                    ++found;
                                    if (!insertable_non_edges(du).empty()) ++extensible;

                                    std::size_t c = 0;
                                    while (c < classes.size() && !are_isomorphic(classes[c], du)) ++c;
                                    if (c == classes.size()) {
                                        classes.push_back(du);
                                        class_size.push_back(0);
                                        class_origin.push_back({});
                                    }
                                    ++class_size[c];
                                    class_origin[c].insert(i);
                                });
                            }
                        }
                    }
                });
            }
        }
        std::cout << "K10 drawing " << i << ": " << found << " drawings of G+u" << std::endl;
    }

    std::cout << "Drawings found: " << total << ", extensible: " << extensible << std::endl;
    std::cout << classes.size() << " drawings of G+u up to isomorphism:" << std::endl;
    for (std::size_t c = 0; c < classes.size(); ++c) {
        const std::string name = "../quasiDrawings/K10_plus4/u/k10_plus4_u_" + std::to_string(c + 1);
        save(classes[c], name);
        std::cout << "  class " << c + 1 << " (saved as " << name << "): " << class_size[c] << " drawings, from K10 drawing(s)";
        for (int i : class_origin[c]) std::cout << " " << i;
        std::cout << std::endl;
    }

    if (total == 0) std::cout << "G+u is not quasiplanar" << std::endl;
    else if (extensible == 0) std::cout << "G+u is maximal (6.5n-20 = " << (6.5 * n - 20) << ", m = 76)" << std::endl;
    else std::cout << "G+u is not maximal" << std::endl;

    return 0;
}

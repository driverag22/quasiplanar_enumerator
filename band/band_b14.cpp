#include "hds_quasiplanar.h"
#include "iso.h"
#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;
typedef std::vector<std::size_t> Vertices;   // vertex lists, and relabellings v -> sigma[v]
typedef std::pair<std::size_t, std::size_t> VPair;

const std::string k10_dir = "../quasiDrawings/K10_all_quasi/";        // the 9 drawings of K10
const std::string b14_dir = "../quasiDrawings/K10_plus4/band/B14/";   // output: the drawings of B_14 (for band_start)

// B_14: K10 on 0..9 plus four vertices,
//   10 ~ 0, 2, 3, 4, 5, 6, 11      11 ~ 1, 3, 4, 7, 8, 9, 10
//   12 ~ 1, 2, 3, 4, 7, 8          13 ~ 0, 2, 3, 5, 6, 9
// In the K10, the twins A = {0, 5, 6} are the neighbors of 10 and 13, the twins C = {1, 7, 8} those of 11 and 12, and
// 2, 3, 4, 9 are neighbors of 10, 12, 13 / of all four / of 10, 11, 12 / of 11, 13, resp.
//
// Draws B_14 in every possible way on top of every quasiplanar drawing of K10, and reports
//   - the automorphisms of B_14,
//   - the drawings of B_14 up to isomorphism (saved to b14_dir, in the labels above),
//   - whether any drawing of B_14 is extensible (B_14 maximal).
// Every drawing of B_14 restricts to a drawing of K10, i.e., to one of the 9 drawings of K10 with its vertices renamed.
// So on every K10 drawing, the program tries every way to give its vertices the roles above: which vertices form A and
// C, and which are 2, 3, 4, 9 (A and C are sets, as permuting twins is an automorphism). Each extra vertex is drawn as
// soon as its neighbors are known: 10 after choosing them, 12 after choosing A and 9, 11 after choosing 2, and 13
// after choosing 3. Run from the repository root.
const std::size_t n = 14;
const std::size_t klim = 2 * n - 7;
const int num_k10 = 9;
typedef Drawing<klim> D;

// the neighbors of 10, 11, 12, 13 with smaller labels
const std::vector<Vertices> extra_neighbors = {
    {0, 2, 3, 4, 5, 6},       // 10
    {1, 3, 4, 7, 8, 9, 10},   // 11
    {1, 2, 3, 4, 7, 8},       // 12
    {0, 2, 3, 5, 6, 9}        // 13
};

// all k-subsets of the given vertices
std::vector<Vertices> subsets(const Vertices& from, std::size_t k) {
    std::vector<Vertices> res;
    std::vector<bool> pick(from.size(), false);
    std::fill(pick.begin(), pick.begin() + k, true);
    do {
        Vertices s;
        for (std::size_t i = 0; i < from.size(); ++i) if (pick[i]) s.push_back(from[i]);
        res.push_back(s);
    } while (std::prev_permutation(pick.begin(), pick.end()));
    return res;
}

// the vertices of from that are not in out (set difference, order kept)
Vertices minus(const Vertices& from, const Vertices& out) {
    Vertices res;
    for (std::size_t v : from) if (std::find(out.begin(), out.end(), v) == out.end()) res.push_back(v);
    return res;
}

// the edges from z to the given vertices
Edges star(std::size_t z, const Vertices& nbrs) {
    Edges res;
    for (std::size_t w : nbrs) res.push_back({w, z});
    return res;
}

// all automorphisms of the graph with the given edges on vertices 0..nv-1 (backtracking over the images)
std::vector<Vertices> automorphisms(const std::set<VPair>& edges, std::size_t nv) {
    auto adjacent = [&](std::size_t a, std::size_t b) { return edges.count(std::minmax(a, b)) > 0; };
    std::vector<std::size_t> degree(nv, 0);
    for (auto [a, b] : edges) { ++degree[a]; ++degree[b]; }

    std::vector<Vertices> res; // relabeling i: v -> res[i][v]
    Vertices image(nv);
    std::vector<bool> used(nv, false);
    std::function<void(std::size_t)> extend = [&](std::size_t v) {
        if (v == nv) { // base case
            res.push_back(image);
            return;
        }
        // try v \mapsto w
        for (std::size_t w = 0; w < nv; ++w) {
            if (used[w] || degree[w] != degree[v]) continue;
            bool ok = true;
            for (std::size_t a = 0; a < v && ok; ++a) ok = adjacent(a, v) == adjacent(image[a], w);
            if (!ok) continue;
            image[v] = w;
            used[w] = true;
            extend(v + 1);
            used[w] = false;
        }
    };
    extend(0);
    return res;
}

// draws edges[i], edges[i+1], ... in all possible ways and calls done(d) for every complete drawing
void draw_edges(D& d, const Edges& edges, std::size_t i, const std::function<void(D&)>& done) {
    if (i == edges.size()) {
        done(d);
        return;
    }
    std::size_t a = edges[i][0];
    std::size_t b = edges[i][1];
    for (HdsPath p = d.first_path(a, b); !p.empty();) {
        d.add_edge(p, b);
        draw_edges(d, edges, i + 1, done);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, b)) break;
    }
}

// non-edges that can be inserted into d
std::vector<VPair> insertable_non_edges(D& d) {
    std::set<VPair> present;
    for (const auto& e : d.edges) present.insert(std::minmax(e.u, e.v));
    std::vector<VPair> res;
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = a + 1; b < n; ++b)
            if (!present.count({a, b}) && !d.first_path(a, b).empty()) res.push_back({a, b});
    return res;
}

// the drawing j with every vertex v renamed sigma[v]
nlohmann::json relabel(nlohmann::json j, const Vertices& sigma) {
    for (auto& st : j["drawing_recipe"]) {
        st["u"] = sigma[st["u"].get<std::size_t>()];
        st["v"] = sigma[st["v"].get<std::size_t>()];
    }
    return j;
}

// saves d as name.json and name.graphml, and fails loudly if a file cannot be written
void save(const D& d, const std::string& name) {
    std::ofstream out(name + ".json"), gml(name + ".graphml");
    out << d.serialize_to_json().dump(4);
    d.graphml_output(gml);
    out.close();
    gml.close();
    if (!out || !gml) throw std::runtime_error("could not write " + name + ".json and .graphml");
}

int main() {
    std::set<VPair> b14_edges;
    for (std::size_t a = 0; a < 10; ++a) // K10 edges
        for (std::size_t b = a + 1; b < 10; ++b) b14_edges.insert({a, b});
    for (std::size_t z = 10; z < n; ++z) // extra edges
        for (std::size_t w : extra_neighbors[z - 10]) b14_edges.insert({w, z});
    const std::vector<Vertices> aut = automorphisms(b14_edges, n);

    // automorphisms of B_14
    std::cout << "B_14 has " << b14_edges.size() << " edges and " << aut.size() << " automorphisms, moving only";
    std::set<std::size_t> moved;
    for (const Vertices& sigma : aut) for (std::size_t v = 0; v < n; ++v) if (sigma[v] != v) moved.insert(v);
    for (std::size_t v : moved) std::cout << " " << v;
    std::cout << std::endl;

    std::size_t total = 0, extensible = 0;
    std::vector<D> classes;              // one drawing of B_14 per isomorphism class (in the labels of its K10 drawing)
    std::vector<Vertices> class_role;    // class_role[c][v]: the vertex of classes[c] that plays vertex v of B_14
    const Vertices k10 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    for (int i = 0; i < num_k10; ++i) {
        std::ifstream in(k10_dir + std::to_string(i) + ".json");
        nlohmann::json j; in >> j; D d(j, n);
        std::size_t found = 0;

        for (const Vertices& n10 : subsets(k10, 6)) {                      
            draw_edges(d, star(10, n10), 0, [&](D& d10) {                               // draw edges incident on vertex 10 (n10)
                for (const Vertices& A : subsets(n10, 3)) {                
                    const Vertices B = minus(n10, A);                      
                    const Vertices rest = minus(k10, n10);                 
                    for (std::size_t p9 : rest) {
                        const Vertices C = minus(rest, {p9});              
                        Vertices n12 = C;
                        n12.insert(n12.end(), B.begin(), B.end());
                        draw_edges(d10, star(12, n12), 0, [&](D& d12) {                 // draw edges incident on vertex 12 (C \cup B)
                            for (std::size_t p2 : B) {
                                const Vertices B34 = minus(B, {p2});       
                                Vertices n11 = C;
                                n11.insert(n11.end(), {B34[0], B34[1], p9, 10});
                                draw_edges(d12, star(11, n11), 0, [&](D& d11) {         // draw edges incident on vertex 11 (C \cup {p3,p4,p9} \cup {10})
                                    for (std::size_t k = 0; k < 2; ++k) {               // two choices for {p3,p4} = B34 = B - {p2}
                                        const std::size_t p3 = B34[k], p4 = B34[1 - k];
                                        Vertices n13 = A;
                                        n13.insert(n13.end(), {p2, p3, p9});
                                        draw_edges(d11, star(13, n13), 0, [&](D& d13) { // Draw edges incident on vertex 13 (A \cup {p2,p3,p9})
                                            assert(d13.verify_quasiplanarity());
                                            ++total;
                                            ++found;
                                            if (!insertable_non_edges(d13).empty()) ++extensible;

                                            std::size_t c = 0;
                                            while (c < classes.size() && !are_isomorphic(classes[c], d13)) ++c;
                                            if (c == classes.size()) {
                                                classes.push_back(d13);
                                                class_role.push_back({A[0], C[0], p2, p3, p4, A[1], A[2], C[1], C[2], p9,
                                                                      10, 11, 12, 13});
                                            }
                                        });
                                    }
                                });
                            }
                        });
                    }
                }
            });
        }
        std::cout << "K10 drawing " << i << ": " << found << " drawings of B_14" << std::endl;
    }

    std::cout << "Drawings found: " << total << ", extensible: " << extensible << std::endl;
    for (std::size_t c = 0; c < classes.size(); ++c) {
        // rename the vertex that plays v to v
        Vertices to_b14(n);
        for (std::size_t v = 0; v < n; ++v) to_b14[class_role[c][v]] = v;
        save(D(relabel(classes[c].serialize_to_json(), to_b14), n), b14_dir + "B14_" + std::to_string(c + 1));
    }
    std::cout << classes.size() << " drawings of B_14 up to isomorphism, saved in " << b14_dir << std::endl;

    if (total == 0) std::cout << "B_14 is not quasiplanar" << std::endl;
    else if (extensible == 0) std::cout << "B_14 is maximal (6.5n-20 = " << (6.5 * n - 20) << ", m = " << b14_edges.size() << ")" << std::endl;
    else std::cout << "B_14 is not maximal" << std::endl;
    return 0;
}

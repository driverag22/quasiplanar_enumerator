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

// H = G + u + x, in the labels of the saved drawings of G + u (../quasiDrawings/K10_plus4/u, from k10_plus4_u.cpp):
//   K10 on 0..9, split into S = {0, 4}, P = {1, 2, 5, 6}, Q = {3, 7, 8, 9}
//   c1 = 10 ~ Q,S,a1    a1 = 11 ~ P,S,c1    a2 = 12 ~ P,S    c2 = 13 ~ Q,S
//   u  = 14 ~ 1, 2 (P), 0 (S), 3 (Q), c1, c2
//   x  = 15 ~ 1 (P, a neighbour of u), 5 (P, not a neighbour of u), 4 (S), 7 (Q), c1, a1
// G + u has four drawings and is maximal (k10_plus4_u.cpp). Every labelled drawing of G + u is therefore one of the
// four saved drawings, relabelled by an automorphism of G + u (up to homeomorphism). Drawing x on all of them in all
// possible ways gives every drawing of H. The program reports
//   - the automorphisms of G + u and of H,
//   - the drawings of H up to isomorphism (saved to ../quasiDrawings/K10_plus4/u_x),
//   - whether any drawing of H is extensible (H maximal),
//   - the edges that are uncrossed in every drawing of H.
const std::size_t n = 16;
const std::size_t klim = 25; // 2n - 7
const int num_drawings = 4;  // drawings of G + u
const std::size_t c1 = 10, a1 = 11, x = 15;
const Edge x_neighbours = {1, 5, 4, 7, c1, a1};

// all automorphisms of the graph with the given edges on vertices 0..nv-1 (backtracking over the images)
std::vector<Edge> automorphisms(const std::set<VPair>& edges, std::size_t nv) {
    auto adjacent = [&](std::size_t a, std::size_t b) { return edges.count(std::minmax(a, b)) > 0; };
    std::vector<std::size_t> degree(nv, 0);
    for (auto [a, b] : edges) { ++degree[a]; ++degree[b]; }

    std::vector<Edge> res;
    Edge image(nv);
    std::vector<bool> used(nv, false);
    std::function<void(std::size_t)> extend = [&](std::size_t v) {
        if (v == nv) {
            res.push_back(image);
            return;
        }
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

// the edges of a drawing, as pairs (smaller, larger)
std::set<VPair> edge_set(const Drawing<klim>& d) {
    std::set<VPair> res;
    for (const auto& e : d.edges) res.insert(std::minmax(e.u, e.v));
    return res;
}

// draws edges[i], edges[i+1], ... in all possible ways and calls done(d) for every complete drawing
void draw_edges(Drawing<klim>& d, const Edges& edges, std::size_t i, const std::function<void(Drawing<klim>&)>& done) {
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

// non-edges that can be inserted into d (d is unchanged afterwards)
std::vector<VPair> insertable_non_edges(Drawing<klim>& d) {
    std::set<VPair> present = edge_set(d);
    std::vector<VPair> res;
    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = a + 1; b < n; ++b) {
            if (present.count({a, b})) continue;
            if (!d.first_path(a, b).empty()) res.push_back({a, b});
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
    std::cout << "K10 + c1,a1,a2,c2,u,x, k = " << klim << ", n = " << n << std::endl;

    // load the drawings of G + u
    std::vector<nlohmann::json> gu;
    for (int t = 1; t <= num_drawings; ++t) {
        std::ifstream input_file("../quasiDrawings/K10_plus4/u/k10_plus4_u_" + std::to_string(t) + ".json");
        if (!input_file.is_open()) {
            std::cerr << "Could not open drawing file " << t << std::endl;
            return 1;
        }
        nlohmann::json j;
        input_file >> j;
        gu.push_back(j);
    }

    // automorphisms of G + u (all four drawings have the same graph)
    const std::set<VPair> gu_edges = edge_set(Drawing<klim>(gu[0], n));
    const std::vector<Edge> aut_gu = automorphisms(gu_edges, n - 1);
    std::cout << "G+u has " << gu_edges.size() << " edges and " << aut_gu.size() << " automorphisms" << std::endl;

    // automorphisms of H
    std::set<VPair> h_edges = gu_edges;
    for (std::size_t v : x_neighbours) h_edges.insert(std::minmax(v, x));
    const std::vector<Edge> aut_h = automorphisms(h_edges, n);
    std::cout << "H has " << h_edges.size() << " edges and " << aut_h.size() << " automorphisms:" << std::endl;
    for (const Edge& p : aut_h) {
        std::cout << "  ";
        for (std::size_t v = 0; v < n; ++v) if (p[v] != v) std::cout << v << "->" << p[v] << " ";
        std::cout << (p == Edge{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} ? "identity" : "") << std::endl;
    }

    std::size_t total = 0, extensible = 0;
    std::vector<Drawing<klim>> classes;      // one drawing of H per isomorphism class
    std::vector<std::size_t> class_size;     // number of drawings found in each class
    std::set<VPair> uncrossed = h_edges;     // edges uncrossed in every drawing found so far
    Edges x_edges;
    for (std::size_t v : x_neighbours) x_edges.push_back({v, x});

    for (int t = 0; t < num_drawings; ++t) {
        std::size_t found = 0;
        for (const Edge& sigma : aut_gu) {
            // the drawing sigma(G+u drawing t): its vertices renamed by sigma
            nlohmann::json j = gu[t];
            for (auto& st : j["drawing_recipe"]) {
                st["u"] = sigma[st["u"].get<std::size_t>()];
                st["v"] = sigma[st["v"].get<std::size_t>()];
            }
            Drawing<klim> d(j, n);

            // ---- H = G + u + x ----
            draw_edges(d, x_edges, 0, [&](Drawing<klim>& dh) {
                assert(dh.verify_quasiplanarity());
                ++total;
                ++found;
                if (!insertable_non_edges(dh).empty()) ++extensible;
                for (const auto& e : dh.edges) if (e.ncr > 0) uncrossed.erase(std::minmax(e.u, e.v));

                std::size_t c = 0;
                while (c < classes.size() && !are_isomorphic(classes[c], dh)) ++c;
                if (c == classes.size()) {
                    classes.push_back(dh);
                    class_size.push_back(0);
                }
                ++class_size[c];
            });
        }
        std::cout << "G+u drawing " << t + 1 << ": " << found << " drawings of H" << std::endl;
    }

    std::cout << "Drawings found: " << total << ", extensible: " << extensible << std::endl;
    std::cout << classes.size() << " drawings of H up to isomorphism:" << std::endl;
    for (std::size_t c = 0; c < classes.size(); ++c) {
        const std::string name = "../quasiDrawings/K10_plus4/u_x/k10_plus4_u_x_" + std::to_string(c + 1);
        save(classes[c], name);
        std::cout << "  class " << c + 1 << " (saved as " << name << "): " << class_size[c] << " drawings" << std::endl;
    }

    std::cout << "Edges uncrossed in every drawing of H:";
    for (auto [a, b] : uncrossed) std::cout << " " << a << "-" << b;
    std::cout << std::endl;

    if (total == 0) std::cout << "H is not quasiplanar" << std::endl;
    else if (extensible == 0) std::cout << "H is maximal (6.5n-20 = " << (6.5 * n - 20) << ", m = " << h_edges.size() << ")" << std::endl;
    else std::cout << "H is not maximal" << std::endl;

    return 0;
}

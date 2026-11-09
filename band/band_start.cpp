// drawings of B_42 have 238 edges, more than the default HDS_MAX_EDGES (96) of the header
#define HDS_MAX_EDGES 256
#include "hds_quasiplanar.h"
#include "iso.h"
#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <set>
#include <string>
#include <vector>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;
typedef std::vector<std::size_t> Vertices;   // vertex lists, and relabellings v -> sigma[v]
typedef std::pair<std::size_t, std::size_t> VPair;

const std::string b14_dir = "../quasiDrawings/K10_plus4/band/B14/";   // the drawings of B_14, up to isomorphism
const std::string b42_dir = "../quasiDrawings/K10_plus4/band/B42/";   // output: the drawings of B_42 (for band_window)

// The band graphs B_n: B_14 is the graph of the saved drawings in b14_dir (K10 on 0..9 plus the vertices 10..13), and
// every further vertex z is joined to the set N^-(z) of its neighbors with smaller labels (earlier_neighbors below),
// so B_{z+1} is B_z plus the vertex z.
// The program computes all drawings of B_15, ..., B_42 and checks that they are maximal:
//   - every labelled drawing of B_14 is one of its saved drawings, relabelled by an automorphism of B_14;
//   - every drawing of B_{z+1} restricts to a drawing of B_z, so drawing z in all possible ways on all drawings of
//     B_z gives all drawings of B_{z+1}.
// For every graph it reports the number of drawings found, how many are extensible, and how many there are up to
// isomorphism. All drawings of B_42 are saved in b42_dir. Run from the repository root.
const std::size_t n_start = 14;
const std::size_t n_end = 42;
// 2 * 42 - 7; for the smaller graphs a larger bound only allows more crossings than any quasiplanar drawing has
const std::size_t klim = 2 * n_end - 7;
typedef Drawing<klim> D;

// N^-(z) for z = 14, ..., 33
const std::vector<Vertices> start_neighbors = {
    {0, 1, 4, 5, 6, 11},         // 14
    {0, 4, 7, 10, 11, 14},       // 15
    {4, 8, 11, 12, 14, 15},      // 16
    {4, 10, 12, 14, 15, 16},     // 17
    {10, 11, 12, 15, 16, 17},    // 18
    {11, 12, 14, 16, 17, 18},    // 19
    {12, 14, 15, 16, 18, 19},    // 20
    {14, 15, 17, 18, 19, 20},    // 21
    {15, 16, 17, 19, 20, 21},    // 22
    {16, 17, 18, 20, 21, 22},    // 23
    {16, 18, 19, 21, 22, 23},    // 24
    {18, 19, 20, 22, 23, 24},    // 25
    {19, 20, 21, 22, 24, 25},    // 26
    {20, 21, 23, 24, 25, 26},    // 27
    {21, 22, 23, 25, 26, 27},    // 28
    {22, 24, 25, 26, 27, 28},    // 29
    {23, 24, 25, 27, 28, 29},    // 30
    {23, 26, 27, 28, 29, 30},    // 31
    {25, 26, 28, 29, 30, 31},    // 32
    {25, 27, 29, 30, 31, 32}     // 33
};

Vertices earlier_neighbors(std::size_t z) {
    if (z < 34) return start_neighbors[z - n_start];
    return {z - 7, z - 6, z - 4, z - 3, z - 2, z - 1};
}

// The automorphisms of B_14: the vertices 0, 5, 6 have the same neighbors, and so do 1, 7, 8; the 36 automorphisms
// permute these two triples and fix all other vertices.
std::vector<Vertices> automorphisms_b14() {
    std::vector<Vertices> res;
    Vertices t1 = {0, 5, 6}, t2 = {1, 7, 8};   // the images of 0, 5, 6 and of 1, 7, 8
    do {
        do {
            Vertices sigma(n_start);
            std::iota(sigma.begin(), sigma.end(), 0);
            sigma[0] = t1[0]; sigma[5] = t1[1]; sigma[6] = t1[2];
            sigma[1] = t2[0]; sigma[7] = t2[1]; sigma[8] = t2[2];
            res.push_back(sigma);
        } while (std::next_permutation(t2.begin(), t2.end()));
    } while (std::next_permutation(t1.begin(), t1.end()));
    return res;
}

// the edges of a drawing, as pairs (smaller, larger)
std::set<VPair> edge_set(const D& d) {
    std::set<VPair> res;
    for (const auto& e : d.edges) res.insert(std::minmax(e.u, e.v));
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

// non-edges among the vertices 0..nv-1 that can be inserted into d
std::vector<VPair> insertable_non_edges(D& d, std::size_t nv) {
    const std::set<VPair> present = edge_set(d);
    std::vector<VPair> res;
    for (std::size_t a = 0; a < nv; ++a)
        for (std::size_t b = a + 1; b < nv; ++b)
            if (!present.count({a, b}) && !d.first_path(a, b).empty()) res.push_back({a, b});
    return res;
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

// the drawings of B_{z+1} found by drawing z on the drawings of B_z
struct Step {
    std::vector<nlohmann::json> drawings;
    std::vector<D> classes;   // one drawing per isomorphism class
    std::size_t extensible = 0;
};

// draws z with its edges to N^-(z) in all possible ways on d (the first edge places z in a face)
void add_vertex(D& d, std::size_t z, Step& step) {
    Edges z_edges;
    for (std::size_t w : earlier_neighbors(z)) z_edges.push_back({w, z});
    draw_edges(d, z_edges, 0, [&](D& dz) {
        assert(dz.verify_quasiplanarity());
        if (!insertable_non_edges(dz, z + 1).empty()) ++step.extensible;
        nlohmann::json j = dz.serialize_to_json();
        step.drawings.push_back(j);
        // are_isomorphic needs a drawing without isolated vertices: keep only the vertices 0..z drawn so far
        j["num_vertices"] = z + 1;
        D drawn(j, z + 1);
        std::size_t c = 0;
        while (c < step.classes.size() && !are_isomorphic(step.classes[c], drawn)) ++c;
        if (c == step.classes.size()) step.classes.push_back(drawn);
    });
}

int main() {
    std::vector<nlohmann::json> b14;
    for (const auto& entry : std::filesystem::directory_iterator(b14_dir)) {
        if (entry.path().extension() != ".json") continue;
        std::ifstream in(entry.path());
        nlohmann::json j;
        in >> j;
        b14.push_back(j);
    }
    const std::vector<Vertices> aut = automorphisms_b14();
    std::cout << "B_14: " << b14.size() << " drawings up to isomorphism, " << aut.size() << " automorphisms" << std::endl;

    // B_15: the vertex 14 on every labelled drawing of B_14
    Step step;
    for (const auto& j : b14)
        for (const Vertices& sigma : aut) {
            D d(relabel(j, sigma), n_end);
            add_vertex(d, n_start, step);
        }
    // every B_n must have a quasiplane drawing, and none of its drawings may be extensible
    bool ok = (!step.drawings.empty() && (step.extensible == 0));
    std::cout << "B_15: " << step.drawings.size() << " drawings, " << step.extensible << " extensible, "
              << step.classes.size() << " up to isomorphism" << std::endl;

    // B_16, ..., B_42
    for (std::size_t z = n_start + 1; z < n_end; ++z) {
        Step next;
        // obtain B_{z+1} by adding z to each drawing of B_{z}
        for (const auto& j : step.drawings) {
            D d(j, n_end);
            add_vertex(d, z, next);
        }
        step = std::move(next);
        ok &= (!step.drawings.empty() && (step.extensible == 0));
        std::cout << "B_" << z + 1 << ": " << step.drawings.size() << " drawings, " << step.extensible << " extensible, "
                  << step.classes.size() << " up to isomorphism" << std::endl;
    }

    for (std::size_t i = 0; i < step.drawings.size(); ++i)
        save(D(step.drawings[i], n_end), b42_dir + "B42_" + std::to_string(i + 1));
    std::cout << (ok ? "all drawings of B_15, ..., B_42 are maximal" : "some B_n has no drawing, or some drawing is extensible")
              << std::endl;
    return ok ? 0 : 1;
}

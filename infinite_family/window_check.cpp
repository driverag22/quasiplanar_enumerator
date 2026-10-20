// drawings of H_3 have 260 edges, more than the default MAX_EDGES (96) of the header
#define HDS_MAX_EDGES 320
#include "hds_quasiplanar.h"
#include <filesystem>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <cassert>

typedef std::pair<std::size_t, std::size_t> VPair;

const std::string chain_dir = "../quasiDrawings/K10_plus4/chain/"; // drawings saved by chain_check (run from the repository root)

// Separation check on the drawings of three copies 0 -J- 1 -J- 2 of the chain, as saved by chain_check in
// chain_dir (file names with three parts, e.g. H1s0_H1s0m_H2s0.json). For every such drawing d:
//   1. containment: the faces of D(copy 1) used by side 0 (copy 0 and the edges M between copies 0 and 1) and the
//      faces used by side 2 (copy 2 and the edges M between copies 1 and 2) are disjoint, and the two shared edges
//      (copy 0 - copy 1 and copy 1 - copy 2) are uncrossed;
//   2. removing side 2 (its edges were drawn last) gives back the drawing of two copies that d was built on;
//   3. s-test: in that drawing of two copies, no curve from an unshared vertex of copy 0 to an isolated point s
//      ends in a face of D(copy 1) that side 2 used in step 1.
// D(copy 1) is the drawing of the edges of copy 1 alone (see copy_faces).
const std::size_t n1 = 16;                  // vertices of H
const std::size_t m1 = 82;                  // edges of H
const std::size_t max_copies = 3;
const std::size_t num_labels = n1 * max_copies;
const std::size_t num_vertices = n1 + (n1 - 2) * (max_copies - 1);
// 2 * 44 - 7, as in chain_check (the type must match the saved drawings). For the s-test on two copies plus s
// (31 vertices) the bound is 2 * 31 - 7 = 55; the larger bound only allows more curves, so a pass still counts.
const std::size_t klim = 2 * num_vertices - 7;

// the joint J: shared vertices (copy c vertex, copy c+1 vertex), as in chain_check
const std::vector<VPair> shared = {{5, 10}, {6, 0}};

typedef Drawing<klim> D;

// label of vertex y (H label) of copy c: a shared vertex keeps the label it has in copy c-1
std::size_t chain_label(std::size_t c, std::size_t y) {
    for (auto [prev, next] : shared)
        if (c > 0 && y == next) return chain_label(c - 1, prev);
    return c * n1 + y;
}

// determine if vertex v is in copy c
bool in_copy(std::size_t c, std::size_t v) {
    for (std::size_t y = 0; y < n1; ++y)
        if (chain_label(c, y) == v) return true;
    return false;
}

// ---- faces of the drawing of one copy ----
//
// d is stored planarized: crossings are extra (dummy) vertices, and the halfedges h, h->next, h->next->next, ...
// go once around a face of d. D(copy c) is d without the edges that are not edges of copy c (both endpoints in
// copy c). Deleting such an edge merges the faces on its two sides, so:
//
//   faces of D(copy c) = faces of d, where the two sides (h and h->twin) of every halfedge h of an edge not in
//                        copy c are joined (union-find).
//
// Every halfedge of d then lies in exactly one face of D(copy c): on its boundary if its edge is in copy c, inside it
// otherwise.
//
// To compare faces of D(copy c) in two different drawings that contain the same D(copy c) (here d with and without
// side 2), a face is named by the points of D(copy c) on its boundary: the vertices of copy c and the crossings
// between two edges of copy c. copy_faces checks that no two faces of D(copy c) get the same name, so the name
// identifies the face.

typedef std::set<std::string> FaceName;

struct CopyFaces {
    std::map<const HdsHalfedge*, std::size_t> face;  // the face of D(copy c) that each halfedge of d lies in
    std::vector<FaceName> name;                      // name of each face, identified by its vertices and crossings, faces numbered 0, 1, ...
};

std::size_t find_root(std::vector<std::size_t>& parent, std::size_t x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
}

std::set<FaceName> all_names(const CopyFaces& cf) { return std::set<FaceName>(cf.name.begin(), cf.name.end()); }

std::string edge_name(const HdsEdge* e) {
    return std::to_string(std::min(e->u, e->v)) + "-" + std::to_string(std::max(e->u, e->v));
}

CopyFaces copy_faces(const D& d, std::size_t c) {
    auto in_c = [&](const HdsEdge* e) { return in_copy(c, e->u) && in_copy(c, e->v); };

    // the faces of d: number the cycles h -> h->next -> ...
    std::map<const HdsHalfedge*, std::size_t> face_of_d;
    std::size_t num_faces_d = 0;
    for (const auto& h : d.halfedges) {
        if (face_of_d.count(&h)) continue;
        const HdsHalfedge* g = &h;
        do {
            face_of_d[g] = num_faces_d;
            g = g->next;
        } while (g != &h);
        ++num_faces_d;
    }

    // join the two sides of every halfedge of an edge not in copy c (union-find)
    std::vector<std::size_t> parent(num_faces_d);
    std::iota(parent.begin(), parent.end(), 0);
    for (const auto& h : d.halfedges)
        if (!in_c(h.edge)) 
            parent[find_root(parent, face_of_d[&h])] = find_root(parent, face_of_d[h.twin]);

    // number the faces of D(copy c) 0, 1, ...
    CopyFaces cf;
    std::map<std::size_t, std::size_t> number;   // root -> face number
    for (const auto& h : d.halfedges) {
        std::size_t r = find_root(parent, face_of_d[&h]);
        if (!number.count(r)) number[r] = number.size();
        cf.face[&h] = number[r];
    }
    cf.name.resize(number.size());

    // names: every halfedge h of copy c is on the boundary of its face, and so is its end point h->vertex, which
    // is a vertex of copy c, or a crossing of h's edge with the edge h->next continues on
    const std::size_t nv = d.vertices.size();
    for (const auto& h : d.halfedges) {
        if (!in_c(h.edge)) continue;
        std::string point;
        if (h.vertex->label < nv) {
            point = std::to_string(h.vertex->label);
        } else {
            const HdsEdge* other = h.next->edge;
            if (!in_c(other)) continue;          // crossing with an edge outside copy c: not a point of D(copy c)
            point = std::min(edge_name(h.edge), edge_name(other)) + "x" + std::max(edge_name(h.edge), edge_name(other));
        }
        cf.name[cf.face[&h]].insert(point);
    }
    // checks:
    // - the edges with both endpoints in copy c are exactly the m1 edges of copy c (no joint edge has both ends there)
    // - Euler for the connected plane graph D(copy c), planarized: V = n1 + X, E = m1 + 2X, so F = m1 - n1 + 2 + X,
    //   where X is the number of crossings between two edges of copy c (this checks the union-find)
    // - no two faces have the same name
    std::size_t edges_c = 0;
    for (const auto& e : d.edges) edges_c += in_c(&e);
    if (edges_c != m1) throw std::runtime_error("copy has a wrong number of edges");
    std::set<std::string> points;
    for (const FaceName& f : cf.name) points.insert(f.begin(), f.end());
    const std::size_t crossings = points.size() - n1;
    if (cf.name.size() != m1 - n1 + 2 + crossings) throw std::runtime_error("Euler's formula fails for D(copy)");
    std::set<FaceName> distinct(cf.name.begin(), cf.name.end());
    if (distinct.size() != cf.name.size()) throw std::runtime_error("two faces of D(copy) have the same name");
    return cf;
}

// ---- the three steps ----

// side of an edge outside copy 1: 0 or 2, the copy of its endpoint outside copy 1 (no edge joins copies 0 and 2)
std::size_t side(const HdsEdge& e) {
    std::size_t out = in_copy(1, e.u) ? e.v : e.u;
    return in_copy(2, out) ? 2 : 0;
}

bool in_side2(const HdsEdge& e) { return !(in_copy(1, e.u) && in_copy(1, e.v)) && side(e) == 2; }

D load(const std::string& file) {
    std::ifstream in(file);
    nlohmann::json j;
    in >> j;
    return D(j, num_labels);
}

int main() {
    std::vector<std::string> files;   // the drawings of three copies: names with three parts
    for (const auto& entry : std::filesystem::directory_iterator(chain_dir)) {
        std::string stem = entry.path().stem().string();
        if (std::count(stem.begin(), stem.end(), '_') == 2) files.push_back(stem);
    }
    std::sort(files.begin(), files.end());

    std::size_t contain_fails = 0, shared_fails = 0, s_fails = 0;
    for (const std::string& name : files) {
        D d = load(chain_dir + name + ".json");

        // 1. containment: the faces of D(copy 1) that the edges of each side lie in
        CopyFaces faces_full = copy_faces(d, 1);
        std::set<FaceName> used_by_side0, used_by_side2;
        for (const auto& h : d.halfedges) {
            if (in_copy(1, h.edge->u) && in_copy(1, h.edge->v)) continue;
            (side(*h.edge) == 2 ? used_by_side2 : used_by_side0).insert(faces_full.name[faces_full.face[&h]]);
        }
        std::size_t common = 0;
        for (const FaceName& f : used_by_side2) common += used_by_side0.count(f);
        if (common > 0) ++contain_fails;

        // the shared edges 5-6 of copy 0 (= 10-0 of copy 1) and 5-6 of copy 1 (= 10-0 of copy 2) are uncrossed
        std::size_t shared_crossings = 0;
        for (const auto& e : d.edges)
            for (std::size_t c = 0; c < 2; ++c)
                if (std::minmax(e.u, e.v) == std::minmax(chain_label(c, 5), chain_label(c, 6))) shared_crossings += e.ncr;
        if (shared_crossings > 0) ++shared_fails;

        // 2. remove side 2 (the last edges drawn) and compare with the drawing of two copies it was built on
        std::size_t removed = 0;
        while (!d.edges.empty() && in_side2(d.edges.back())) {
            d.remove_edge();
            ++removed;
        }
        for (const auto& e : d.edges)
            if (in_side2(e)) throw std::runtime_error("an edge of side 2 was not drawn last: " + name);
        const std::string base = name.substr(0, name.find_last_of('_'));
        nlohmann::json mine = nlohmann::json::parse(d.serialize_to_json()["drawing_recipe"].dump());
        std::ifstream in(chain_dir + base + ".json");
        nlohmann::json theirs;
        in >> theirs;
        if (mine != theirs["drawing_recipe"]) throw std::runtime_error("removing side 2 does not give " + base);
        // D(copy 1) is untouched, so its faces keep their names
        if (all_names(copy_faces(d, 1)) != all_names(faces_full)) throw std::runtime_error("face names changed: " + name);

        // 3. s-test: s = the first label with no edge (label 16, unused since copy 1's vertex 0 is labelled 6); draw every curve x-s from
        // an unshared vertex x of copy 0 and look up the face of D(copy 1) that s lies in
        std::size_t s = 0;
        while (s < num_labels && d.vertices[s].halfedge) ++s;
        assert(s < num_labels);
        std::size_t curves = 0, hits = 0;
        for (std::size_t y = 0; y < n1; ++y) {
            std::size_t x = chain_label(0, y);
            if (in_copy(1, x)) continue;
            // test all paths from x in copy 0 to new isolated vertex s
            for (HdsPath p = d.first_path(x, s); !p.empty();) {
                d.add_edge(p, s);
                ++curves;
                CopyFaces faces_curve = copy_faces(d, 1);
                if (all_names(faces_curve) != all_names(faces_full)) throw std::runtime_error("face names changed: " + name);
                // check if vertex s is placed in face of copy 2
                if (used_by_side2.count(faces_curve.name[faces_curve.face[d.vertices[s].halfedge]])) ++hits;
                p = d.edges.back().built;
                d.remove_edge(); // undo edge, go to next one
                if (!d.next_path(p, s)) break;
            }
        }
        // some edge from copy 0 was able to reach a face of copy 2
        if (hits > 0) ++s_fails;

        std::cout << name << ": D(copy 1) faces used by side 0: " << used_by_side0.size() << ", by side 2: " << used_by_side2.size()
                  << ", common " << common << ", crossings on shared edges " << shared_crossings << " | removed " << removed << " edges | s-test: " << curves
                  << " curves, ending in a face of side 2: " << hits << std::endl;
    }
    std::cout << "drawings of three copies: " << files.size() << ", containment fails " << contain_fails
              << ", crossed shared edges " << shared_fails << ", s-test fails " << s_fails << std::endl;
    return 0;
}

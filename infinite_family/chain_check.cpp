// drawings of H_3 have 260 edges, more than the default MAX_EDGES (96) of the header
#define HDS_MAX_EDGES 320
#include "hds_quasiplanar.h"
#include <filesystem>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <cassert>

typedef std::pair<std::size_t, std::size_t> VPair;

const std::string h_dir = "../quasiDrawings/K10_plus4/u_x/";     // the drawings of H
const std::string chain_dir = "../quasiDrawings/K10_plus4/chain/"; // output

// Chains of copies of H = G+u+x (16 vertices, 82 edges; its four drawings from k10_plus4_u_x.cpp).
// Consecutive copies are glued by the joint J: vertex 5 of copy c is vertex 10 of copy c+1, vertex 6 of copy c is
// vertex 0 of copy c+1 (so the uncrossed edge 5-6 of copy c is the uncrossed edge 10-0 of copy c+1), plus the edges M.
//
// The program enumerates every drawing of the chain with two and with three copies and checks that each one is
// saturated (no missing edge can be inserted):
//   - copy 1 runs through every labelled drawing of H: sigma(H_t) for the 4 drawings H_t and the 2 automorphisms
//     sigma of H (mirroring copy 1 = mirroring the whole drawing, which changes nothing);
//   - every further copy is drawn so that its own drawing is sigma(H_t) or its mirror image, for every t and sigma,
//     guided by the planar map of that target drawing (see add_copy).
// All drawings found are saved in chain_dir (input for window_check). Run from the repository root.
const std::size_t n1 = 16;                  // vertices of H
const std::size_t num_drawings = 4;         // drawings of H, up to its automorphisms and mirroring
const std::size_t max_copies = 3;
// labels: copy c (0-based) uses c * 16 + (0..15), see chain_label(); its two shared vertices take the labels of the previous
// copy, so 2 labels per further copy stay unused (isolated): three copies have 48 labels but 16 + 14 + 14 = 44 vertices
const std::size_t num_labels = n1 * max_copies;
const std::size_t num_vertices = n1 + (n1 - 2) * (max_copies - 1);
// 2 * 44 - 7; for two copies a larger bound than 2 * 30 - 7 only allows more crossings than any quasiplanar drawing
// has, so it changes no result
const std::size_t klim = 2 * num_vertices - 7;

// the joint J: shared vertices (copy c vertex, copy c+1 vertex) and the edges M (copy c vertex, copy c+1 vertex)
const std::vector<VPair> shared = {{5, 10}, {6, 0}};
const std::vector<VPair> M = {{7, 14}, {8, 3}, {8, 14}, {11, 2}, {11, 3}, {11, 4}, {15, 2}, {15, 14}};

// Aut(H) = {id, 8<->9}
const std::vector<std::vector<std::size_t>> automorphisms = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 10, 11, 12, 13, 14, 15}};

typedef Drawing<klim> D;

// label of vertex y (H label) of copy c: a shared vertex keeps the label it has in copy c-1
std::size_t chain_label(std::size_t c, std::size_t y) {
    for (auto [prev, next] : shared)
        if (c > 0 && y == next) return chain_label(c - 1, prev);
    return c * n1 + y;
}

// ---- planar map of a sub-drawing: everything needed to compare the drawing of one copy with a target ----
// For every kept edge a-b (a < b): the kept edges it crosses, in order from a to b, each with the direction in which
// it is crossed; for every vertex: the cyclic order of its kept edges. Two drawings of the same edges with equal maps
// are the same drawing (up to an orientation-preserving homeomorphism of the sphere).
struct Map {
    std::map<VPair, std::vector<std::pair<VPair, int>>> cross;
    std::map<std::size_t, std::vector<std::size_t>> rotation;
    bool operator==(const Map& o) const { return cross == o.cross && rotation == o.rotation; }
};

// cyclic order written from its smallest element
std::vector<std::size_t> normalize_cycle(std::vector<std::size_t> c) {
    std::rotate(c.begin(), std::min_element(c.begin(), c.end()), c.end());
    return c;
}

// the map of the edges of d with keep(u, v), written in the labels name(v)
Map extract_map(const D& d, const std::function<bool(std::size_t, std::size_t)>& keep,
                const std::function<std::size_t(std::size_t)>& name) {
    Map m;
    const std::size_t nv = d.vertices.size();
    // walk every kept edge from its endpoint u to v; "forward" halfedges point along that walk
    std::map<const HdsHalfedge*, bool> forward;
    std::map<const HdsEdge*, std::vector<const HdsHalfedge*>> walk;
    for (const auto& h : d.halfedges) {
        const HdsEdge* e = h.edge;
        if (!keep(e->u, e->v) || h.twin->vertex->label != e->u) continue;   // h is the segment leaving u
        std::vector<const HdsHalfedge*> segments;
        for (const HdsHalfedge* g = &h;; g = g->next->twin->next) {        // straight on through each crossing
            segments.push_back(g);
            forward[g] = true;
            forward[g->twin] = false;
            if (g->vertex->label < nv) break;                              // reached v (a real vertex)
        }
        walk[e] = segments;
    }
    for (auto& [e, segments] : walk) {
        std::size_t a = name(e->u), b = name(e->v);
        std::vector<std::pair<VPair, int>> list;
        for (std::size_t s = 0; s + 1 < segments.size(); ++s) {
            const HdsHalfedge* other = segments[s]->next;                  // leaves the crossing along the other edge
            const HdsEdge* f = other->edge;
            if (!keep(f->u, f->v)) continue;
            std::size_t c = name(f->u), dd = name(f->v);
            int direction = forward.at(other) ? 1 : 0;
            if (c > dd) { std::swap(c, dd); direction ^= 1; }
            list.push_back({{c, dd}, direction});
        }
        if (a > b) {                                                       // write it from the smaller end
            std::swap(a, b);
            std::reverse(list.begin(), list.end());
            for (auto& x : list) x.second ^= 1;
        }
        m.cross[{a, b}] = list;
    }
    for (std::size_t v = 0; v < nv; ++v) {
        const HdsHalfedge* start = d.vertices[v].halfedge;
        if (!start) continue;
        std::vector<std::size_t> cycle;
        const HdsHalfedge* h = start;
        do {
            std::size_t other = h->edge->u == v ? h->edge->v : h->edge->u;
            if (keep(h->edge->u, h->edge->v)) cycle.push_back(name(other));
            h = h->twin->prev;
        } while (h != start);
        if (!cycle.empty()) m.rotation[name(v)] = normalize_cycle(cycle);
    }
    return m;
}

// the mirror image of a map: every crossing direction flips, every rotation reverses
Map mirror(const Map& m) {
    Map r = m;
    for (auto& [e, list] : r.cross)
        for (auto& x : list) x.second ^= 1;
    for (auto& [v, cycle] : r.rotation) {
        std::reverse(cycle.begin(), cycle.end());
        cycle = normalize_cycle(cycle);
    }
    return r;
}

// ---- drawing edges ----

// Paths from first_path / next_path are quasiplanar by construction (find_crossing never lets the new edge cross two
// edges that cross each other), so no extra test is needed while drawing; every complete drawing is checked once
// with verify_quasiplanarity() (see add_copy).

// draws the edges es[i], es[i+1], ... in all possible ways and calls done(d) for every complete drawing
void draw_edges(D& d, const std::vector<VPair>& es, std::size_t i, const std::function<void(D&)>& done) {
    if (i == es.size()) {
        done(d);
        return;
    }
    auto [a, b] = es[i];
    for (HdsPath p = d.first_path(a, b); !p.empty();) {
        d.add_edge(p, b);
        draw_edges(d, es, i + 1, done);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, b)) break;
    }
}

// can a-b be inserted into d, i.e. is there any path for it?
bool insertable(D& d, std::size_t a, std::size_t b) {
    return !d.first_path(a, b).empty();
}

// ---- adding a copy ----

// Adds copy c (c >= 1) to the drawing d of copies 0..c-1 in every possible way and calls done for every result.
// The drawing of copy c alone must be one of the labelled drawings sigma(H_t) or its mirror image. For each choice,
// the edges of copy c are drawn one at a time (in a fixed order, each new vertex reached from vertices already
// drawn); after each edge the drawing of copy c so far must equal the same part of the target drawing. The edges M
// to copy c - 1 are drawn as soon as both endpoints exist. Every drawing of the chain restricts to such a target,
// and so is found.
void add_copy(D& d, std::size_t c, const std::vector<nlohmann::json>& targets,
              const std::function<void(D&, std::size_t t, std::size_t sigma, bool mirrored)>& done) {
    // get vertices and edges of copy c in chain labels
    std::vector<std::size_t> copy(n1);
    for (std::size_t y = 0; y < n1; ++y) copy[y] = chain_label(c, y);
    std::vector<VPair> joint;   // M between copies c-1 and c
    for (auto [x, y] : M) joint.push_back({chain_label(c - 1, x), chain_label(c, y)});

    std::set<VPair> copy_edges;   // edges of copy c in chain labels
    {
        D h(targets[0], n1);
        for (const auto& e : h.edges) copy_edges.insert(std::minmax(copy[e.u], copy[e.v]));
    }

    // fix order of the copy-c vertices and edges, the same for every target (all H_t have the same graph, and every
    // sigma maps it to itself): the shared vertices are drawn already (with the shared edge); then the vertex with
    // most edges to drawn vertices, together with those edges follow in order (the first of them places it)
    std::set<std::size_t> drawn;
    for (auto [prev, next] : shared) drawn.insert(copy[next]);
    std::map<VPair, std::size_t> step_of;   // edge -> number of edges of copy c drawn with it
    for (auto e : copy_edges) if (drawn.count(e.first) && drawn.count(e.second)) step_of[e] = 0;
    std::vector<VPair> steps;
    while (drawn.size() < n1) {
        std::size_t best = 0, best_count = 0;
        for (std::size_t v : copy) { // find vertex with most edges to drawn vertices
            if (drawn.count(v)) continue;
            std::size_t count = 0;
            for (std::size_t w : drawn) count += copy_edges.count(std::minmax(v, w));
            if (count > best_count) { best = v; best_count = count; }
        }
        if (best_count == 0) throw std::runtime_error("copy is not connected to the drawn vertices");
        for (std::size_t w : drawn)
            if (copy_edges.count(std::minmax(best, w))) {
                // add edges corresponding to best vertex
                steps.push_back({w, best});
                step_of[std::minmax(best, w)] = steps.size();
            }
        drawn.insert(best);
    }

    // iterate over the 4 drawings, the 2 automorphisms and mirrored or not
    for (std::size_t t = 0; t < num_drawings; ++t) {
        for (std::size_t s = 0; s < automorphisms.size(); ++s) {
            for (bool mir : {false, true}) {
                // target vertex v is copy c's vertex sigma(v)
                auto to_chain = [&](std::size_t v) { return copy[automorphisms[s][v]]; };
                D target(targets[t], n1);

                // target maps after each step (in chain labels, mirrored if needed)
                // after each step, we check with corresponding target map to see if drawing is valid or not
                std::vector<Map> want(steps.size() + 1);
                for (std::size_t k = 1; k <= steps.size(); ++k) {
                    Map m = extract_map(target,
                                        [&](std::size_t u, std::size_t v) { return step_of.at(std::minmax(to_chain(u), to_chain(v))) <= k; },
                                        to_chain);
                    want[k] = mir ? mirror(m) : m;
                }

                std::set<VPair> copy_drawn;   // copy-c edges drawn so far (the shared edge is drawn already)
                for (auto& [e, k] : step_of) if (k == 0) copy_drawn.insert(e);
                auto keep = [&](std::size_t u, std::size_t v) { return copy_drawn.count(std::minmax(u, v)) > 0; };
                auto same = [](std::size_t v) { return v; };

                // step k: draw steps[k], compare with the target, then the edges of M that have become drawable
                std::vector<bool> joint_drawn(joint.size(), false);
                std::function<void(std::size_t)> step = [&](std::size_t k) {
                    if (k == steps.size()) {
                        assert(d.verify_quasiplanarity());
                        done(d, t, s, mir);
                        return;
                    }
                    auto [a, b] = steps[k];
                    copy_drawn.insert(std::minmax(a, b));
                    // try all ways to draw edge (a,b)
                    for (HdsPath p = d.first_path(a, b); !p.empty();) {
                        d.add_edge(p, b);
                        // make sure the induced drawing of copy c matches the target
                        if (extract_map(d, keep, same) == want[k + 1]) {
                            std::vector<VPair> now;   // edges of M with both endpoints drawn, not drawn yet
                            std::vector<std::size_t> idx;
                            for (std::size_t i = 0; i < joint.size(); ++i) // check if new joint edges can be added
                                if (!joint_drawn[i] && d.vertices[joint[i].first].halfedge && d.vertices[joint[i].second].halfedge) {
                                    now.push_back(joint[i]);
                                    idx.push_back(i);
                                }
                            // add joint edges, if there are any
                            for (auto i : idx) joint_drawn[i] = true;
                            draw_edges(d, now, 0, [&](D&) { step(k + 1); });
                            for (auto i : idx) joint_drawn[i] = false;
                        }
                        p = d.edges.back().built;
                        d.remove_edge();
                        if (!d.next_path(p, b)) break;
                    }
                    copy_drawn.erase(std::minmax(a, b));
                };
                step(0);
            }
        }
    }
}

// the missing edges between copy c and the earlier copies that can be inserted into d
std::vector<VPair> insertable_pairs(D& d, std::size_t c) {
    std::set<std::size_t> earlier, vertices_c;
    for (std::size_t k = 0; k < c; ++k) for (std::size_t y = 0; y < n1; ++y) earlier.insert(chain_label(k, y));
    for (std::size_t y = 0; y < n1; ++y) 
        // ensure we dont take vertices shared with previous drawing twice
        if (!earlier.count(chain_label(c, y))) vertices_c.insert(chain_label(c, y));
    std::set<VPair> present;
    for (const auto& e : d.edges) present.insert(std::minmax(e.u, e.v));
    std::vector<VPair> res;
    for (std::size_t a : earlier)
        for (std::size_t b : vertices_c)
            if (!present.count({a, b}) && insertable(d, a, b)) res.push_back({a, b});
    return res;
}

void save(const D& d, const std::string& name) {
    std::ofstream(chain_dir + name + ".json") << d.serialize_to_json().dump(4);
}

int main() {
    std::filesystem::create_directories(chain_dir);
    std::vector<nlohmann::json> targets;
    for (std::size_t t = 1; t <= num_drawings; ++t) { // loading 4 drawings of H = G+u+x
        std::ifstream in(h_dir + "k10_plus4_u_x_" + std::to_string(t) + ".json");
        nlohmann::json j;
        in >> j;
        targets.push_back(j);
    }

    std::size_t two = 0, two_ext = 0, three = 0, three_ext = 0;
    // copy 0 runs through the 4 drawings of H, together with their automorphisms (id, 8 <-> 9)
    for (std::size_t t1 = 0; t1 < num_drawings; ++t1) {
        for (std::size_t s1 = 0; s1 < automorphisms.size(); ++s1) {
            // copy 0 = sigma(H_t1): the target drawing with its vertices renamed by sigma
            nlohmann::json j = targets[t1];
            for (auto& st : j["drawing_recipe"]) {
                st["u"] = automorphisms[s1][st["u"].get<std::size_t>()];
                st["v"] = automorphisms[s1][st["v"].get<std::size_t>()];
            }
            D d(j, num_labels);
            const std::string name1 = "H" + std::to_string(t1 + 1) + "s" + std::to_string(s1);
            std::size_t found2 = 0, found3 = 0;

            // add copy 2
            add_copy(d, 1, targets, [&](D& d2, std::size_t t2, std::size_t s2, bool m2) {
                ++two;
                ++found2;
                const std::string name2 = name1 + "_H" + std::to_string(t2 + 1) + "s" + std::to_string(s2) + (m2 ? "m" : "");
                save(d2, name2);
                if (!insertable_pairs(d2, 1).empty()) {
                    ++two_ext;
                    std::cout << "  extensible: " << name2 << std::endl;
                }
                // add copy 3
                add_copy(d2, 2, targets, [&](D& d3, std::size_t t3, std::size_t s3, bool m3) {
                    ++three;
                    ++found3;
                    const std::string name3 = name2 + "_H" + std::to_string(t3 + 1) + "s" + std::to_string(s3) + (m3 ? "m" : "");
                    save(d3, name3);
                    if (!insertable_pairs(d3, 2).empty()) {
                        ++three_ext;
                        std::cout << "  extensible: " << name3 << std::endl;
                    }
                });
            });
            std::cout << "copy 1 = " << name1 << ": " << found2 << " drawings of two copies, " << found3
                      << " of three copies" << std::endl;
        }
    }
    std::cout << "two copies: " << two << " drawings, extensible " << two_ext << std::endl;
    std::cout << "three copies: " << three << " drawings, extensible " << three_ext << std::endl;
    return 0;
}

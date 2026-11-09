// drawings of B_42 have 238 edges, more than the default HDS_MAX_EDGES (96) of the header
#define HDS_MAX_EDGES 256
#include "hds_quasiplanar.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

typedef std::pair<std::size_t, std::size_t> VPair;

const std::string b42_dir = "../quasiDrawings/K10_plus4/band/B42/";   // the drawings of B_42, from band_start
const std::size_t num_drawings = 2;   // band_start finds exactly two drawings of B_42 (Lemma Small cases (b))

// Window check for the band graphs B_n: from vertex 34 on, a new vertex z is joined to z-7, z-6, z-4, z-3, z-2, z-1.
// The window of a drawing of B_n is the drawing of the edges with an endpoint in n-8, ..., n-1 (their other endpoints
// are n-15, ..., n-1); its hidden faces are the faces that the other edges lie in.
// For every drawing of B_42, with its window W drawn alone and the new vertex z = 42, the program checks:
//   1. z can be drawn on W in exactly one way, and then z and its edges lie in no hidden face;
//   2. no vertex of W other than the neighbors of z can be joined to z;
//   3. s-test: no curve from z to an isolated vertex s ends in a hidden face (every vertex of B_42 outside W lies in one);
//   4. the window of B_43, which is part of W + z, is W shifted by one, and so are its hidden faces: those that contain
//      a hidden face of W or an edge of W that left the window.
// The checks only involve W and its hidden faces, so by induction they hold for every n >= 42 (as window_check.cpp
// does for the chain of copies of H). Run from the repository root.
const std::size_t n = 42;
const std::size_t z = n;                 // the new vertex
const std::size_t s = n + 1;             // an isolated vertex: the end of curves from z
const std::size_t num_labels = n + 2;
// the curves are drawn on W (48 edges) and z (6 edges); a curve crosses every edge at most once, so a bound of at
// least 54 finds them all, for every n; 2 * 44 - 7 also bounds the crossings of the edges of B_42
const std::size_t klim = 2 * num_labels - 7;
typedef Drawing<klim> D;

const std::vector<std::size_t> z_neighbors = {z - 7, z - 6, z - 4, z - 3, z - 2, z - 1};

// is u-v an edge of the window of B_m (both endpoints below m, one of them in m-8, ..., m-1)?
bool in_window(std::size_t m, std::size_t u, std::size_t v) {
    return std::max(u, v) < m && std::max(u, v) + 8 >= m;
}

// vertex labels counted from the first vertex m-15 of the window of B_m, so that windows can be compared
std::size_t window_label(std::size_t m, std::size_t v) {
    return v + 15 - m;
}

// ---- faces of the window ----
//
// As copy_faces in window_check.cpp: d is stored planarized, and deleting an edge outside the window merges the faces
// on its two sides, so the faces of the window of B_m inside d are the faces of d, joined across every halfedge of an
// edge outside the window (union-find). A face is named by the points of the window on its boundary: its vertices and
// the crossings of two of its edges, in window labels. No two faces may get the same name, so a name identifies a face.
typedef std::set<std::string> FaceName;

// struct to store the faces of window W
struct WindowFaces {
    std::map<const HdsHalfedge*, std::size_t> face; // face[h]: face of the window containing halfedge h of d
    std::vector<FaceName> name;                     // name[f]: the name (encoding) of window face f, as a sequence of strings
    std::set<std::size_t> met;                      // "met" window faces (containing an edge not in the window)
};

// union-find for merged faces
std::size_t find_root(std::vector<std::size_t>& parent, std::size_t x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
}

// compute window faces of d with vertex z=m
WindowFaces window_faces(const D& d, std::size_t m) {
    auto inside = [&](const HdsEdge* e) { return in_window(m, e->u, e->v); };

    // get faces of d: number the cycles h -> h->next -> ...
    std::map<const HdsHalfedge*, std::size_t> face_of_d;
    std::size_t num_faces_d = 0;
    for (const auto& h : d.halfedges) {
        if (face_of_d.count(&h)) continue;
        const HdsHalfedge* g = &h;
        do { // walk fully around face
            face_of_d[g] = num_faces_d;
            g = g->next;
        } while (g != &h);
        ++num_faces_d;
    }

    // join the two sides of every halfedge of an edge outside the window
    std::vector<std::size_t> parent(num_faces_d);
    std::iota(parent.begin(), parent.end(), 0);
    for (const auto& h : d.halfedges)
        if (!inside(h.edge)) 
            parent[find_root(parent, face_of_d[&h])] = find_root(parent, face_of_d[h.twin]);

    // number the faces of the window 0, 1, ...
    WindowFaces wf;
    std::map<std::size_t, std::size_t> number;   // root -> face number
    for (const auto& h : d.halfedges) {
        std::size_t r = find_root(parent, face_of_d[&h]);
        if (!number.count(r)) {
            std::size_t k = number.size();
            number[r] = k;
        }
        wf.face[&h] = number[r];
        if (!inside(h.edge)) wf.met.insert(number[r]);
    }
    wf.name.resize(number.size());

    // names: the end point of every halfedge h of the window is on the boundary of h's face;
    // it is a vertex, or a crossing of h and h->next
    auto edge_name = [&](const HdsEdge* e) {
        return std::to_string(window_label(m, std::min(e->u, e->v))) + "-" + std::to_string(window_label(m, std::max(e->u, e->v)));
    };
    const std::size_t nv = d.vertices.size();
    // build encodings (names) for each face of the window
    for (const auto& h : d.halfedges) {
        if (!inside(h.edge)) continue;
        std::string point;
        if (h.vertex->label < nv) // a vertex, simply add it
            point = std::to_string(window_label(m, h.vertex->label));
        else { // a crossing
            const HdsEdge* other = h.next->edge;
            if (!inside(other)) continue;   // crossing with an edge outside the window: not a point of the window
            // we add "x" to signal a crossing
            point = std::min(edge_name(h.edge), edge_name(other)) + "x" + std::max(edge_name(h.edge), edge_name(other));
        }
        wf.name[wf.face[&h]].insert(point);
    }

    // check 1: all names of faces of window (wf.name) are distint
    if (std::set<FaceName>(wf.name.begin(), wf.name.end()).size() != wf.name.size())
        throw std::runtime_error("two faces of a window have the same name");

    // check 2: Euler's formula for the planarized window 
    //    V = vertices + X, E = edges + 2X, so F = E - V + 2 = edges - vertices + 2 + X, 
    // with X the crossings of two window edges.
    std::set<std::size_t> vertices;
    std::set<std::string> points;
    std::size_t edges = 0;
    for (const auto& e : d.edges) // number of vertices and edges
        if (inside(&e)) { 
            ++edges; 
            vertices.insert(e.u); vertices.insert(e.v); 
        }
    for (const FaceName& f : wf.name) // number of vertices and crossings (V)
        points.insert(f.begin(), f.end());
    const std::size_t crossings = points.size() - vertices.size();
    if (wf.name.size() != edges - vertices.size() + 2 + crossings) throw std::runtime_error("Euler's formula fails for a window");

    return wf;
}

std::set<FaceName> all_names(const WindowFaces& wf) { return std::set<FaceName>(wf.name.begin(), wf.name.end()); }

// ---- planar map of the window (as in chain_check.cpp) ----
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
        for (std::size_t i = 0; i + 1 < segments.size(); ++i) {
            const HdsHalfedge* other = segments[i]->next;                  // leaves the crossing along the other edge
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

// ---- drawing ----

// draws edges[i], edges[i+1], ... in all possible ways and calls done(d) for every complete drawing
void draw_edges(D& d, const std::vector<VPair>& edges, std::size_t i, const std::function<void(D&)>& done) {
    if (i == edges.size()) {
        done(d);
        return;
    }
    auto [a, b] = edges[i];
    for (HdsPath p = d.first_path(a, b); !p.empty();) {
        d.add_edge(p, b);
        draw_edges(d, edges, i + 1, done);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, b)) break;
    }
}

// Draws the window of the drawing d of B_42 alone into w (empty): its edges one at a time, in an order in which every
// edge has an endpoint drawn already, each along the route for which the map of the edges drawn so far equals their
// map in d (as chain_check.cpp draws a copy along a target drawing).
void draw_window(const D& d, D& w) {
    std::vector<VPair> order;   // the window edges, every one with an endpoint in an earlier edge
    std::set<std::size_t> reached;
    std::vector<VPair> todo;
    for (const auto& e : d.edges) if (in_window(n, e.u, e.v)) todo.push_back({e.u, e.v});
    while (!todo.empty())
        for (std::size_t i = 0; i < todo.size(); ++i)
            if (order.empty() || reached.count(todo[i].first) || reached.count(todo[i].second)) {
                order.push_back(todo[i]);
                reached.insert(todo[i].first);
                reached.insert(todo[i].second);
                todo.erase(todo.begin() + i);
                break;
            }

    // keep track of drawn edges
    std::set<VPair> drawn = {std::minmax(order[0].first, order[0].second)};
    auto is_drawn = [&](std::size_t u, std::size_t v) { return drawn.count(std::minmax(u, v)) > 0; };
    // keep, name params for extract_map on w
    auto all = [](std::size_t, std::size_t) { return true; };
    auto same = [](std::size_t v) { return v; };
    w.add_first_edge(order[0].first, order[0].second);
    // recursive function: step(k) draws order[k], order[k+1], ..., 
    // returns true if it gets through all the edges
    std::function<bool(std::size_t)> step = [&](std::size_t k) {
        if (k == order.size()) return true;
        auto [a, b] = order[k];
        if (!w.vertices[a].halfedge) std::swap(a, b);   // draw the edge from its endpoint that is drawn already
        drawn.insert(std::minmax(a, b));
        const Map want = extract_map(d, is_drawn, same);
        // for each way of drawing edge k, recurse with step(k+1)
        for (HdsPath p = w.first_path(a, b); !p.empty();) {
            w.add_edge(p, b);
            if (extract_map(w, all, same) == want && step(k + 1)) return true;
            p = w.edges.back().built;
            w.remove_edge();
            if (!w.next_path(p, b)) break;
        }
        drawn.erase(std::minmax(a, b));
        return false;
    };
    if (!step(1)) throw std::runtime_error("could not draw the window alone");
}

int main() {
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(b42_dir))
        if (entry.path().extension() == ".json") files.push_back(entry.path().string());
    std::sort(files.begin(), files.end());
    // make sure we find expected number of drawings 
    if (files.size() != num_drawings)
        throw std::runtime_error("expected " + std::to_string(num_drawings) + " drawings of B_42 in " + b42_dir + ", found " +
                                 std::to_string(files.size()));

    std::size_t fails = 0;
    for (const std::string& file : files) {
        std::ifstream in(file);
        nlohmann::json j;
        in >> j;
        D d(j, num_labels);

        // W, its faces and its hidden faces, inside the drawing of B_42
        const WindowFaces faces_d = window_faces(d, n);
        std::set<FaceName> hidden;
        for (std::size_t f : faces_d.met) hidden.insert(faces_d.name[f]);

        // W alone, which must have the same faces
        D w(num_labels);
        draw_window(d, w);
        if (all_names(window_faces(w, n)) != all_names(faces_d)) throw std::runtime_error("W alone has other faces: " + file);

        std::vector<VPair> z_edges;
        for (std::size_t v : z_neighbors) z_edges.push_back({v, z});
        std::size_t placements = 0, in_hidden = 0, reaching = 0, curves = 0, curves_hidden = 0;
        bool shifted = false;
        draw_edges(w, z_edges, 0, [&](D& wz) {
            // Window Lemma (a), placements: ways to draw z on W (exp. 1)
            ++placements;
            // Window Lemma (a), in_hidden: number of hidden faces of W that z and its edges lie in (exp. 0)
            const WindowFaces faces_wz = window_faces(wz, n);
            for (std::size_t f : faces_wz.met) in_hidden += hidden.count(faces_wz.name[f]);

            // Window Lemma (b), reaching: number of vertices of W that z can be joined to (exp. 0)
            for (std::size_t v = n - 15; v < n; ++v)
                if (std::find(z_neighbors.begin(), z_neighbors.end(), v) == z_neighbors.end() && !wz.first_path(v, z).empty())
                    ++reaching;

            // Window Lemma (c), curves_hidden: number of curves z-s with s in a hidden face
            for (HdsPath p = wz.first_path(z, s); !p.empty();) {
                // s-test: draw every curve z-s and look up the face of W that s lies in)
                wz.add_edge(p, s);
                ++curves;
                const WindowFaces faces_s = window_faces(wz, n);
                if (hidden.count(faces_s.name[faces_s.face.at(wz.vertices[s].halfedge)])) ++curves_hidden;
                p = wz.edges.back().built;
                wz.remove_edge();
                if (!wz.next_path(p, s)) break;
            }

            // Window Lemma (d): the window of B_43: its hidden faces contain a hidden face of W (untouched by z, see 1) or an edge of
            // W that left the window; compare with W, both in window labels
            const WindowFaces faces_43 = window_faces(wz, n + 1);
            std::set<FaceName> hidden_43;
            for (const auto& h : wz.halfedges) {
                const std::size_t f = faces_43.face.at(&h);
                if (faces_43.met.count(f) || hidden.count(faces_wz.name[faces_wz.face.at(&h)])) hidden_43.insert(faces_43.name[f]);
            }
            auto window_42 = [](std::size_t u, std::size_t v) { return in_window(n, u, v); };
            auto window_43 = [](std::size_t u, std::size_t v) { return in_window(n + 1, u, v); };
            auto label_42 = [](std::size_t v) { return window_label(n, v); };
            auto label_43 = [](std::size_t v) { return window_label(n + 1, v); };
            shifted = extract_map(wz, window_43, label_43) == extract_map(d, window_42, label_42)
                      && all_names(faces_43) == all_names(faces_d) && hidden_43 == hidden;
        });

        const bool ok = placements == 1 && in_hidden == 0 && reaching == 0 && curves_hidden == 0 && shifted;
        fails += !ok;
        std::cout << file << ": W " << w.edges.size() << " edges, " << faces_d.name.size() << " faces, " << hidden.size()
                  << " hidden | 1. ways to draw z " << placements << ", in hidden faces " << in_hidden
                  << " | 2. other vertices reaching z " << reaching << " | 3. s-test: " << curves
                  << " curves, ending in a hidden face " << curves_hidden << " | 4. window of B_43 = W shifted: "
                  << (shifted ? "yes" : "no") << std::endl;
    }
    std::cout << "drawings of B_42: " << files.size() << ", window check fails " << fails << std::endl;
    return fails == 0 ? 0 : 1;
}

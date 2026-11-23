// drawings of B_60 have 6 * 60 - 14 = 346 edges, more than the default HDS_MAX_EDGES (96) of the header
#define HDS_MAX_EDGES 384
#include "hds_quasiplanar.h"
#include "iso.h"
#include "../validation/common.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

// Extension of band_start beyond B_42: starting from the two saved labelled drawings of B_42 (b42_dir), the vertices
// 42, 43, ..., 59 are drawn one at a time in all possible ways on every drawing found (exactly as band_start does),
// which gives all drawings of B_43, ..., B_60 (every drawing of B_{z+1} restricts to a drawing of B_z). For every n it
// reports the number of drawings, how many are extensible, and how many there are up to isomorphism, and it checks:
//   - swap: the drawings are two, and they are related by the relabelling exchanging 5 and 6 (orientation-preserving,
//     label-preserving comparison by Canonizer::labelled of validation/common.h);
//   - regular: the window of every drawing of B_n (induced by n-10, ..., n-1), its cells and its hidden cells (cells
//     met by the rest of B_n), computed directly from the full drawing and written in window labels v -> v - (n - 10)
//     (i.e. relabelled by v -> v - (n - 42)), coincide with those of the drawing of B_42 it descends from;
//   - F^+: the hidden cells of W_{n+1}, computed directly, are the cells of W_{n+1} that contain a hidden cell of W_n
//     or meet an edge of W_n \ W_{n+1} (Lemma 11(d)), both inside the full drawing of B_{n+1}.
// The window code (cells by union-find, names, map) is that of band_window.cpp. Run from the repository root.

typedef std::vector<std::size_t> Vertices;
typedef std::pair<std::size_t, std::size_t> VPair;

const std::string b42_dir = "../quasiDrawings/K10_plus4/band/B42/";   // the drawings of B_42, from band_start
const std::size_t n_start = 42;
const std::size_t n_end = 60;
const std::size_t klim = 2 * n_end - 7;
typedef Drawing<klim> D;

Vertices earlier_neighbors(std::size_t z) { return {z - 7, z - 6, z - 4, z - 3, z - 2, z - 1}; }

// ---- as in band_start.cpp ----

std::set<VPair> edge_set(const D& d) {
    std::set<VPair> res;
    for (const auto& e : d.edges) res.insert(std::minmax(e.u, e.v));
    return res;
}

nlohmann::json relabel(nlohmann::json j, const Vertices& sigma) {
    for (auto& st : j["drawing_recipe"]) {
        st["u"] = sigma[st["u"].get<std::size_t>()];
        st["v"] = sigma[st["v"].get<std::size_t>()];
    }
    return j;
}

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

std::vector<VPair> insertable_non_edges(D& d, std::size_t nv) {
    const std::set<VPair> present = edge_set(d);
    std::vector<VPair> res;
    for (std::size_t a = 0; a < nv; ++a)
        for (std::size_t b = a + 1; b < nv; ++b)
            if (!present.count({a, b}) && !d.first_path(a, b).empty()) res.push_back({a, b});
    return res;
}

// ---- as in band_window.cpp (window of B_m: induced by m-10, ..., m-1) ----

const std::size_t window_size = 10;

bool in_window(std::size_t m, std::size_t u, std::size_t v) {
    return std::max(u, v) < m && std::min(u, v) + window_size >= m;
}

std::size_t window_label(std::size_t m, std::size_t v) { return v + window_size - m; }

typedef std::set<std::string> FaceName;

struct WindowFaces {
    std::map<const HdsHalfedge*, std::size_t> face;
    std::vector<FaceName> name;
    std::set<std::size_t> met;
};

std::size_t find_root(std::vector<std::size_t>& parent, std::size_t x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
}

WindowFaces window_faces(const D& d, std::size_t m) {
    auto inside = [&](const HdsEdge* e) { return in_window(m, e->u, e->v); };
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
    std::vector<std::size_t> parent(num_faces_d);
    std::iota(parent.begin(), parent.end(), 0);
    for (const auto& h : d.halfedges)
        if (!inside(h.edge)) parent[find_root(parent, face_of_d[&h])] = find_root(parent, face_of_d[h.twin]);
    WindowFaces wf;
    std::map<std::size_t, std::size_t> number;
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
    auto edge_name = [&](const HdsEdge* e) {
        return std::to_string(window_label(m, std::min(e->u, e->v))) + "-" + std::to_string(window_label(m, std::max(e->u, e->v)));
    };
    const std::size_t nv = d.vertices.size();
    for (const auto& h : d.halfedges) {
        if (!inside(h.edge)) continue;
        std::string point;
        if (h.vertex->label < nv)
            point = std::to_string(window_label(m, h.vertex->label));
        else {
            const HdsEdge* other = h.next->edge;
            if (!inside(other)) continue;
            point = std::min(edge_name(h.edge), edge_name(other)) + "x" + std::max(edge_name(h.edge), edge_name(other));
        }
        wf.name[wf.face[&h]].insert(point);
    }
    if (std::set<FaceName>(wf.name.begin(), wf.name.end()).size() != wf.name.size())
        throw std::runtime_error("two faces of a window have the same name");
    std::set<std::size_t> vertices;
    std::set<std::string> points;
    std::size_t edges = 0;
    for (const auto& e : d.edges)
        if (inside(&e)) {
            ++edges;
            vertices.insert(e.u);
            vertices.insert(e.v);
        }
    for (const FaceName& f : wf.name) points.insert(f.begin(), f.end());
    const std::size_t crossings = points.size() - vertices.size();
    if (wf.name.size() != edges - vertices.size() + 2 + crossings) throw std::runtime_error("Euler's formula fails for a window");
    return wf;
}

std::set<FaceName> all_names(const WindowFaces& wf) { return std::set<FaceName>(wf.name.begin(), wf.name.end()); }

struct Map {
    std::map<VPair, std::vector<std::pair<VPair, int>>> cross;
    std::map<std::size_t, std::vector<std::size_t>> rotation;
    bool operator==(const Map& o) const { return cross == o.cross && rotation == o.rotation; }
};

std::vector<std::size_t> normalize_cycle(std::vector<std::size_t> c) {
    std::rotate(c.begin(), std::min_element(c.begin(), c.end()), c.end());
    return c;
}

Map extract_map(const D& d, const std::function<bool(std::size_t, std::size_t)>& keep,
                const std::function<std::size_t(std::size_t)>& name) {
    Map m;
    const std::size_t nv = d.vertices.size();
    std::map<const HdsHalfedge*, bool> forward;
    std::map<const HdsEdge*, std::vector<const HdsHalfedge*>> walk;
    for (const auto& h : d.halfedges) {
        const HdsEdge* e = h.edge;
        if (!keep(e->u, e->v) || h.twin->vertex->label != e->u) continue;
        std::vector<const HdsHalfedge*> segments;
        for (const HdsHalfedge* g = &h;; g = g->next->twin->next) {
            segments.push_back(g);
            forward[g] = true;
            forward[g->twin] = false;
            if (g->vertex->label < nv) break;
        }
        walk[e] = segments;
    }
    for (auto& [e, segments] : walk) {
        std::size_t a = name(e->u), b = name(e->v);
        std::vector<std::pair<VPair, int>> list;
        for (std::size_t i = 0; i + 1 < segments.size(); ++i) {
            const HdsHalfedge* other = segments[i]->next;
            const HdsEdge* f = other->edge;
            if (!keep(f->u, f->v)) continue;
            std::size_t c = name(f->u), dd = name(f->v);
            int direction = forward.at(other) ? 1 : 0;
            if (c > dd) { std::swap(c, dd); direction ^= 1; }
            list.push_back({{c, dd}, direction});
        }
        if (a > b) {
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

// ---- windows of a drawing of B_m, in window labels ----

struct WindowData {
    Map map;
    std::set<FaceName> cells, hidden;
};

WindowData window_data(const D& d, std::size_t m) {
    WindowData w;
    w.map = extract_map(d, [m](std::size_t u, std::size_t v) { return in_window(m, u, v); },
                        [m](std::size_t v) { return window_label(m, v); });
    const WindowFaces wf = window_faces(d, m);
    w.cells = all_names(wf);
    for (std::size_t f : wf.met) w.hidden.insert(wf.name[f]);
    return w;
}

// F^+ inside the drawing d of B_{m+1}: the cells of W_{m+1} that contain a hidden cell of W_m (hidden_m, from the
// drawing of B_m) or meet an edge of W_m \ W_{m+1}
std::set<FaceName> f_plus(const D& d, std::size_t m, const std::set<FaceName>& hidden_m) {
    const WindowFaces fa = window_faces(d, m), fb = window_faces(d, m + 1);
    std::set<FaceName> res;
    for (const auto& h : d.halfedges) {
        const bool left = in_window(m, h.edge->u, h.edge->v) && !in_window(m + 1, h.edge->u, h.edge->v);
        if (left || hidden_m.count(fa.name[fa.face.at(&h)])) res.insert(fb.name[fb.face.at(&h)]);
    }
    return res;
}

// ---- one step ----

struct Item {
    nlohmann::json j;
    std::size_t root;              // index of the drawing of B_42 it descends from
    std::set<FaceName> hidden;     // hidden cells of its window, computed directly
};

struct Step {
    std::vector<Item> drawings;
    std::vector<D> classes;
    std::size_t extensible = 0;
};

void add_vertex(D& d, std::size_t z, const Item& parent, Step& step) {
    Edges z_edges;
    for (std::size_t w : earlier_neighbors(z)) z_edges.push_back({w, z});
    draw_edges(d, z_edges, 0, [&](D& dz) {
        assert(dz.verify_quasiplanarity());
        if (!insertable_non_edges(dz, z + 1).empty()) ++step.extensible;
        nlohmann::json j = dz.serialize_to_json();
        step.drawings.push_back({j, parent.root, {}});
        j["num_vertices"] = z + 1;
        D drawn(j, z + 1);
        std::size_t c = 0;
        while (c < step.classes.size() && !are_isomorphic(step.classes[c], drawn)) ++c;
        if (c == step.classes.size()) step.classes.push_back(drawn);
    });
}

double secs(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

int main() {
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(b42_dir))
        if (entry.path().extension() == ".json") files.push_back(entry.path().string());
    std::sort(files.begin(), files.end());

    Canonizer<D> cz;
    Vertices swap56(n_end);
    std::iota(swap56.begin(), swap56.end(), 0);
    std::swap(swap56[5], swap56[6]);
    // are the drawings two, and is the second the first with 5 and 6 exchanged (labelled, orientation preserved)?
    auto swap_check = [&](const std::vector<Item>& ds) {
        if (ds.size() != 2) return std::string("not two drawings");
        const Code c0 = cz.labelled(D(ds[0].j, n_end)), c1 = cz.labelled(D(ds[1].j, n_end));
        const Code s0 = cz.labelled(D(relabel(ds[0].j, swap56), n_end));
        const Code s1 = cz.labelled(D(relabel(ds[1].j, swap56), n_end));
        std::string r = (s0 == c1 && s1 == c0) ? "yes" : "NO";
        if (c0 == c1) r += " (labelled-equal drawings!)";
        return r;
    };

    // B_42
    Step step;
    std::vector<WindowData> root_window;
    for (const std::string& file : files) {
        std::ifstream in(file);
        nlohmann::json j;
        in >> j;
        D d(j, n_end);
        root_window.push_back(window_data(d, n_start));
        step.drawings.push_back({d.serialize_to_json(), root_window.size() - 1, root_window.back().hidden});
        std::cout << "B_42 drawing " << root_window.size() << " (" << file << "): " << d.edges.size() << " edges, "
                  << d.crossings.size() << " crossings, window " << root_window.back().cells.size() << " cells, "
                  << root_window.back().hidden.size() << " hidden" << std::endl;
    }
    std::cout << "B_42: " << step.drawings.size() << " drawings; 5<->6 swap: " << swap_check(step.drawings) << std::endl;
    bool ok = step.drawings.size() == 2;

    const auto t_all = std::chrono::steady_clock::now();
    for (std::size_t z = n_start; z < n_end; ++z) {
        const auto t0 = std::chrono::steady_clock::now();
        Step next;
        for (const Item& it : step.drawings) {
            D d(it.j, n_end);
            add_vertex(d, z, it, next);
        }
        const double t_draw = secs(t0);
        // regularity and F^+ for every drawing of B_{z+1}
        std::size_t regular = 0, fplus_ok = 0, crossings = 0, edges = 0;
        for (Item& it : next.drawings) {
            D d(it.j, n_end);
            crossings = d.crossings.size();
            edges = d.edges.size();
            const WindowData w = window_data(d, z + 1);
            it.hidden = w.hidden;
            const WindowData& r = root_window[it.root];
            regular += (w.map == r.map && w.cells == r.cells && w.hidden == r.hidden);
            const Item* parent = nullptr;   // the drawing of B_z it was obtained from: same root
            for (const Item& p : step.drawings)
                if (p.root == it.root) parent = &p;
            fplus_ok += (f_plus(d, z, parent->hidden) == w.hidden);
        }
        std::set<std::size_t> roots;
        for (const Item& it : next.drawings) roots.insert(it.root);
        const std::string sw = swap_check(next.drawings);
        step = std::move(next);
        const bool good = step.drawings.size() == 2 && roots.size() == 2 && step.extensible == 0 && step.classes.size() == 1 &&
                          sw == "yes" && regular == step.drawings.size() && fplus_ok == step.drawings.size();
        ok &= good;
        std::cout << "B_" << z + 1 << ": " << step.drawings.size() << " drawings (" << edges << " edges, " << crossings
                  << " crossings), " << step.extensible << " extensible, " << step.classes.size()
                  << " up to isomorphism | 5<->6 swap: " << sw << " | regular " << regular << "/" << step.drawings.size()
                  << " | F^+ = hidden " << fplus_ok << "/" << step.drawings.size() << " | " << t_draw << " s drawing, "
                  << secs(t0) << " s step, " << secs(t_all) << " s total, peak " << peak_mem_mb() << " MB"
                  << (good ? "" : "  <-- UNEXPECTED") << std::endl;
    }
    std::cout << (ok ? "all checks passed for B_43, ..., B_60" : "some check failed") << std::endl;
    return ok ? 0 : 1;
}

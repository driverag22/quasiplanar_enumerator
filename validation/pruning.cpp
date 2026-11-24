// Empirical check of the pruning test may_complete of dodecahedron/dodecahedron.cpp (k = 3).
//
// Build: g++ -O3 -Wall -std=c++17 -Iinclude validation/pruning.cpp -o validation/.bin/pruning
// Run:   validation/.bin/pruning a <instance>        (a) same drawings with and without pruning
//            instances: K7free (K7, only edge 01 fixed), K8 (star at 0 fixed), K44 (star 0->B fixed),
//                       K9star3, K9tri, K9path3 (K9 minus 3 edges, star at 0 fixed)
//        validation/.bin/pruning b <instance>        (b) look-ahead edges (as in step 2 of dodecahedron.cpp)
//            instances: K8v7 (G = K8, G' = G minus the edges at 7), K44v3 (G = K4,4, G' = G minus the edges at
//                       vertex 3 of side A), K7pend (G = K7 + vertex 7 ~ 0,1,2 + vertex 8 ~ 3,7; G' = K7)
//        validation/.bin/pruning c <seconds>         (c) frame F of dodecahedron.cpp without pruning (time-limited)
//
// may_complete, ball, faces_at, Faces and draw_edges are copied from dodecahedron/dodecahedron.cpp; draw_edges gets
// a switch to turn the pruning off and counts search nodes (edges added).
#include "hds_kplanar.h"
#include "iso.h"
#include "common.h"
#include <set>

const std::size_t klim = 3;
typedef Drawing<klim> D;
typedef std::vector<std::size_t> Vertices;
typedef std::vector<bool> FaceSet;

// ---------------- copied from dodecahedron/dodecahedron.cpp (unchanged) ----------------
bool is_drawn(const D& d, std::size_t v) { return d.vertices[v].halfedge != nullptr; }

struct Faces {
    std::vector<std::size_t> face;
    std::vector<std::vector<const HdsHalfedge*>> doors;

    explicit Faces(const D& d) : face(d.halfedges.size(), SIZE_MAX) {
        for (const auto& h : d.halfedges) {
            if (face[h.label] != SIZE_MAX) continue; // face already processed
            doors.emplace_back();
            const HdsHalfedge* g = &h;
            do { // walk around face, setting labels and storing doors
                face[g->label] = doors.size() - 1;
                if (g->edge->ncr < klim) doors.back().push_back(g);
                g = g->next;
            } while (g != &h);
        }
    }
    std::size_t size() const { return doors.size(); }
};

FaceSet faces_at(const D& d, const Faces& faces, std::size_t v) {
    FaceSet res(faces.size(), false);
    const HdsHalfedge* h = d.vertices[v].halfedge;   // points to v
    do { // walk around halfedges pointing to v
        res[faces.face[h->label]] = true;
        h = h->next->twin;
    } while (h != d.vertices[v].halfedge);
    return res;
}

FaceSet ball(const Faces& faces, const FaceSet& start, std::size_t steps, std::size_t avoid = SIZE_MAX) {
    FaceSet reached = start;
    std::vector<std::size_t> layer;   // the faces reached in the last step
    for (std::size_t f = 0; f < faces.size(); ++f)
        if (start[f]) layer.push_back(f);
    for (std::size_t s = 0; s < steps && !layer.empty(); ++s) {
        std::vector<std::size_t> next;
        for (std::size_t f : layer) // take every "door" halfedge on face f
            for (const HdsHalfedge* h : faces.doors[f]) {
                const std::size_t g = faces.face[h->twin->label];
                if (reached[g] || h->edge->u == avoid || h->edge->v == avoid) continue;
                reached[g] = true;
                next.push_back(g);
            }
        layer = next;
    }
    return reached;
}

bool may_complete(D& d, const Edges& edges, std::size_t i) {
    const Faces faces(d);
    std::vector<FaceSet> possible(d.vertices.size(), FaceSet(faces.size(), true));
    Edges links;   // the edges between two undrawn vertices
    for (std::size_t k = i; k < edges.size(); ++k) {
        const std::size_t u = edges[k][0], v = edges[k][1];
        if (is_drawn(d, u) && is_drawn(d, v)) { // both present
            if (d.first_path(u, v).empty()) return false; // edge cannot be added
        } else if (is_drawn(d, u) || is_drawn(d, v)) { // one present
            const std::size_t a = (is_drawn(d, u) ? u : v), x = ((a == u) ? v : u);
            const FaceSet reach = ball(faces, faces_at(d, faces, a), klim, a);
            for (std::size_t f = 0; f < faces.size(); ++f) possible[x][f] = possible[x][f] && reach[f];
        } else // neither present
            links.push_back(edges[k]);
    }
    for (bool changed = true; changed;) {
        changed = false;
        for (const Edge& l : links)
            for (std::size_t s = 0; s < 2; ++s) {
                const FaceSet reach = ball(faces, possible[l[s]], klim);
                for (std::size_t f = 0; f < faces.size(); ++f)
                    if (possible[l[1 - s]][f] && !reach[f]) {
                        possible[l[1 - s]][f] = false;
                        changed = true;
                    }
            }
    }
    for (std::size_t v = 0; v < d.vertices.size(); ++v)
        if (!is_drawn(d, v) && std::find(possible[v].begin(), possible[v].end(), true) == possible[v].end())
            return false;
    return true;
}

// draw_edges of dodecahedron.cpp, with a switch for the pruning and a node counter
uint64_t g_nodes = 0;
double g_deadline = 1e30;
bool g_timeout = false;
void draw_edges(D& d, const Edges& edges, std::size_t i, std::size_t j, const std::function<void(D&)>& done,
                bool prune) {
    if (i == j) {
        done(d);
        return;
    }
    std::size_t a = edges[i][0];
    std::size_t b = edges[i][1];
    for (HdsPath p = d.first_path(a, b); !p.empty();) {
        d.add_edge(p, b);
        ++g_nodes;
        if ((g_nodes & 0xfff) == 0 && now_s() > g_deadline) g_timeout = true;
        if (!g_timeout && (!prune || may_complete(d, edges, i + 1))) draw_edges(d, edges, i + 1, j, done, prune);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, b)) break;
    }
}
// ---------------------------------------------------------------------------------------

// an instance: n vertices, the edges drawn before the search (first one with add_first_edge, the others with
// first_path), the edges to draw, and optionally look-ahead-only edges
struct Instance {
    std::size_t n;
    Edges fixed, draw, look;
};

Instance make(const std::string& name) {
    Instance I;
    auto kn_star = [&](std::size_t n, const std::set<Edge>& miss) {
        I.n = n;
        for (std::size_t v = 1; v < n; ++v) I.fixed.push_back({0, v});
        for (std::size_t v = 2; v < n; ++v)
            for (std::size_t u = 1; u < v; ++u)
                if (!miss.count({u, v})) I.draw.push_back({u, v});
    };
    if (name == "K7free") {   // K7, only the edge 01 fixed, vertices placed one by one
        I.n = 7;
        I.fixed.push_back({0, 1});
        for (std::size_t v = 2; v < 7; ++v)
            for (std::size_t u = 0; u < v; ++u) I.draw.push_back({u, v});
    } else if (name == "K8") kn_star(8, {});
    else if (name == "K9star3") kn_star(9, {{1, 2}, {1, 3}, {1, 4}});
    else if (name == "K9tri") kn_star(9, {{1, 2}, {1, 3}, {2, 3}});
    else if (name == "K9path3") kn_star(9, {{1, 2}, {2, 3}, {3, 4}});
    else if (name == "K44" || name == "K44v3") {   // A = 0..3, B = 4..7, star 0 -> B fixed
        I.n = 8;
        for (std::size_t j = 4; j < 8; ++j) I.fixed.push_back({0, j});
        for (std::size_t x = 1; x < 4; ++x)
            for (std::size_t j = 4; j < 8; ++j) (name == "K44v3" && x == 3 ? I.look : I.draw).push_back({j, x});
    } else if (name == "K8v7") {   // G' = K7 on 0..6 (star at 0 fixed), look-ahead: edges at 7
        I.n = 8;
        for (std::size_t v = 1; v < 7; ++v) I.fixed.push_back({0, v});
        for (std::size_t v = 2; v < 7; ++v)
            for (std::size_t u = 1; u < v; ++u) I.draw.push_back({u, v});
        for (std::size_t u = 0; u < 7; ++u) I.look.push_back({u, 7});
    } else if (name == "K7pend") {   // G' = K7 (star fixed), look-ahead: 7 ~ 0,1,2 and 8 ~ 3,7
        I.n = 9;
        for (std::size_t v = 1; v < 7; ++v) I.fixed.push_back({0, v});
        for (std::size_t v = 2; v < 7; ++v)
            for (std::size_t u = 1; u < v; ++u) I.draw.push_back({u, v});
        I.look = {{0, 7}, {1, 7}, {2, 7}, {3, 8}, {7, 8}};
    } else throw std::runtime_error("unknown instance " + name);
    return I;
}

D start_drawing(const Instance& I) {
    D d(I.n);
    d.add_first_edge(I.fixed[0][0], I.fixed[0][1]);
    for (std::size_t i = 1; i < I.fixed.size(); ++i) {
        HdsPath p = d.first_path(I.fixed[i][0], I.fixed[i][1]);
        if (p.empty()) throw std::runtime_error("fixed edges");
        d.add_edge(p, I.fixed[i][1]);
    }
    return d;
}

// enumerate the drawings of I.draw (with I.look as look-ahead edges if look) and return their labelled codes
std::set<Key> run(const Instance& I, bool prune, bool look, uint64_t& found, uint64_t& nodes, double& secs,
                  const std::function<void(D&)>& extra = nullptr) {
    Edges edges = I.draw;
    if (look) edges.insert(edges.end(), I.look.begin(), I.look.end());
    D d = start_drawing(I);
    std::set<Key> res;
    Canonizer<D> cz;
    found = 0;
    g_nodes = 0;
    const double t0 = now_s();
    draw_edges(d, edges, 0, I.draw.size(), [&](D& f) {
        ++found;
        res.insert(code_str(cz.labelled(f)));
        if (extra) extra(f);
    }, prune);
    secs = now_s() - t0;
    nodes = g_nodes;
    return res;
}

void part_a(const std::string& name) {
    Instance I = make(name);
    uint64_t f0, n0, f1, n1;
    double s0, s1;
    std::set<Key> A = run(I, false, false, f0, n0, s0);
    std::set<Key> B = run(I, true, false, f1, n1, s1);
    std::size_t common = 0;
    for (const Key& k : B) common += A.count(k);
    std::cout << "RESULT a " << name << " edges=" << I.fixed.size() + I.draw.size() << " | no pruning: found=" << f0
              << " distinct=" << A.size() << " nodes=" << n0 << " time=" << s0 << "s | pruning: found=" << f1
              << " distinct=" << B.size() << " nodes=" << n1 << " time=" << s1 << "s | common=" << common
              << (A == B ? " SAME" : " DIFFERENT") << std::endl;
}

void part_b(const std::string& name) {
    Instance I = make(name);
    // drawings of G' with the look-ahead edges (pruning on, as in step 2), and without look-ahead (pruning on/off)
    uint64_t fl, nl, fp, np_, fn, nn;
    double sl, sp, sn;
    std::set<Key> L = run(I, true, true, fl, nl, sl);
    std::set<Key> Pr = run(I, true, false, fp, np_, sp);
    std::set<Key> N = run(I, false, false, fn, nn, sn);
    // all drawings of G, without pruning, restricted to G' (the look-ahead edges are the last edges drawn)
    Instance G = I;
    G.draw.insert(G.draw.end(), I.look.begin(), I.look.end());
    G.look.clear();
    std::set<Key> restr;
    Canonizer<D> cz;
    uint64_t fg, ng;
    double sg;
    run(G, false, false, fg, ng, sg, [&](D& g) {
        D r(g);
        for (std::size_t k = 0; k < I.look.size(); ++k) r.remove_edge();
        restr.insert(code_str(cz.labelled(r)));
    });
    std::size_t missing = 0;
    for (const Key& k : restr) missing += !L.count(k);
    bool subset_n = true;
    for (const Key& k : L) subset_n = subset_n && N.count(k);
    std::cout << "RESULT b " << name << " | G' drawings: no pruning " << N.size() << " (nodes " << nn << ", " << sn
              << "s), pruning " << Pr.size() << " (nodes " << np_ << "), pruning+look-ahead " << L.size() << " (nodes "
              << nl << ", " << sl << "s)" << (subset_n ? " subset-of-all" : " NOT-SUBSET")
              << " | G drawings (no pruning): " << fg << " (nodes " << ng << ", " << sg << "s), distinct restrictions "
              << restr.size() << ", restrictions missing from look-ahead set: " << missing
              << (missing == 0 ? " OK" : " FAIL") << std::endl;
}

// ---- part c: the frame F of dodecahedron.cpp (dodecahedron + pentagon diagonals), copied edge order ----
std::vector<Vertices> pentagons() {
    std::vector<Vertices> res = {{0, 1, 2, 3, 4}, {15, 16, 17, 18, 19}};
    for (std::size_t i = 0; i < 5; ++i) {
        const std::size_t j = (i + 1) % 5;
        res.push_back({i, j, 5 + j, 10 + i, 5 + i});
        res.push_back({15 + i, 15 + j, 10 + j, 5 + j, 10 + i});
    }
    return res;
}
Edge sorted_edge(std::size_t a, std::size_t b) { return {std::min(a, b), std::max(a, b)}; }
Edges frame_edges() {
    std::set<Edge> de;
    for (const Vertices& p : pentagons())
        for (std::size_t i = 0; i < 5; ++i) de.insert(sorted_edge(p[i], p[(i + 1) % 5]));
    Edges res(de.begin(), de.end());
    for (const Vertices& p : pentagons())
        for (std::size_t i = 0; i < 5; ++i) res.push_back(sorted_edge(p[i], p[(i + 2) % 5]));
    std::stable_sort(res.begin(), res.end(), [](const Edge& e, const Edge& f) { return e[1] < f[1]; });
    return res;
}

void part_c(double seconds) {
    const Edges edges = frame_edges();
    // count the drawings of the first s edges with and without pruning, for growing s, within the time limit
    for (std::size_t s = 2; s <= 20; ++s) {
        for (int prune = 1; prune >= 0; --prune) {
            D d(20);
            d.add_first_edge(edges[0][0], edges[0][1]);
            uint64_t cnt = 0;
            g_nodes = 0;
            g_timeout = false;
            const double t0 = now_s();
            g_deadline = t0 + seconds;
            draw_edges(d, edges, 1, s, [&](D&) { ++cnt; }, prune);
            std::cout << "RESULT c first " << s << " edges of F, " << (prune ? "pruning   " : "no pruning")
                      << ": drawings=" << cnt << " nodes=" << g_nodes << " time=" << now_s() - t0 << "s"
                      << (g_timeout ? " (TIMEOUT, partial)" : "") << std::endl;
            if (g_timeout) return;
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " a|b <instance> | c <seconds>\n";
        return 2;
    }
    const std::string part = argv[1];
    if (part == "a") part_a(argv[2]);
    else if (part == "b") part_b(argv[2]);
    else if (part == "c") part_c(std::stod(argv[2]));
    std::cout << "peak_mem=" << peak_mem_mb() << "MB" << std::endl;
    return 0;
}

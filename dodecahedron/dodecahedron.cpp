#include "hds_kplanar.h"
#include "iso.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <thread>
#include <vector>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;
typedef std::vector<std::size_t> Vertices;
typedef std::vector<bool> FaceSet;   // FaceSet[f]: is the face f of a drawing in the set

const std::size_t klim = 3;
typedef Drawing<klim> D;
const std::string frame_dir = "../quasiDrawings/dodecahedron/frame/";         // output: the drawings of F
const std::string skeleton_dir = "../quasiDrawings/dodecahedron/skeleton/";   // output: the drawings of S
const std::size_t split = 13;      // step 1 first draws the first split edges of F, then continues in parallel threads
const std::size_t num_threads = 12;

//   Step 1: all drawings of the frame F, i.e. D and the diagonals (20 vertices, 90 edges). They are counted up to
//     isomorphism and saved in frame_dir.
//   Step 2: for every drawing of F, find all drawings of the skeleton S, i.e. T_1 without the vertices s_{p,i} (F with the
//     hubs h_p, the vertices t_{p,i} and the edges p_i t_{p,i} and t_{p,i} h_p), in which the vertices s_{p,i} can
//     still be added, s.t. every drawing of T_1 contains one of them. 
//     They are counted up to isomorphism, and saved in skeleton_dir.

// Both steps use draw_edges, which draws a list of edges one after the other in all possible ways. After each edge,
// may_complete checks that the edges after it can still be drawn, which cuts most branches early.
// Run from the repository root (about half an hour with 12 threads).

// labels
const std::size_t num_pentagons = 12;
const std::size_t n_frame = 20;      // vertices of F
const std::size_t n_skeleton = 92;   // vertices of S
const std::size_t n_t1 = 152;        // vertices of T_1
std::size_t hub(std::size_t p) { return 20 + p; } // hub h_p
std::size_t top(std::size_t p, std::size_t i) { return 32 + 5 * p + i; } // t_{p,i}
std::size_t pin(std::size_t p, std::size_t i) { return 92 + 5 * p + i; } // s_{p,i}

// the 12 pentagons of D, each with its vertices p_0, ..., p_4 in cyclic order
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

// the 30 edges of D
Edges dodecahedron_edges() {
    std::set<Edge> res;
    for (const Vertices& p : pentagons())
        for (std::size_t i = 0; i < 5; ++i) res.insert(sorted_edge(p[i], p[(i + 1) % 5]));
    return Edges(res.begin(), res.end());
}

// the 90 edges ab (a < b) of the frame F
Edges frame_edges() {
    Edges res = dodecahedron_edges();
    for (const Vertices& p : pentagons()) // diagonals
        for (std::size_t i = 0; i < 5; ++i) res.push_back(sorted_edge(p[i], p[(i + 2) % 5]));
    std::stable_sort(res.begin(), res.end(), [](const Edge& e, const Edge& f) { return e[1] < f[1]; });
    return res;
}

// the edges of S that are not in F: p_i t_{p,i} and t_{p,i} h_p for every pentagon p and every i
Edges skeleton_edges() {
    Edges res;
    const std::vector<Vertices> ps = pentagons();
    for (std::size_t p = 0; p < num_pentagons; ++p)
        for (std::size_t i = 0; i < 5; ++i) {
            res.push_back({ps[p][i], top(p, i)});
            res.push_back({top(p, i), hub(p)});
        }
    return res;
}

// the edges of T_1 that are not in S: p_i s_{p,i} and s_{p,i} t_{p,i} for every pentagon p and every i
Edges pinned_edges() {
    Edges res;
    const std::vector<Vertices> ps = pentagons();
    for (std::size_t p = 0; p < num_pentagons; ++p)
        for (std::size_t i = 0; i < 5; ++i) {
            res.push_back({ps[p][i], pin(p, i)});
            res.push_back({pin(p, i), top(p, i)});
        }
    return res;
}

bool is_drawn(const D& d, std::size_t v) { return d.vertices[v].halfedge != nullptr; }

// the faces of a drawing: face[h.label] is the face of the halfedge h, and doors[f] are the halfedges around the face f
// whose edge has a crossing left (an edge can leave f through them)
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

// the faces at the drawn vertex v
FaceSet faces_at(const D& d, const Faces& faces, std::size_t v) {
    FaceSet res(faces.size(), false);
    const HdsHalfedge* h = d.vertices[v].halfedge;   // points to v
    do { // walk around halfedges pointing to v
        res[faces.face[h->label]] = true;
        h = h->next->twin;
    } while (h != d.vertices[v].halfedge);
    return res;
}

// Heuristic to discard drawings early
// BFS to determine the faces reached from the faces in start by crossing at most steps edges 
// that have a crossing left and are not incident to the vertex avoid
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

// can the edges edges[i], edges[i+1], ... still be drawn? 
// Each of them between two drawn vertices needs a way to be drawn. 
// Each undrawn vertex x needs a face to be put in: one that the edge from every drawn neighbour a can reach
// (within klim crossings and not crossing the edges at a), and from which every undrawn neighbour y has a face of its own
// within klim crossings.
// The possible faces of the undrawn vertices are narrowed down until nothing changes.
// The answer may be yes for edges that cannot be drawn, but never no for edges that can.
bool may_complete(D& d, const Edges& edges, std::size_t i) {
    const Faces faces(d);
    std::vector<FaceSet> possible(d.vertices.size(), FaceSet(faces.size(), true));
    Edges links;   // the edges between two undrawn vertices
    for (std::size_t k = i; k < edges.size(); ++k) {
        const std::size_t u = edges[k][0], v = edges[k][1];
        if (is_drawn(d, u) && is_drawn(d, v)) { // both present
            if (d.first_path(u, v).empty()) return false; // edge cannot be added
        } else if (is_drawn(d, u) || is_drawn(d, v)) { // one present
            // a is the present vertex, x the non-present vertex
            const std::size_t a = (is_drawn(d, u) ? u : v), x = ((a == u) ? v : u);
            const FaceSet reach = ball(faces, faces_at(d, faces, a), klim, a);
            // mark faces reachable from a as possible for x (with &)
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
                    // remove possible faces of 1-s that s cannot reach
                    if (possible[l[1 - s]][f] && !reach[f]) {
                        possible[l[1 - s]][f] = false;
                        changed = true;
                    }
            }
    }
    for (std::size_t v = 0; v < d.vertices.size(); ++v)
        // for undrawn vertices, return false if they have no possible face
        if (!is_drawn(d, v) && std::find(possible[v].begin(), possible[v].end(), true) == possible[v].end())
            return false;
    return true;
}

// draws edges[i], ..., edges[j-1] (each from its first endpoint, which is drawn by then) in all possible ways, and
// calls done(d) for every drawing of them in which the edges edges[j], edges[j+1], ... may still be drawn
void draw_edges(D& d, const Edges& edges, std::size_t i, std::size_t j, const std::function<void(D&)>& done) {
    if (i == j) {
        done(d);
        return;
    }
    std::size_t a = edges[i][0];
    std::size_t b = edges[i][1];
    for (HdsPath p = d.first_path(a, b); !p.empty();) {
        d.add_edge(p, b);
        // only recurse if the may_complete() heuristic passes
        if (may_complete(d, edges, i + 1)) draw_edges(d, edges, i + 1, j, done);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, b)) break;
    }
}

// all drawings of F: the drawings of the first split edges, extended by num_threads threads; 
// found[i] are the drawings of F that contain partial[i], so every drawing of F is found exactly once
std::vector<nlohmann::json> frame_drawings() {
    const Edges edges = frame_edges();
    D d(n_frame);
    d.add_first_edge(edges[0][0], edges[0][1]);
    std::vector<D> partial; // partial edges on first split=13 edges
    draw_edges(d, edges, 1, split, [&](D& p) { partial.push_back(p); });
    // the copy constructor of Drawing redirects all pointers into the copy, so each copy is a drawing of its own
    std::cout << partial.size() << " drawings of the first " << split << " edges of F" << std::endl;

    // each thread repeatedly takes the next partial drawing that nobody has taken yet and extends it in place; 
    // the threads share nothing else, and each writes only to its own found[i]
    std::vector<std::vector<nlohmann::json>> found(partial.size());
    std::atomic<std::size_t> next(0); // atomic counter shared by all threads
    auto work = [&]() {
        for (std::size_t i = next++; i < partial.size(); i = next++)
            draw_edges(partial[i], edges, split, edges.size(), 
                    [&](D& f) { found[i].push_back(f.serialize_to_json()); });
    };
    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < num_threads; ++t) threads.emplace_back(work);
    for (std::thread& t : threads) t.join();

    std::vector<nlohmann::json> res;
    for (const auto& f : found) res.insert(res.end(), f.begin(), f.end());
    return res;
}

// the numbers of crossings of the diagonals, of the edges t_{p,i} h_p and of the other edges in the drawing s of S,
// e.g. "diagonals 3, edges t_{p,i} h_p 1, other edges 0" if all edges of a kind have the same number of crossings
std::string crossing_numbers(const D& s) {
    const Edges de = dodecahedron_edges();
    const std::set<Edge> dodecahedron(de.begin(), de.end());
    std::set<std::size_t> diagonals, spokes, others;   // the numbers of crossings that occur for each kind of edge
    for (const auto& e : s.edges) {
        const Edge uv = sorted_edge(e.u, e.v);
        if (uv[1] < n_frame && !dodecahedron.count(uv)) diagonals.insert(e.ncr);         // both ends in D, not an edge of D
        else if (uv[0] >= hub(0) && uv[0] < hub(num_pentagons)) spokes.insert(e.ncr);   // a hub and a t_{p,i}
        else others.insert(e.ncr);                                                      // edges of D and p_i t_{p,i}
    }
    auto list = [](const std::set<std::size_t>& x) {
        std::string res;
        for (std::size_t c : x) res += (res.empty() ? "" : " or ") + std::to_string(c);
        return res;
    };
    return "diagonals " + list(diagonals) + ", edges t_{p,i} h_p " + list(spokes) + ", other edges " + list(others);
}

void save(const D& d, const std::string& name) {
    std::ofstream out(name + ".json"), gml(name + ".graphml");
    out << d.serialize_to_json().dump(4);
    d.graphml_output(gml);
    out.close();
    gml.close();
    if (!out || !gml) throw std::runtime_error("could not write " + name + ".json and .graphml");
}

int main() {
    const auto start = std::chrono::steady_clock::now();
    std::filesystem::create_directories(frame_dir);
    std::filesystem::create_directories(skeleton_dir);

    // step 1
    const std::vector<nlohmann::json> frames = frame_drawings(); // all frame drawings
    std::vector<D> classes; // drawings of F up to isomorphism
    for (std::size_t i = 0; i < frames.size(); ++i) {
        D f(frames[i], n_frame);
        std::size_t c = 0;
        while (c < classes.size() && !are_isomorphic(classes[c], f)) ++c;
        if (c == classes.size()) classes.push_back(f);
        std::cout << "drawing " << i + 1 << " of F: " << f.crossings.size() << " crossings, class " << c + 1 << std::endl;
        save(f, frame_dir + "frame_" + std::to_string(i + 1));
    }
    std::cout << frames.size() << " drawings of F, " << classes.size() << " up to isomorphism" << std::endl;

    // step 2: the edges of S that are not in F, then the edges at the vertices s_{p,i}, which only cut branches
    Edges edges = skeleton_edges();
    const std::size_t num_skeleton_edges = edges.size();
    for (const Edge& e : pinned_edges()) edges.push_back(e); // for look-ahead only, pinned_edges not actually drawn
    std::vector<D> skeleton_classes;   // one drawing of S per isomorphism class
    std::size_t num_skeletons = 0;
    for (std::size_t i = 0; i < frames.size(); ++i) {
        D d(frames[i], n_t1);   // all vertices of T_1, only those of F drawn
        std::size_t found = 0;
        draw_edges(d, edges, 0, num_skeleton_edges, [&](D& d_s) {
            nlohmann::json j = d_s.serialize_to_json();
            j["num_vertices"] = n_skeleton;   // without the vertices s_{p,i}, which are not drawn
            const D s(j, n_skeleton);
            std::size_t c = 0;
            while (c < skeleton_classes.size() && !are_isomorphic(skeleton_classes[c], s)) ++c;
            if (c == skeleton_classes.size()) skeleton_classes.push_back(s);
            ++found;
            std::cout << "drawing " << found << " of S from drawing " << i + 1 << " of F: " << s.crossings.size()
                      << " crossings (" << crossing_numbers(s) << "), class " << c + 1 << std::endl;
            save(s, skeleton_dir + "skeleton_" + std::to_string(i + 1) + "_" + std::to_string(found));
        });
        std::cout << "drawing " << i + 1 << " of F: " << found << " drawings of S" << std::endl;
        num_skeletons += found;
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << num_skeletons << " drawings of S, " << skeleton_classes.size() << " up to isomorphism (" << secs
              << " s)" << std::endl;
    return num_skeletons > 0 ? 0 : 1;   // T_1 has a drawing, finding none is an error
}

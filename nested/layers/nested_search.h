// Search for the 3-plane drawings of one layer of the nested cycle graph: given a drawing whose last ("active")
// cycle is uncrossed or crossable, add the next layer (a new cycle and the matching to the active cycle) in all
// possible ways. Used by nested_layers.h. This is the part of iterative_nested_split.h (the exploration version,
// with options for extra edge classes etc.) that the runs of the thesis use.
//
//   Pass 1: the new cycle is crossable. Every complete drawing with a potential final face (a face that a further
//           uncrossed cycle can still be drawn in) is an intermediate drawing: it is reduced to its relevant part
//           (as in Hoffmann et al.) and reported with the canonical code of the reduced drawing.
//   Pass 2: the new cycle is uncrossed. Every complete drawing is reported.
//
// Pruning (all sound, i.e. only branches without a completion are cut): after every edge, every remaining edge
// must still be routable on its own (reachability check, which also reorders the remaining edges, fewest options
// first); at the checkpoints the remaining edges must be routable together (integral multicommodity flow
// relaxation). In Pass 1 additionally: after every edge, a potential final face must still be possible and every
// unplaced vertex must still have a face it can be put in; paths that put a new vertex into an impossible face are
// skipped.
//
// Splitting over tasks: the paths of the first split_depth edges that survive the pruning ("prefixes", "units") are
// numbered 0, 1, 2, ... in DFS order, identically in every task; a task only searches below the units with
// index % task_count == task_index. Finished units are listed in <output_dir>/done_units.txt and skipped when a
// killed task is restarted.

#ifndef NESTED_SEARCH_H
#define NESTED_SEARCH_H

#include "hds_kplanar.h"
#include <algorithm>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <boost/functional/hash.hpp>

namespace nested_search {

    typedef std::vector<std::size_t> Edge;  // {u, v} or {u, v, prescribed crossings in Pass 2}

    // Max flow for threshold tests "is the max flow from source to sink at least limit?": flow(source, sink,
    // limit) returns min(max flow, limit) using at most limit BFS augmenting paths (the limits are small).
    struct ThresholdFlow {
        std::vector<int> to, cap, orig;          // arc a: to[a], residual cap[a]; a^1 is its reverse
        std::vector<std::vector<int>> out;       // arcs leaving each node
        std::vector<int> parent;                 // BFS: arc used to reach a node
        std::vector<int> queue;

        explicit ThresholdFlow(std::size_t num_nodes) : out(num_nodes), parent(num_nodes) {}

        void add_edge(int from, int to_, int capacity) {
            out[from].push_back(static_cast<int>(to.size()));
            to.push_back(to_); orig.push_back(capacity);
            out[to_].push_back(static_cast<int>(to.size()));
            to.push_back(from); orig.push_back(0);
        }

        int flow(int source, int sink, int limit) {
            cap = orig;
            int f = 0;
            while (f < limit) {
                std::fill(parent.begin(), parent.end(), -1);
                parent[source] = -2;
                queue.assign(1, source);
                for (std::size_t qi = 0; qi < queue.size() && parent[sink] == -1; ++qi)
                    for (int a : out[queue[qi]])
                        if (cap[a] > 0 && parent[to[a]] == -1) {
                            parent[to[a]] = a;
                            queue.push_back(to[a]);
                        }
                if (parent[sink] == -1) break; // no augmenting path: f is the max flow
                int push = limit - f;
                for (int v = sink; v != source; v = to[parent[v] ^ 1]) push = std::min(push, cap[parent[v]]);
                for (int v = sink; v != source; v = to[parent[v] ^ 1]) {
                    cap[parent[v]] -= push;
                    cap[parent[v] ^ 1] += push;
                }
                f += push;
            }
            return f;
        }
    };

    // ----------------------------------------------------------------------------------------------------------
    // Canonical codes
    //
    // The next layer is attached to the active cycle vertices by label (0, ..., c-1 in a reduced drawing), so two
    // drawings are interchangeable only under isomorphisms that map the active cycle onto itself by a symmetry of
    // the layer edges. Mirroring is always allowed.

    typedef std::vector<std::size_t> CycleSym; // sym[i] = image of active cycle vertex i
    typedef std::vector<std::uint32_t> IsoCode;

    struct IsoCodeHash {
        std::size_t operator()(const IsoCode& c) const { return boost::hash_range(c.begin(), c.end()); }
    };

    // A kill (time limit) may cut off the last line of an append-only output file: remove such a line, so that the
    // next run of the task starts on a fresh line. Its record is written again, since its unit was not finished.
    inline void drop_partial_line(const std::string& path) {
        std::ifstream f(path, std::ios::ate | std::ios::binary);
        if (!f) return;
        const std::streamoff size = f.tellg();
        std::streamoff keep = size; // length up to and including the last '\n'
        while (keep > 0) {
            f.seekg(keep - 1);
            if (f.get() == '\n') break;
            --keep;
        }
        f.close();
        if (keep < size) std::filesystem::resize_file(path, static_cast<std::uintmax_t>(keep));
    }

    // All rotations/reflections sym of the old cycle 0, ..., c-1 such that relabelling the old endpoints by sym and
    // the new cycle by some rotation/reflection maps the layer edges (including the capacity entries) onto
    // themselves. The result is a group: all 24 dihedral maps for C12, the 7 even rotations and the 7
    // parity-swapping reflections for C14.
    inline std::vector<CycleSym> gadget_symmetries(
            std::size_t c, const std::function<std::vector<Edge>(std::size_t)>& builder) {
        const std::size_t nm = c; // offset of the new cycle, any offset >= c works
        const std::vector<Edge> edges = builder(nm);
        for (const Edge& e : edges)
            for (std::size_t i = 0; i < 2; ++i)
                if (e[i] >= nm + c) throw std::runtime_error("gadget_symmetries: layer edge out of range");

        typedef std::tuple<std::size_t, std::size_t, std::size_t> Key;
        auto edge_set = [&](const auto& map) {
            std::set<Key> s;
            for (const Edge& e : edges) {
                std::size_t a = map(e[0]), b = map(e[1]);
                s.emplace(std::min(a, b), std::max(a, b), e.size() > 2 ? e[2] + 1 : 0);
            }
            return s;
        };
        auto dihedral = [c](std::size_t i, std::size_t shift, bool refl) {
            return refl ? (shift + c - i) % c : (shift + i) % c;
        };
        const std::set<Key> original = edge_set([](std::size_t x) { return x; });

        std::vector<CycleSym> syms;
        for (int refl_old = 0; refl_old < 2; ++refl_old)
            for (std::size_t s_old = 0; s_old < c; ++s_old) {
                CycleSym sym(c);
                for (std::size_t i = 0; i < c; ++i) sym[i] = dihedral(i, s_old, refl_old);
                bool found = false;
                for (int refl_new = 0; refl_new < 2 && !found; ++refl_new)
                    for (std::size_t s_new = 0; s_new < c && !found; ++s_new) {
                        auto map = [&](std::size_t x) { return x < c ? sym[x] : nm + dihedral(x - nm, s_new, refl_new); };
                        found = (edge_set(map) == original);
                    }
                if (found) syms.push_back(std::move(sym));
            }
        return syms;
    }

    // Encode d by a BFS from halfedge h0 (pointing to the cycle vertex that sym maps to 0), with cycle vertex i
    // labelled sym[i] and all other vertices/crossings labelled c, c+1, ... in order of discovery; rotations are
    // reversed if mirror is set. Every entry records the neighbour's label, whether it is a crossing, and the
    // crossing count of the edge. The code is compared with best while it is built: stop as soon as it is larger,
    // replace best if it ends up smaller.
    // PRE: d is connected, apart from isolated vertices.
    template <int kplane>
    void encode_if_smaller(const Drawing<kplane>& d, const CycleSym& sym, bool mirror,
            const HdsHalfedge* h0, IsoCode& best, IsoCode& cur) {
        const std::uint32_t SEP = UINT32_MAX; // end of a rotation
        const std::size_t c = sym.size();
        const std::size_t nv = d.vertices.size();
        const std::size_t n = nv + d.crossings.size();
        const std::size_t unset = n + c;

        std::vector<std::size_t> phi(n, unset);
        for (std::size_t i = 0; i < c; ++i) phi[i] = sym[i];
        std::vector<char> queued(n, 0);
        std::size_t next_label = c, nqueued = 1;

        cur.clear();
        bool smaller = best.empty();
        auto emit = [&](std::size_t x) {
            if (!smaller) {
                if (x > best[cur.size()]) return false;
                if (x < best[cur.size()]) smaller = true;
            }
            cur.push_back(static_cast<std::uint32_t>(x));
            return true;
        };

        // isolated (reserved) vertices carry no structure, only their number is recorded
        std::size_t isolated = 0;
        for (const auto& v : d.vertices) if (!v.halfedge) ++isolated;

        if (!emit(nv) || !emit(d.crossings.size()) || !emit(d.edges.size()) || !emit(isolated)) return;
        // without fixed cycle labels (canonical_code_free) the start vertex gets the first label
        if (phi[h0->vertex->label] == unset) phi[h0->vertex->label] = next_label++;
        std::queue<const HdsHalfedge*> q;
        q.push(h0);
        queued[h0->vertex->label] = 1;
        while (!q.empty()) {
            const HdsHalfedge* s = q.front();
            q.pop();
            const HdsHalfedge* t = s;
            do {
                std::size_t w = t->twin->vertex->label;
                if (phi[w] == unset) phi[w] = next_label++;
                if (!queued[w]) { queued[w] = 1; ++nqueued; q.push(t->twin); }
                if (!emit(2 * phi[w] + (w >= nv ? 1 : 0)) || !emit(t->edge->ncr)) return;
                t = (mirror ? t->twin->prev : t->next->twin);
            } while (t != s);
            if (!emit(SEP)) return;
        }
        if (nqueued + isolated != n) throw std::runtime_error("canonical_code: disconnected drawing");
        if (smaller) best.swap(cur);
    }

    // Canonical code of d with respect to the gadget symmetries syms: the smallest encoding over all sym in syms,
    // both orientations and all start halfedges pointing to the cycle vertex that sym maps to 0. Two drawings get
    // the same code iff there is an isomorphism of their planarizations (possibly mirrored) that preserves
    // vertex/crossing type and crossing counts and acts on the active cycle as an element of syms.
    template <int kplane>
    IsoCode canonical_code(const Drawing<kplane>& d, const std::vector<CycleSym>& syms) {
        IsoCode best, cur;
        for (const CycleSym& sym : syms) {
            std::size_t u0 = std::find(sym.begin(), sym.end(), 0) - sym.begin();
            const HdsHalfedge* h = d.vertices[u0].halfedge;
            if (!h) throw std::runtime_error("canonical_code: isolated cycle vertex");
            for (int mirror = 0; mirror < 2; ++mirror) {
                const HdsHalfedge* s = h;
                do {
                    encode_if_smaller(d, sym, mirror, s, best, cur);
                    s = s->next->twin;
                } while (s != h);
            }
        }
        return best;
    }

    // Canonical code for strong isomorphism without fixed labels: equal codes iff the planarizations are isomorphic
    // (possibly mirrored), preserving vertex/crossing type and crossing counts. Used for the complete drawings.
    template <int kplane>
    IsoCode canonical_code_free(const Drawing<kplane>& d) {
        IsoCode best, cur;
        const CycleSym none;
        for (const auto& v : d.vertices) {
            if (!v.halfedge) continue;
            for (int mirror = 0; mirror < 2; ++mirror) {
                const HdsHalfedge* s = v.halfedge;
                do {
                    encode_if_smaller(d, none, mirror, s, best, cur);
                    s = s->next->twin;
                } while (s != v.halfedge);
            }
        }
        return best;
    }

    // For the reduced drawing: a halfedge pointing to the source vertex of the relevant edge piece containing e,
    // and the number of crossings this piece gets in the reduced drawing (vact[v] == vact.size(): v irrelevant)
    template <std::size_t klim, typename F, typename A, typename V>
    inline std::pair<const HdsHalfedge*, std::size_t> prev_active(const HdsHalfedge* e, const F& face, const A& active, const V& vact) {
        if (active[face[e->label]] >= 0) {
            std::size_t c = (active[face[e->twin->label]] < 0 ? klim : e->edge->ncr);
            do e = e->prev; while (vact[e->vertex->label] == vact.size());
            return std::make_pair(e, c);
        }
        e = e->twin;
        if (active[face[e->label]] < 0) throw std::runtime_error("no active face");
        while (vact[e->vertex->label] == vact.size()) e = e->next;
        return std::make_pair(e, (active[face[e->twin->label]] < 0 ? klim : e->edge->ncr));
    }

    // ----------------------------------------------------------------------------------------------------------

    template <std::size_t klim = 3>
    class LayerSearcher {
    public:
        struct Config {
            std::size_t cycle_size = 12;
            // the layer edges for a new cycle starting at label nm: the new cycle (with klim as third entry, i.e.
            // uncrossable in Pass 2) and the matching from the active cycle 0, ..., cycle_size-1
            std::function<std::vector<Edge>(std::size_t nm)> layer_edges;
            std::vector<std::size_t> early_prune_checkpoints; // multicommodity flow check after these edges
            bool pass1 = true;                                // false: Pass 2

            // the drawing to extend: nullptr for the base drawing (D1 with an uncrossable star inside), or a reduced
            // intermediate drawing; in both, the active cycle is 0, ..., cycle_size-1 and there are no crossings
            const Drawing<klim>* input_drawing = nullptr;

            // splitting over tasks (split_depth = 0: no split, search everything, no files)
            std::size_t split_depth = 0, task_index = 0, task_count = 1;
            std::string output_dir;                       // for done_units.txt

            // Pass 1: every intermediate drawing, reduced (new cycle 0, ..., cycle_size-1), with its canonical code
            // (called for every one, also for isomorphic ones)
            std::function<void(const Drawing<klim>&, const IsoCode&)> intermediate_callback;
            // Pass 2: every complete drawing
            std::function<void(const Drawing<klim>&)> solution_callback;
        };

        explicit LayerSearcher(Config config) : config_(std::move(config)) {}

        // runs the search; returns the number of units (prefixes) at split_depth, 1 without a split
        std::size_t run();

        static std::size_t label_faces(const Drawing<klim>& d, std::vector<int>& face) {
            face.assign(d.halfedges.size(), -1);
            std::size_t num_faces = 0;
            for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                if (face[i->label] != -1) continue;
                const HdsHalfedge* j = &*i;
                do { face[j->label] = static_cast<int>(num_faces); j = j->next; } while (j->label != i->label);
                ++num_faces;
            }
            return num_faces;
        }

    private:
        static constexpr std::size_t mcf_node_budget = 200000; // backtracking nodes before giving up (no prune)
        static constexpr std::size_t mcf_max_paths = 4096;     // commodities with more candidate paths are dropped
        enum class McfResult { Feasible, Infeasible, Unknown };

        Drawing<klim> create_base_drawing() const;
        static std::vector<int> incident_faces(const Drawing<klim>& d, const std::vector<int>& face, std::size_t v);
        McfResult check_remaining_edges_mcf(const Drawing<klim>& d, const std::vector<Edge>& edges, std::size_t next) const;
        bool check_remaining_edges_reachability(const Drawing<klim>& d, std::vector<Edge>& edges, std::size_t next) const;
        bool placements_still_possible(const Drawing<klim>& d, const std::vector<Edge>& edges, std::size_t next,
                std::vector<int>* face_out = nullptr, std::size_t* num_faces_out = nullptr,
                std::map<std::size_t, std::vector<char>>* domains_out = nullptr) const;
        bool final_face_still_possible(const Drawing<klim>& d, std::size_t nm) const;
        void intermediate(const Drawing<klim>& d, std::size_t nm, const std::vector<CycleSym>& syms) const;

        Config config_;
    };

    // D1 (uncrossable) with an uncrossable star inside, at vertex cycle_size
    template <std::size_t klim>
    Drawing<klim> LayerSearcher<klim>::create_base_drawing() const {
        const std::size_t c = config_.cycle_size;
        Drawing<klim> d(c + 1);
        std::vector<HdsHalfedge*> cycle(c, nullptr);
        cycle[1] = d.add_first_edge(0, 1, klim);
        for (std::size_t i = 2; i < c; ++i) cycle[i] = d.add_edge(HdsPath({cycle[i - 1], nullptr}), i, klim);
        cycle[0] = d.add_edge(HdsPath({cycle[c - 1], cycle[1]->twin}), 0, klim);
        auto e = d.add_edge(HdsPath({cycle[0], nullptr}), c, klim);
        for (std::size_t i = c - 1; i > 0; --i) e = d.add_edge(HdsPath({cycle[i], e}), c, klim);
        return d;
    }

    // all faces incident to vertex v (empty if v is isolated)
    template <std::size_t klim>
    std::vector<int> LayerSearcher<klim>::incident_faces(const Drawing<klim>& d, const std::vector<int>& face, std::size_t v) {
        std::vector<int> faces;
        if (v >= d.vertices.size() || !d.vertices[v].halfedge) return faces;
        const HdsHalfedge* start = d.vertices[v].halfedge;
        const HdsHalfedge* h = start;
        do {
            const int f = face[h->label];
            if (f >= 0 && std::find(faces.begin(), faces.end(), f) == faces.end()) faces.push_back(f);
            h = h->next->twin;
        } while (h != start);
        return faces;
    }

    // Integral multicommodity flow relaxation for the remaining edges: every remaining edge (commodity) picks one
    // sequence of crossed drawn edges, such that no drawn edge exceeds its remaining capacity klim - ncr. Crossings
    // among the remaining edges are ignored, so Infeasible proves that no completion exists. Unknown: the node
    // budget was hit.
    template <std::size_t klim>
    typename LayerSearcher<klim>::McfResult LayerSearcher<klim>::check_remaining_edges_mcf(
            const Drawing<klim>& d, const std::vector<Edge>& edges, std::size_t next) const {
        std::vector<int> face;
        const std::size_t num_faces = label_faces(d, face);
        if (num_faces == 0) return McfResult::Unknown;

        std::vector<int> cap(d.edges.size(), 0); // remaining crossing capacity per drawn edge
        for (const auto& e : d.edges) cap[e.label] = static_cast<int>(klim) - static_cast<int>(e.ncr);

        // dual graph: face -> (adjacent face, crossed edge)
        std::vector<std::vector<std::pair<int, const HdsEdge*>>> dual(num_faces);
        for (const auto& h : d.halfedges) {
            if (h.edge->ncr >= klim) continue;
            const int f1 = face[h.label], f2 = face[h.twin->label];
            if (f1 != f2) dual[f1].push_back({f2, h.edge});
        }

        typedef std::vector<std::size_t> McfPath; // sorted labels of crossed edges
        std::vector<std::vector<McfPath>> paths;  // candidate paths per commodity
        std::vector<char> is_target(num_faces), in_path(num_faces);
        McfPath cur;

        for (std::size_t i = next; i < edges.size(); ++i) {
            const std::size_t u = edges[i][0], v = edges[i][1];
            const std::vector<int> u_faces = incident_faces(d, face, u), v_faces = incident_faces(d, face, v);
            if (u_faces.empty() || v_faces.empty()) continue; // unplaced endpoint: can go anywhere, drop

            std::fill(is_target.begin(), is_target.end(), 0);
            for (int f : v_faces) is_target[f] = 1;
            bool free_edge = false;
            for (int f : u_faces) if (is_target[f]) free_edge = true;
            if (free_edge) continue; // can be drawn without crossings, uses no capacity

            std::fill(in_path.begin(), in_path.end(), 0);
            for (int f : u_faces) in_path[f] = 1; // paths through another face of u are dominated
            std::vector<McfPath> cp;
            bool overflow = false;
            auto dfs = [&](auto& self, int f, std::size_t depth) -> void {
                if (depth >= klim) return;
                for (const auto& [nb, e] : dual[f]) {
                    if (overflow) return;
                    if (in_path[nb]) continue;
                    // as in Drawing::find_crossing: no crossing of edges sharing an endpoint, no edge crossed twice
                    if (e->u == u || e->v == u || e->u == v || e->v == v) continue;
                    if (std::find(cur.begin(), cur.end(), e->label) != cur.end()) continue;
                    cur.push_back(e->label);
                    if (is_target[nb]) { // any extension would be a dominated superset
                        cp.push_back(cur);
                        std::sort(cp.back().begin(), cp.back().end());
                        if (cp.size() > mcf_max_paths) overflow = true;
                    } else {
                        in_path[nb] = 1;
                        self(self, nb, depth + 1);
                        in_path[nb] = 0;
                    }
                    cur.pop_back();
                }
            };
            for (int f : u_faces) {
                cur.clear();
                dfs(dfs, f, 0);
                if (overflow) break;
            }
            if (overflow) continue; // too many options, drop commodity
            if (cp.empty()) return McfResult::Infeasible;

            // dedupe and remove dominated (superset) paths
            std::sort(cp.begin(), cp.end(), [](const McfPath& a, const McfPath& b) {
                return a.size() != b.size() ? a.size() < b.size() : a < b;
            });
            cp.erase(std::unique(cp.begin(), cp.end()), cp.end());
            std::vector<McfPath> kept;
            for (const auto& p : cp) {
                bool dominated = false;
                for (const auto& q : kept) {
                    if (q.size() >= p.size()) break; // kept is ordered by size
                    if (std::includes(p.begin(), p.end(), q.begin(), q.end())) { dominated = true; break; }
                }
                if (!dominated) kept.push_back(p);
            }
            paths.push_back(std::move(kept));
        }

        // backtracking, always on the commodity with the fewest fitting paths
        const std::size_t nc = paths.size();
        std::vector<char> assigned(nc, 0);
        std::size_t nodes = 0;
        bool aborted = false;
        auto fits = [&](const McfPath& p) {
            for (std::size_t l : p) if (cap[l] <= 0) return false;
            return true;
        };
        auto solve = [&](auto& self, std::size_t remaining) -> bool {
            if (remaining == 0) return true;
            if (++nodes > mcf_node_budget) { aborted = true; return false; }
            std::size_t best = nc, best_cnt = SIZE_MAX;
            for (std::size_t i = 0; i < nc; ++i) {
                if (assigned[i]) continue;
                std::size_t cnt = 0;
                for (const auto& p : paths[i]) {
                    if (fits(p)) ++cnt;
                    if (cnt >= best_cnt) break;
                }
                if (cnt == 0) return false;
                if (cnt < best_cnt) { best_cnt = cnt; best = i; }
            }
            assigned[best] = 1;
            for (const auto& p : paths[best]) {
                if (!fits(p)) continue;
                for (std::size_t l : p) --cap[l];
                if (self(self, remaining - 1)) return true;
                for (std::size_t l : p) ++cap[l];
                if (aborted) return false;
            }
            assigned[best] = 0;
            return false;
        };
        if (solve(solve, nc)) return McfResult::Feasible;
        return aborted ? McfResult::Unknown : McfResult::Infeasible;
    }

    // Every remaining edge with both endpoints placed must have a path with at most klim crossings (on its own).
    // Also reorders the remaining edges: fewest paths (1, 2, 3) first, the others keep their order. The result
    // depends on the current order, i.e. on the search history (see the split snapshot in run()).
    template <std::size_t klim>
    bool LayerSearcher<klim>::check_remaining_edges_reachability(const Drawing<klim>& d, std::vector<Edge>& edges, std::size_t next) const {
        std::vector<int> face;
        const std::size_t num_faces = label_faces(d, face);
        if (num_faces == 0) return false;

        // dual graph over crossable edges: face -> (adjacent face, crossed edge)
        std::vector<std::vector<std::pair<int, const HdsEdge*>>> dual(num_faces);
        for (const auto& h : d.halfedges)
            if (h.edge->ncr < klim && face[h.label] != face[h.twin->label])
                dual[face[h.label]].push_back({face[h.twin->label], h.edge});

        std::vector<std::size_t> crossed; // edges crossed on the current DFS path
        std::vector<std::pair<std::size_t, Edge>> counted; // (number of paths, capped at 4; edge)
        std::vector<bool> is_target(num_faces), in_path(num_faces);
        for (std::size_t i = next; i < edges.size(); ++i) {
            const std::size_t u = edges[i][0], v = edges[i][1];
            const std::vector<int> u_faces = incident_faces(d, face, u), v_faces = incident_faces(d, face, v);
            // an edge with an unplaced endpoint can always be drawn (the endpoint goes into a suitable face)
            if (u_faces.empty() || v_faces.empty()) { counted.push_back({4, edges[i]}); continue; }

            std::fill(is_target.begin(), is_target.end(), false);
            for (int f : v_faces) is_target[f] = true;
            std::fill(in_path.begin(), in_path.end(), false);
            std::size_t path_count = 0;
            auto dfs = [&](auto& self, int f, std::size_t depth) -> void {
                if (path_count > 3) return;
                if (is_target[f]) { ++path_count; return; }
                if (depth >= klim) return;
                in_path[f] = true;
                for (const auto& [nb, e] : dual[f]) {
                    if (in_path[nb]) continue;
                    // as in Drawing::find_crossing: no crossing of edges sharing an endpoint, no edge crossed twice
                    if (e->u == u || e->v == u || e->u == v || e->v == v) continue;
                    if (std::find(crossed.begin(), crossed.end(), e->label) != crossed.end()) continue;
                    crossed.push_back(e->label);
                    self(self, nb, depth + 1);
                    crossed.pop_back();
                    if (path_count > 3) break;
                }
                in_path[f] = false;
            };
            for (int f : u_faces) {
                dfs(dfs, f, 0);
                if (path_count > 3) break;
            }
            if (path_count == 0) return false;
            counted.push_back({path_count, edges[i]});
        }
        std::stable_sort(counted.begin(), counted.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        for (std::size_t i = 0; i < counted.size(); ++i) edges[next + i] = counted[i].second;
        return true;
    }

    // Placement check for the remaining edges (constraint propagation). Every unplaced vertex x will lie in some
    // current face. A remaining edge with at most L crossings between x and a placed vertex p forces x into a face
    // within dual distance L of p (over crossable edges not incident to p); an edge between two unplaced vertices
    // forces their faces within distance L of each other. The possible faces of the unplaced vertices are narrowed
    // until nothing changes; false if some vertex has none left. Optional outputs (for the placement filter): the
    // face labels of d, the number of faces and the possible faces of the unplaced vertices.
    template <std::size_t klim>
    bool LayerSearcher<klim>::placements_still_possible(const Drawing<klim>& d, const std::vector<Edge>& edges,
            std::size_t next, std::vector<int>* face_out, std::size_t* num_faces_out,
            std::map<std::size_t, std::vector<char>>* domains_out) const {
        std::vector<int> face;
        const std::size_t num_faces = label_faces(d, face);
        if (face_out) *face_out = face;
        if (num_faces_out) *num_faces_out = num_faces;
        std::vector<std::vector<std::pair<int, const HdsEdge*>>> dual(num_faces);
        for (const auto& h : d.halfedges)
            if (h.edge->ncr < klim && face[h.label] != face[h.twin->label])
                dual[face[h.label]].push_back({face[h.twin->label], h.edge});

        // faces within distance limit of the start set, not crossing edges incident to avoid (SIZE_MAX: any)
        std::vector<int> dist(num_faces);
        std::vector<int> queue;
        auto ball = [&](const std::vector<char>& start, std::size_t limit, std::size_t avoid) {
            std::fill(dist.begin(), dist.end(), -1);
            queue.clear();
            for (std::size_t f = 0; f < num_faces; ++f) if (start[f]) { dist[f] = 0; queue.push_back(static_cast<int>(f)); }
            for (std::size_t qi = 0; qi < queue.size(); ++qi) {
                const int f = queue[qi];
                if (static_cast<std::size_t>(dist[f]) >= limit) continue;
                for (const auto& [g, e] : dual[f]) {
                    if (dist[g] >= 0) continue;
                    if (avoid != SIZE_MAX && (e->u == avoid || e->v == avoid)) continue;
                    dist[g] = dist[f] + 1;
                    queue.push_back(g);
                }
            }
            std::vector<char> r(num_faces);
            for (std::size_t f = 0; f < num_faces; ++f) r[f] = dist[f] >= 0;
            return r;
        };

        std::map<std::size_t, std::vector<char>> domain; // unplaced vertex -> possible faces
        auto placed = [&](std::size_t x) { return d.vertices[x].halfedge != nullptr; };
        struct Link { std::size_t x, y, limit; };
        std::vector<Link> links; // remaining edges between two unplaced vertices
        for (std::size_t i = next; i < edges.size(); ++i) {
            const Edge& e = edges[i];
            const std::size_t pcr = (!config_.pass1 && e.size() == 3) ? e[2] : 0;
            const std::size_t limit = klim >= pcr ? klim - pcr : 0;
            const bool pu = placed(e[0]), pv = placed(e[1]);
            if (pu && pv) continue; // handled by the reachability check
            for (std::size_t x : {e[0], e[1]})
                if (!placed(x) && !domain.count(x)) domain[x].assign(num_faces, 1);
            if (!pu && !pv) { links.push_back({e[0], e[1], limit}); continue; }
            const std::size_t p = pu ? e[0] : e[1], x = pu ? e[1] : e[0];
            std::vector<char> start(num_faces, 0);
            for (int f : incident_faces(d, face, p)) start[f] = 1;
            const std::vector<char> r = ball(start, limit, p);
            auto& dx = domain[x];
            bool any = false;
            for (std::size_t f = 0; f < num_faces; ++f) { dx[f] = dx[f] && r[f]; any = any || dx[f]; }
            if (!any) return false;
        }
        for (bool changed = true; changed;) {
            changed = false;
            for (const Link& l : links)
                for (int side = 0; side < 2; ++side) {
                    auto& from = domain[side ? l.y : l.x];
                    auto& to = domain[side ? l.x : l.y];
                    const std::vector<char> r = ball(from, l.limit, SIZE_MAX);
                    bool any = false;
                    for (std::size_t f = 0; f < num_faces; ++f) {
                        if (to[f] && !r[f]) { to[f] = 0; changed = true; }
                        any = any || to[f];
                    }
                    if (!any) return false;
                }
        }
        if (domains_out) *domains_out = std::move(domain);
        return true;
    }

    // Pass 1, once all new cycle vertices nm, ..., nm+c-1 are placed: can some face of a completion still be a
    // potential final face (flow >= c from the new cycle vertices)? Faces of a completion refine the current faces,
    // so contracting them maps every flow of the leaf network onto this network; an edge piece with x crossings
    // can be split into j pieces later, which together carry at most bound(x) = max_j j*(klim-x-j+1) per direction.
    // False only if no completion has a potential final face, so pruning on false is sound.
    template <std::size_t klim>
    bool LayerSearcher<klim>::final_face_still_possible(const Drawing<klim>& d, std::size_t nm) const {
        const std::size_t c = config_.cycle_size;
        for (std::size_t i = 0; i < c; ++i)
            if (!d.vertices[nm + i].halfedge) return true; // not all placed: no test
        std::vector<int> face;
        const std::size_t num_faces = label_faces(d, face);
        auto bound = [](std::size_t x) {
            int best = 0;
            for (int j = 1; j <= static_cast<int>(klim) - static_cast<int>(x) + 1; ++j)
                best = std::max(best, j * (static_cast<int>(klim) - static_cast<int>(x) - j + 1));
            return best;
        };
        const std::size_t source = num_faces + c;
        ThresholdFlow net(source + 1);
        std::vector<int> in_cap(num_faces, 0); // cheap upper bound on the flow into each face
        for (std::size_t i = 0; i < c; ++i) net.add_edge(static_cast<int>(source), static_cast<int>(num_faces + i), 1);
        for (const auto& h : d.halfedges) {
            const std::size_t w = h.vertex->label;
            if (w >= nm && w < nm + c) {
                net.add_edge(static_cast<int>(num_faces + (w - nm)), face[h.label], 1);
                ++in_cap[face[h.label]];
            }
            const int b = bound(h.edge->ncr);
            if (b > 0 && face[h.label] != face[h.twin->label]) {
                net.add_edge(face[h.label], face[h.twin->label], b);
                in_cap[face[h.twin->label]] += b;
            }
        }
        for (std::size_t f = 0; f < num_faces; ++f)
            if (in_cap[f] >= static_cast<int>(c) && net.flow(static_cast<int>(source), static_cast<int>(f), static_cast<int>(c)) >= static_cast<int>(c))
                return true;
        return false;
    }

    // A complete Pass 1 drawing d (new cycle nm, ..., nm+c-1): if it has a potential final face, reduce it to its
    // relevant part (faces that a further layer can use: possible final faces, active, passive and transit faces,
    // as in Hoffmann et al., with the distances adapted to klim = 3), with an uncrossable star in every large
    // irrelevant region, and report it with the canonical code of the reduced drawing.
    template <std::size_t klim>
    void LayerSearcher<klim>::intermediate(const Drawing<klim>& d, std::size_t nm, const std::vector<CycleSym>& syms) const {
        const std::size_t c = config_.cycle_size;
        std::vector<int> face;
        const std::size_t num_faces = label_faces(d, face);
        if (num_faces == 0) return;

        // flow network: faces 0, ..., num_faces-1, new cycle vertices num_faces, ..., num_faces+c-1, source
        const std::size_t source = num_faces + c;
        ThresholdFlow net(source + 1);
        std::vector<std::vector<std::size_t>> dualg(num_faces);
        for (std::size_t i = 0; i < c; ++i) net.add_edge(static_cast<int>(source), static_cast<int>(num_faces + i), 1);
        const std::size_t cycle_end = nm + c;
        for (const auto& h : d.halfedges) {
            if (h.vertex->label >= nm && h.vertex->label < cycle_end)
                net.add_edge(static_cast<int>(num_faces + (h.vertex->label - nm)), face[h.label], 1);
            const int remaining = static_cast<int>(klim) - static_cast<int>(h.edge->ncr);
            if (remaining > 0) {
                net.add_edge(face[h.label], face[h.twin->label], remaining);
                dualg[face[h.label]].push_back(face[h.twin->label]);
            }
        }
        auto flow_at_least = [&](std::size_t s, std::size_t t, int limit) {
            return net.flow(static_cast<int>(s), static_cast<int>(t), limit) >= limit;
        };

        // potential final faces: flow >= c from the new cycle vertices
        std::vector<int> pff;
        for (std::size_t f = 0; f < num_faces; ++f)
            if (flow_at_least(source, f, static_cast<int>(c))) pff.push_back(static_cast<int>(f));
        if (pff.empty()) return;

        // face classes: 3 = potential final face, 2 = active, 1 = passive, 0 = transit, -1 = irrelevant
        std::vector<int> active(num_faces, -1);
        for (int f : pff) active[f] = 3;
        // active: flow >= 3 to a potential final face
        for (std::size_t i = 0; i < num_faces; ++i) {
            if (active[i] != -1) continue;
            for (int f : pff)
                if (flow_at_least(i, f, 3)) { active[i] = 2; break; }
        }
        std::vector<const HdsHalfedge*> fhedge(num_faces, nullptr); // a boundary halfedge of every face
        for (const auto& h : d.halfedges) if (!fhedge[face[h.label]]) fhedge[face[h.label]] = &h;

        auto dual_distance = [&](auto in_set) { // dual distances (over crossable edges) from a set of faces
            std::vector<int> dist(num_faces, -1);
            std::vector<std::size_t> queue;
            for (std::size_t f = 0; f < num_faces; ++f) if (in_set(f)) { dist[f] = 0; queue.push_back(f); }
            for (std::size_t qi = 0; qi < queue.size(); ++qi)
                for (std::size_t g : dualg[queue[qi]])
                    if (dist[g] < 0) { dist[g] = dist[queue[qi]] + 1; queue.push_back(g); }
            return dist;
        };
        const int kmax = static_cast<int>(klim);
        const std::vector<int> dist_active = dual_distance([&](std::size_t f) { return active[f] >= 2; });

        // passive: incident to a new cycle vertex v, with a dual path of length <= klim to an active face whose
        // first edge is not incident to v (an edge from v to the next cycle has at most klim crossings)
        for (std::size_t i = 0; i < num_faces; ++i) {
            if (active[i] != -1) continue;
            const HdsHalfedge* h = fhedge[i];
            do {
                const std::size_t v = h->vertex->label;
                if (v >= nm && v < cycle_end)
                    for (const HdsHalfedge* f = h->next->next; f != h; f = f->next)
                        if (f->edge->ncr < klim && f->vertex->label != v && f->twin->vertex->label != v) {
                            const int dfn = dist_active[face[f->twin->label]];
                            if (dfn >= 0 && dfn <= kmax - 1) { active[i] = 1; break; }
                        }
                h = h->next;
            } while (active[i] == -1 && h != fhedge[i]);
        }

        // transit: on a dual path of length <= klim from an active/passive face to an active face
        const std::vector<int> dist_ap = dual_distance([&](std::size_t f) { return active[f] >= 1; });
        for (std::size_t i = 0; i < num_faces; ++i)
            if (active[i] == -1 && dist_ap[i] >= 1 && dist_active[i] >= 1 && dist_ap[i] + dist_active[i] <= kmax)
                active[i] = 0;

        // relevant vertices and crossings (incident to >= 2 relevant faces, or new cycle vertices, which get the
        // labels 0, ..., c-1); possibly relevant ones (one relevant face) are kept only to avoid lenses
        std::vector<const HdsHalfedge*> relv(c, nullptr), prelv;
        auto classify = [&](const HdsHalfedge* start, std::size_t label) -> bool {
            const HdsHalfedge* e = start;
            const HdsHalfedge* a = e;
            std::size_t nirf = 0; // incident relevant faces
            do {
                if (active[face[e->label]] >= 0) { ++nirf; a = e; }
                e = e->next->twin;
            } while (e != start);
            if (label >= nm && label < cycle_end) {
                if (nirf == 0) return false;
                relv[label - nm] = a;
            } else if (nirf >= 2) relv.push_back(a);
            else if (nirf == 1) prelv.push_back(a);
            return true;
        };
        for (const auto& v : d.vertices)
            if (v.halfedge && !classify(v.halfedge, v.label)) return;
        for (const auto& x : d.crossings) classify(x.halfedge, SIZE_MAX);

        std::vector<std::size_t> fsize(num_faces, 0); // number of relevant vertices on each face
        const std::size_t total = d.vertices.size() + d.crossings.size();
        std::vector<std::size_t> vact(total, total); // old label -> label in the reduced drawing (total: irrelevant)
        std::size_t vic = 0;
        for (const HdsHalfedge* r : relv) {
            const HdsHalfedge* j = r;
            do { ++fsize[face[j->label]]; j = j->next->twin; } while (j != r);
            vact[r->vertex->label] = vic++;
        }
        for (const HdsHalfedge* r : prelv)
            if (fsize[face[r->label]] < 3) { // would become a lens: keep the vertex
                ++fsize[face[r->label]];
                relv.push_back(r);
                vact[r->vertex->label] = vic++;
            }

        // build the reduced drawing nd, merging connected irrelevant regions into "black" faces
        Drawing<klim> nd(relv.size() + 1);
        std::deque<const HdsHalfedge*> todo(1, relv[0]);
        std::vector<int> done(total, -1); // -1: not seen, 0: queued and in nd, 1: incident edges in nd
        std::vector<HdsHalfedge*> black;
        while (!todo.empty()) {
            const HdsHalfedge* i = todo.front();
            todo.pop_front();
            const std::size_t u = i->vertex->label;
            if (active[face[i->label]] < 0) throw std::runtime_error("irrelevant face");
            HdsHalfedge* last = nd.vertices[vact[u]].halfedge;
            if (nd.edges.empty()) {
                const auto p = prev_active<klim>(i, face, active, vact);
                const std::size_t v = p.first->vertex->label;
                last = nd.add_first_edge(vact[u], vact[v], p.second)->twin;
                if (last->vertex->label != 0) throw std::runtime_error("wrong first edge in reduced drawing");
                todo.push_back(p.first);
                done[v] = 0;
            } else if (last == nullptr) {
                todo.push_back(i);
                continue;
            }
            // u is connected in nd and last points to it: add its remaining edges in rotation order
            done[u] = 1;
            const std::size_t wnd = last->twin->vertex->label; // neighbour of u in nd
            std::size_t w = 0;
            const HdsHalfedge* wi = i;
            do {
                w = prev_active<klim>(wi, face, active, vact).first->vertex->label;
                if (vact[w] == wnd) break;
                wi = wi->next->twin;
            } while (wi != i);
            if (vact[w] != wnd) throw std::runtime_error("no w");
            for (const HdsHalfedge* j = wi->next->twin; j != wi; j = j->next->twin) {
                if (active[face[j->label]] < 0 && active[face[j->twin->label]] < 0) continue;
                const auto p = prev_active<klim>(j, face, active, vact);
                const std::size_t v = p.first->vertex->label;
                if (last->next->vertex->label == vact[v]) { last = last->next->twin; continue; } // already there
                HdsPath path(2, last);
                if (!nd.find_target(path, vact[v])) throw std::runtime_error("cannot find v");
                HdsHalfedge* ne = nd.add_edge(path, vact[v], p.second);
                if (active[face[j->label]] < 0) black.push_back(ne->twin);
                else if (active[face[j->twin->label]] < 0) black.push_back(ne);
                last = last->next->twin;
                if (done[v] < 0) { done[v] = 0; todo.push_back(p.first); }
            }
        }

        // an uncrossable star in every black face with at least 4 sides
        std::vector<bool> blackdone(nd.halfedges.size(), false);
        std::vector<HdsHalfedge*> large_black_faces;
        for (HdsHalfedge* b : black) {
            if (blackdone[b->label]) continue;
            HdsHalfedge* j = b;
            std::size_t sides = 0;
            do { blackdone[j->label] = true; ++sides; j = j->next; } while (j != b);
            if (sides >= 4) large_black_faces.push_back(b);
        }
        if (large_black_faces.size() > 1) nd.add_vertices(large_black_faces.size() - 1);
        std::size_t star = relv.size();
        for (HdsHalfedge* start : large_black_faces) {
            HdsHalfedge* j = start;
            HdsHalfedge* x = nd.add_edge(HdsPath({j, nullptr}), star, klim);
            for (;;) {
                j = j->next->twin->next; // jump over the spoke
                if (j == start) break;
                nd.add_edge(HdsPath({j, x}), star, klim);
            }
            ++star;
        }
        if (config_.intermediate_callback) config_.intermediate_callback(nd, canonical_code(nd, syms));
    }

    template <std::size_t klim>
    std::size_t LayerSearcher<klim>::run() {
        if (!config_.layer_edges) throw std::runtime_error("LayerSearcher: layer_edges not set");
        const bool pass1 = config_.pass1;
        const bool split = config_.split_depth > 0;
        const std::vector<CycleSym> syms = gadget_symmetries(config_.cycle_size, config_.layer_edges);

        // units finished by an earlier (killed) run of this task
        const std::string done_file = config_.output_dir + "/done_units.txt";
        std::set<std::size_t> done_units;
        if (split) {
            std::filesystem::create_directories(config_.output_dir);
            drop_partial_line(done_file);
            std::ifstream in(done_file);
            for (std::size_t u; in >> u;) done_units.insert(u);
        }
        std::size_t unit = 0;              // next unit number
        std::size_t open_unit = SIZE_MAX;  // owned unit being searched; complete once the DFS reaches the next one
        auto close_unit = [&]() {
            if (open_unit == SIZE_MAX) return;
            std::ofstream(done_file, std::ios::app) << open_unit << "\n";
            open_unit = SIZE_MAX;
        };
        auto owns = [&](std::size_t u) { return u % config_.task_count == config_.task_index && !done_units.count(u); };

        // the drawing with the vertices of the new cycle appended (only possible in a drawing without crossings)
        Drawing<klim> d = config_.input_drawing ? *config_.input_drawing : create_base_drawing();
        if (!d.crossings.empty()) throw std::runtime_error("LayerSearcher: the input drawing has crossings");
        const std::size_t nm = d.vertices.size(); // first new cycle vertex
        d.add_vertices(config_.cycle_size);
        std::vector<Edge> edges = config_.layer_edges(nm);

        // Pass 1 uses all checks after every edge; Pass 2 only the reachability and flow checks at the checkpoints
        auto is_checkpoint = [&](std::size_t i) {
            return std::find(config_.early_prune_checkpoints.begin(), config_.early_prune_checkpoints.end(), i)
                != config_.early_prune_checkpoints.end();
        };
        auto pcr_of = [&](const Edge& e) -> int { return (!pass1 && e.size() == 3) ? static_cast<int>(e[2]) : 0; };

        // placement filter (Pass 1): when an edge places a new vertex, skip every path that lands in a face the
        // vertex cannot be in (placements_still_possible for the later edges). Cached per position, valid as long
        // as the edge before it was not re-placed.
        struct PlacementFilter { bool valid = false; std::uint64_t stamp = 0; std::vector<int> face; std::vector<char> allowed; };
        std::vector<PlacementFilter> pfilter(edges.size());
        std::vector<std::uint64_t> placed_stamp(edges.size(), 0);
        std::uint64_t stamp_counter = 0;
        auto landing_ok = [&](std::size_t pos, const HdsPath& path, std::size_t target) -> bool {
            if (!pass1 || d.vertices[target].halfedge) return true;
            PlacementFilter& f = pfilter[pos];
            const std::uint64_t base = pos == 0 ? 0 : placed_stamp[pos - 1];
            if (!f.valid || f.stamp != base) {
                std::size_t num_faces = 0;
                std::map<std::size_t, std::vector<char>> domains;
                const bool possible = placements_still_possible(d, edges, pos + 1, &f.face, &num_faces, &domains);
                auto it = domains.find(target);
                if (!possible) f.allowed.assign(num_faces, 0);
                else if (it == domains.end()) f.allowed.assign(num_faces, 1);
                else f.allowed = it->second;
                f.valid = true;
                f.stamp = base;
            }
            // the new vertex lands in the face of path[0] (no crossing) or behind the last crossing
            const HdsHalfedge* h = path.size() == 2 ? path[0] : path[path.size() - 2]->twin;
            return f.allowed[f.face[h->label]];
        };
        auto next_path_filtered = [&](HdsPath& path, std::size_t target, int pcr, std::size_t pos) {
            while (d.next_path(path, target, pcr))
                if (landing_ok(pos, path, target)) return true;
            return false;
        };
        // does the drawing d, with the edge at position i just added, fail one of the checks?
        auto prunes = [&](std::size_t i) -> bool {
            if (i + 1 >= edges.size()) return false;
            const bool checkpoint = is_checkpoint(i);
            if (!pass1 && !checkpoint) return false;
            if (!check_remaining_edges_reachability(d, edges, i + 1)) return true;
            if (pass1 && !placements_still_possible(d, edges, i + 1)) return true;
            if (checkpoint && check_remaining_edges_mcf(d, edges, i + 1) == McfResult::Infeasible) return true;
            if (pass1 && i >= 1 && !final_face_still_possible(d, nm)) return true;
            return false;
        };

        // the order of the edges when an owned unit was entered: the reordering by the reachability check depends
        // on the search history, so it is undone when the DFS returns above the split depth, so that every task
        // numbers the units alike
        std::vector<Edge> split_snapshot;
        bool split_snapshot_valid = false;

        // depth-first search over the paths of the edges, in order
        auto e = edges.begin();
        for (;;) {
            std::size_t v = (*e)[1];
            int pcr = pcr_of(*e);
            HdsPath p = d.first_path((*e)[0], v, pcr);
            {
                const std::size_t pos = static_cast<std::size_t>(e - edges.begin());
                if (!p.empty() && !landing_ok(pos, p, v) && !next_path_filtered(p, v, pcr, pos)) p.clear();
            }
            if (p.empty()) {
BACKUP:
                for (;;) { // remove the last edge and try its next path
                    if (e == edges.begin()) goto FINISH;
                    --e;
                    v = (*e)[1];
                    pcr = pcr_of(*e);
                    p = d.edges.back().built;
                    d.remove_edge();
                    if (next_path_filtered(p, v, pcr, static_cast<std::size_t>(e - edges.begin()))) break;
                }
            }
            d.add_edge(p, v, pcr);

PLACED: // edge e has just been added
            {
                const std::size_t i = static_cast<std::size_t>(e - edges.begin());
                if (split_snapshot_valid && i + 1 <= config_.split_depth) { // back above the split depth
                    edges = split_snapshot;
                    split_snapshot_valid = false;
                }
                placed_stamp[i] = ++stamp_counter;
                while (prunes(i)) { // try the next paths of edge e
                    p = d.edges.back().built;
                    d.remove_edge();
                    if (!next_path_filtered(p, v, pcr, i)) goto BACKUP;
                    d.add_edge(p, v, pcr);
                    placed_stamp[i] = ++stamp_counter;
                }

                // split: only search below the units owned by this task
                if (split && i + 1 == config_.split_depth) {
                    close_unit(); // reaching the next unit completes the previous one
                    if (owns(unit)) {
                        open_unit = unit++;
                        split_snapshot = edges;
                        split_snapshot_valid = true;
                    } else {
                        ++unit;
                        p = d.edges.back().built;
                        d.remove_edge();
                        if (next_path_filtered(p, v, pcr, i)) {
                            d.add_edge(p, v, pcr);
                            goto PLACED;
                        }
                        goto BACKUP;
                    }
                }
            }

            if (++e == edges.end()) { // all edges drawn
                if (pass1) intermediate(d, nm, syms);
                else if (config_.solution_callback) config_.solution_callback(d);
                goto BACKUP;
            }
        }
FINISH:
        close_unit();
        return split ? unit : 1;
    }

} // namespace nested_search

#endif // NESTED_SEARCH_H

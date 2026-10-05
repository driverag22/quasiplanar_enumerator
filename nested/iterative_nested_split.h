#ifndef NESTED_CYCLE_SEARCH_HPP
#define NESTED_CYCLE_SEARCH_HPP

#include "hds_kplanar.h"
#include "iso.h"
#include <iostream>
#include <vector>
#include <deque>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <chrono>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <memory>
#include <map>
#include <nlohmann/json.hpp>
#include <boost/functional/hash.hpp>

// BGL Flow
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

namespace nested_cycle_build {

    // Max-Flow Dual Network Structure
    struct DualNetwork {
        typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;

        typedef boost::adjacency_list<
            boost::vecS, boost::vecS, boost::directedS,
            boost::no_property,
            boost::property<boost::edge_capacity_t, int,
            boost::property<boost::edge_residual_capacity_t, int,
            boost::property<boost::edge_reverse_t, Traits::edge_descriptor>>>
                > Graph;

        Graph G;
        DualNetwork(std::size_t num_nodes) : G(num_nodes) {}

        void add_edge(int from, int to, int capacity) {
            auto c_map = boost::get(boost::edge_capacity, G);
            auto r_map = boost::get(boost::edge_reverse, G);
            const auto e = boost::add_edge(from, to, G).first;
            const auto rev_e = boost::add_edge(to, from, G).first;
            c_map[e] = capacity;
            c_map[rev_e] = 0; // residual reverse capacity
            r_map[e] = rev_e;
            r_map[rev_e] = e;
        }

        int flow(int source, int sink) {return boost::push_relabel_max_flow(G, source, sink);}
    };

    // Same network as DualNetwork, for threshold tests "is the max flow from source to sink at
    // least limit?": flow(source, sink, limit) returns min(max flow, limit) using at most limit
    // BFS augmenting paths, which is much cheaper than a full max flow for the small limits
    // used here (cycle_size, 3).
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

    typedef std::vector<std::size_t> Edge;
    typedef std::vector<Edge> Edges;

    // Gadget-aware isomorphism
    //
    // The next extension step attaches the local edges to the active cycle vertices by label
    // (0,...,cycle_size-1). So two partial drawings are interchangeable only under isomorphisms
    // that map the active cycle onto itself by a permutation sym that is a symmetry of the local
    // edges. Such an isomorphism maps every extension of one drawing to an extension of the other.
    // Mirroring is always allowed (the mirror of an extension is an extension of the mirror).

    typedef std::vector<std::size_t> CycleSym; // sym[i] = image of active cycle vertex i
    typedef std::vector<std::uint32_t> IsoCode;

    struct IsoCodeHash {
        std::size_t operator()(const IsoCode& c) const { return boost::hash_range(c.begin(), c.end()); }
    };

    // All rotations/reflections sym of the old cycle 0,...,c-1 such that relabelling the old
    // endpoints by sym and the new cycle by some rotation/reflection maps the local edge set
    // (including the capacity entries) onto itself. The result is a group and contains the identity.
    //
    // Formally, computes the list of relabellings σ of cycle 0,...,c-1 s.t.:
    //   - σ(i) = i+s [rotation] or  σ(i) = s-i (mod c) [reflection]
    //   - some rot/refl π of the next cycle exists s.t. applying σ to old cycle and π to new maps 
    //     the local edge builder's edge list onto itself
    inline std::vector<CycleSym> gadget_symmetries(
            std::size_t c, const std::function<std::vector<Edge>(std::size_t)>& builder) {
        const std::size_t nm = c; // offset of the new cycle, any offset >= c works
        const std::vector<Edge> edges = builder(nm);
        for (const Edge& e : edges)
            for (std::size_t i = 0; i < 2; ++i)
                if (e[i] >= nm + c) throw std::runtime_error("gadget_symmetries: local edge out of range");

        typedef std::tuple<std::size_t, std::size_t, std::size_t> Key;
        auto edge_set = [&](const auto& map) { // returns edge set through given map
            std::set<Key> s;
            for (const Edge& e : edges) {
                std::size_t a = map(e[0]), b = map(e[1]);
                s.emplace(std::min(a, b), std::max(a, b), e.size() > 2 ? e[2] + 1 : 0);
            }
            return s;
        };
        auto dihedral = [c](std::size_t i, std::size_t shift, bool refl) {
            // rotation or reflection
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
                        auto map = [&](std::size_t x) {
                            // x < c takes old map, else take new map
                            return x < c ? sym[x] : nm + dihedral(x - nm, s_new, refl_new);
                        };
                        found = (edge_set(map) == original);
                    }
                if (found) syms.push_back(std::move(sym));
            }
        return syms;
    }

    // Pairwise test: true iff there is an isomorphism d1 -> d2 (possibly mirrored) of the
    // planarizations that preserves vertex/crossing type and ncr, and acts on the active
    // cycle vertices 0,...,c-1 as one of the given symmetries.
    // PRE: Both drawings are connected apart from isolated vertices, and have at least c = sym.size() vertices.
    template <int kplane>
    bool are_isomorphic_gadget(const Drawing<kplane>& d1, const Drawing<kplane>& d2,
            const std::vector<CycleSym>& syms) {
        const std::size_t nv = d1.vertices.size();
        if (nv != d2.vertices.size()) return false;
        if (d1.crossings.size() != d2.crossings.size()) return false;
        if (d1.edges.size() != d2.edges.size()) return false;
        const std::size_t n = nv + d1.crossings.size();

        // try to extend the fixed map on the cycle, aligning halfedge s1 (to vertex 0) with s2
        auto try_map = [&](const CycleSym& sym, bool mirror, const HdsHalfedge* s1, const HdsHalfedge* s2) {
            std::vector<std::size_t> phi(n, n), psi(n, n); // d1 -> d2 and inverse
            for (std::size_t i = 0; i < sym.size(); ++i) { phi[i] = sym[i]; psi[sym[i]] = i; }
            std::vector<char> queued(n, 0);
            queued[0] = 1;
            std::queue<std::pair<const HdsHalfedge*, const HdsHalfedge*>> q;
            q.emplace(s1, s2);
            while (!q.empty()) {
                auto [a1, a2] = q.front();
                q.pop();
                const HdsHalfedge* t1 = a1;
                const HdsHalfedge* t2 = a2;
                do {
                    if (t1->edge->ncr != t2->edge->ncr) return false;
                    std::size_t i1 = t1->twin->vertex->label;
                    std::size_t i2 = t2->twin->vertex->label;
                    if ((i1 >= nv) != (i2 >= nv)) return false; // vertices to vertices, crossings to crossings
                    if (phi[i1] == n) {
                        if (psi[i2] != n) return false;
                        phi[i1] = i2;
                        psi[i2] = i1;
                    } else if (phi[i1] != i2) return false;
                    if (!queued[i1]) { queued[i1] = 1; q.emplace(t1->twin, t2->twin); }
                    t1 = (mirror ? t1->twin->prev : t1->next->twin);
                    t2 = t2->next->twin;
                } while (t1 != a1 && t2 != a2);
                if (t1 != a1 || t2 != a2) return false; // degrees differ
            }
            for (std::size_t z = 0; z < n; ++z) {
                if (z < nv && !d1.vertices[z].halfedge) continue; // isolated, no structure
                if (phi[z] == n || psi[phi[z]] != z)
                    throw std::runtime_error("are_isomorphic_gadget: disconnected drawing");
            }
            return true;
        };

        // isolated vertices (e.g. the reserved star vertex of an nd without large black faces)
        // can be mapped onto each other arbitrarily, only their number has to agree
        auto isolated = [](const Drawing<kplane>& d) {
            std::size_t k = 0;
            for (const auto& v : d.vertices) if (!v.halfedge) ++k;
            return k;
        };
        if (isolated(d1) != isolated(d2)) return false;

        for (const CycleSym& sym : syms)
            for (int mirror = 0; mirror < 2; ++mirror) {
                const HdsHalfedge* s1 = d1.vertices[0].halfedge;
                const HdsHalfedge* h2 = d2.vertices[sym[0]].halfedge;
                if (!s1 || !h2) throw std::runtime_error("are_isomorphic_gadget: isolated cycle vertex");
                const HdsHalfedge* s2 = h2;
                do {
                    if (try_map(sym, mirror, s1, s2)) return true;
                    s2 = s2->next->twin;
                } while (s2 != h2);
            }
        return false;
    }

    // Encode d by a BFS from halfedge h0 (pointing to the cycle vertex that sym maps to 0), with
    // cycle vertex i labelled sym[i] and all other vertices/crossings labelled c, c+1, ... in
    // order of discovery; rotations are reversed if mirror is set. The code is compared with best
    // while it is built: stop as soon as it is larger, replace best if it ends up smaller.
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
        // helper function to add to cur code
        auto emit = [&](std::size_t x) {
            if (!smaller) {
                if (x > best[cur.size()]) return false;
                if (x < best[cur.size()]) smaller = true;
            }
            cur.push_back(static_cast<std::uint32_t>(x));
            return true;
        };

        // isolated vertices (e.g. the reserved star vertex of an nd without large black faces)
        // carry no structure, only their number is recorded
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
                // 1/0 bit determines crossing or not, so even labels are vertices, odd are crossings
                if (!emit(2 * phi[w] + (w >= nv ? 1 : 0)) || !emit(t->edge->ncr)) return;
                t = (mirror ? t->twin->prev : t->next->twin);
            } while (t != s);
            if (!emit(SEP)) return;
        }
        if (nqueued + isolated != n) throw std::runtime_error("canonical_code: disconnected drawing");
        if (smaller) best.swap(cur);
    }

    // Canonical code of d with respect to the gadget symmetries syms.
    //
    // For each sym in syms, each orientation (mirror = 0/1), and each halfedge h
    // pointing to the cycle vertex u0 with sym[u0] == 0, gets encoding 
    // from encode_if_smaller.
    // The code is the lexicographically smallest of these encodings.
    //
    // Two drawings get the same code iff there is an isomorphism of their
    // planarizations (possibly mirrored) that preserves vertex/crossing type and
    // ncr and acts on the active cycle 0,...,c-1 as an element of syms. 
    // The "if" direction requires syms to be the group returned by gadget_symmetries.
    //
    // Cost: at most 2 * |syms| * deg(u0) BFS runs, each O(#halfedges); most runs
    // stop after a few entries because the code is compared with the current
    // best while it is built.
    //
    // PRE: d is connected apart from isolated vertices, and has at least c = syms[0].size() vertices.
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

    // Canonical code for strong isomorphism without fixed labels (as are_isomorphic in iso.h):
    // equal codes iff the planarizations are isomorphic (possibly mirrored), preserving vertex/crossing
    // type and ncr. Used to dedup complete (Pass 2) drawings. Minimum over all start halfedges pointing
    // to a vertex and both orientations.
    // PRE: d is connected apart from isolated vertices.
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

    template <std::size_t klim>
    inline void check_and_terminate_if_invalid(const Drawing<klim>& d, const std::string& filename = "../quasiDrawings/failExample.graphml") {
        // 1. Check for self-loops (u == v)
        for (const auto& edge : d.edges) {
            if (edge.u == edge.v) {
                std::cout << "\n[INVALID DRAWING] Self-loop detected at vertex " << edge.u << "!" << std::endl;
                std::ofstream of_graphml(filename);
                d.graphml_output(of_graphml);
                of_graphml.close();
                std::cerr << "Saved invalid drawing to " << filename << ". Terminating execution." << std::endl;
                std::exit(1);
            }
        }

        // 2. Check for parallel edges (duplicate {u, v} pairs)
        std::set<std::pair<std::size_t, std::size_t>> seen_edges;
        for (const auto& edge : d.edges) {
            std::size_t u = std::min(edge.u, edge.v);
            std::size_t v = std::max(edge.u, edge.v);
            if (seen_edges.count({u, v}) > 0) {
                std::cout << "\n[INVALID DRAWING] Parallel edge detected between vertices " << u << " and " << v << "!" << std::endl;
                std::ofstream of_graphml(filename);
                d.graphml_output(of_graphml);
                of_graphml.close();
                std::cerr << "Saved invalid drawing to " << filename << ". Terminating execution." << std::endl;
                std::exit(1);
            }
            seen_edges.insert({u, v});
        }
    }

    // Returns a halfedge pointing to the desired source vertex, and the number 
    // of crossings the edge should have in the reduced drawing.
    //
    // Used when a halfedge has an irrelevant endpoint:
    //   vact[v] == vact.size()   iff   v is irrelevant
    template <std::size_t klim, typename F, typename A, typename V>
        inline std::pair<const HdsHalfedge*, std::size_t> prev_active(
                const HdsHalfedge* e,
                const F& face,
                const A& active,
                const V& vact) 
        {
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

    template <std::size_t klim = 3>
        class NestedCycleSearcher {
            public:
                // Queue node wrapping drawing and explicit depth level
                struct SearchNode {
                    Drawing<klim> drawing;
                    std::size_t depth;
                };
                // Configuration options for running a search instance
                struct Config {
                    std::size_t cycle_size = 12;
                    std::size_t max_cycles = 2;
                    std::string output_dir = "./output_drawings";
                    bool export_files = true;
                    bool verbose = true;
                    bool enable_early_pruning = true;

                    // integral multicommodity flow check at the early prune checkpoints
                    bool enable_mcf_pruning = true;
                    std::size_t mcf_node_budget = 200000; // backtracking nodes before giving up (no prune)
                    std::size_t mcf_max_paths = 4096;     // commodities with more candidate paths are dropped
                    bool mcf_log_progress = true;         // periodically print running MCF prune count

                    // debug: cross-check every canonical-code lookup against pairwise are_isomorphic_gadget
                    // (and count matches that only the label-agnostic are_isomorphic of iso.h would find)
                    bool verify_iso = false;

                    // debug: also run Boost's push_relabel_max_flow on every threshold flow test of the
                    // leaf and stop if the results ever disagree
                    bool verify_threshold_flow = false;

                    // Splitting the search over independent processes (e.g. a cluster job array).
                    // With split_mode set, run() only processes the given parents (in order), does not
                    // extend the drawings it finds, and exports them (one JSON per line in
                    // output_dir/children.jsonl, a single file to keep the file count low) for a
                    // later merge. Parents with at least one valid child (before dedup) are appended to
                    // output_dir/extensible.txt.
                    bool split_mode = false;
                    std::vector<std::string> input_parents;  // JSON files, or "base" for the base drawing
                    bool run_pass1 = true, run_pass2 = true; // split_mode only (replaces the Pass 2 gating)

                    // Prefix ownership: the paths of edge split_depth-1 that survive pruning are numbered
                    // 0,1,2,... over all parents and passes (identically in every process), and only those
                    // with index % task_count == task_index are searched further. split_depth == 0: one
                    // unit per (parent, pass). task_count == 0: only count the units, search nothing.
                    std::size_t split_depth = 0;
                    std::size_t task_index = 0;
                    std::size_t task_count = 1;

                    // split_mode memory: children are not kept in solutions (only their codes, for local
                    // dedup), and the codes are dropped once there are more than split_dedup_limit of
                    // them; this only causes extra duplicate exports, which the merge removes.
                    // split_keep_children keeps them in solutions anyway (for tests).
                    std::size_t split_dedup_limit = 100000;

                    // print a progress line (counters) every progress_seconds; 0 = off
                    double progress_seconds = 0;
                    bool split_keep_children = false;
                    bool split_export_children = true; // false: don't write children.jsonl (benchmarks)

                    // custom edge index checkpoints where early flow checks run.
                    std::vector<std::size_t> early_prune_checkpoints;

                    // run the reachability check (not MCF) also after every edge with index >= this one;
                    // remaining edges with an unplaced endpoint are skipped there. Sound (same relaxation as
                    // at the checkpoints), it just finds dead ends earlier. SIZE_MAX = off.
                    std::size_t reach_every_edge_from = SIZE_MAX;

                    // Pass 1 only: after every edge with index >= this one, once all new cycle vertices are
                    // placed, prune if no face can still become a possible final face (see
                    // final_face_still_possible). SIZE_MAX = off.
                    std::size_t final_face_every_edge_from = SIZE_MAX;

                    // after every edge with index >= this one: prune if some unplaced vertex has no face left
                    // that all its remaining edges can reach (see placements_still_possible). SIZE_MAX = off.
                    std::size_t placement_every_edge_from = SIZE_MAX;

                    // draw the gadget edges that connect to the old drawing (the matching) before the others
                    // (the new cycle). The set of completions does not depend on the edge order.
                    bool matching_edges_first = false;

                    // placement filter: when an edge places a new vertex, skip (without adding the edge or
                    // running any check) every path that lands in a face outside the vertex's possible faces
                    // (placements_still_possible for the remaining edges). Sound for the same reason; it only
                    // avoids generating placements that would be pruned right after being added.
                    bool placement_filter = false;

                    // Generator function for local edges added during extension step.
                    // Takes parent vertex offset (nm) and cycle size, returns vector of edges {u, v, [capacity]}.
                    std::function<std::vector<Edge>(std::size_t nm)> local_edges_builder;
                };

                // Search summary stats
                struct SearchResult {
                    std::vector<SearchNode> solutions;
                    std::size_t discarded_count = 0;
                    std::size_t full_solution_count = 0;
                    std::size_t total_processed = 0;
                    std::size_t pruned_early_count = 0;
                    std::size_t pruned_mcf_count = 0;   // subset of pruned_early_count caused by the MCF check
                    std::size_t pruned_final_face_count = 0; // subset caused by final_face_still_possible
                    std::size_t pruned_placement_count = 0;  // subset caused by placements_still_possible
                    std::size_t filtered_paths_count = 0;   // paths skipped by the placement filter
                    std::size_t leaf_count = 0;          // complete Pass 1 drawings reached
                    std::size_t units_done = 0;          // split_mode: units completed in this run
                    std::size_t mcf_unknown_count = 0;  // MCF checks that hit the node budget
                    double mcf_seconds = 0;             // total time spent in MCF checks
                    double iso_seconds = 0;             // total time spent in isomorphism dedup
                    std::size_t gadget_symmetry_count = 0;
                    // verify_iso only: new drawings that the label-agnostic test would have discarded
                    std::size_t label_agnostic_only_count = 0;
                    std::size_t split_units = 0;        // units enumerated (owned or not)
                    std::size_t children_exported = 0;  // split_mode only
                };

                enum class McfResult { Feasible, Infeasible, Unknown };

                // Constructors
                NestedCycleSearcher() = default;
                explicit NestedCycleSearcher(Config config) : config_(std::move(config)) {}


                Drawing<klim> create_base_drawing() const {
                    std::size_t cycle_size = config_.cycle_size;
                    // cycle_size vertices for C_{cycle_size} + 1 dummy star vertex = cycle_size+1 vertices
                    Drawing<klim> d(cycle_size + 1);
                    std::vector<HdsHalfedge*> cycle(cycle_size, nullptr);

                    // build uncrossable C_{cycle_size} cycle by setting ncr=klim=3
                    cycle[1] = d.add_first_edge(0, 1, klim);
                    for (std::size_t i = 2; i < cycle_size; ++i) {
                        cycle[i] = d.add_edge(HdsPath({cycle[i - 1], nullptr}), i, klim);
                    }
                    cycle[0] = d.add_edge(HdsPath({cycle[cycle_size - 1], cycle[1]->twin}), 0, klim);

                    // block one side of the cycle using an uncrossable star at dummy vertex cycle_size
                    std::size_t star_center = cycle_size;
                    auto e = d.add_edge(HdsPath({cycle[0], nullptr}), star_center, klim);
                    for (std::size_t i = cycle_size - 1; i > 0; --i) {
                        e = d.add_edge(HdsPath({cycle[i], e}), star_center, klim);
                    }

                    return d;
                }

                // map halfedges to face indices, returns number of faces
                static std::size_t label_faces(const Drawing<klim>& d, std::vector<int>& face) {
                    face.assign(d.halfedges.size(), -1);
                    std::size_t num_faces = 0;
                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                        if (face[i->label] != -1) continue;
                        const HdsHalfedge* j = &*i;
                        do { // walk around face
                            face[j->label] = static_cast<int>(num_faces);
                            j = j->next;
                        } while (j->label != i->label);
                        ++num_faces;
                    }
                    return num_faces;
                }

                // all face IDs incident to a vertex (empty if the vertex is isolated)
                static std::vector<int> incident_faces(
                        const Drawing<klim>& d,
                        const std::vector<int>& face,
                        std::size_t v_label) {
                    std::vector<int> faces;
                    if (v_label >= d.vertices.size()) return faces;
                    // halfedge pointing to v_label
                    auto start_h = d.vertices[v_label].halfedge;
                    if (!start_h) return faces;

                    auto curr_h = start_h;
                    do {
                        int f = face[curr_h->label];
                        // make sure face is valid and not already in faces
                        if (f >= 0) if (std::find(faces.begin(), faces.end(), f) == faces.end())
                            faces.push_back(f);
                        curr_h = curr_h->next->twin; // walk around halfedges incident to vertex
                    } while (curr_h != start_h);
                    return faces;
                }

                // Integral multicommodity flow relaxation for the remaining edges: every remaining edge
                // (commodity) picks one sequence of crossed drawn edges, such that no drawn edge exceeds
                // its remaining capacity klim - ncr. Crossings among the remaining edges are ignored, so
                // Infeasible is a proof that no completion exists. Unknown means the node budget was hit.
                McfResult check_remaining_edges_mcf(
                        const Drawing<klim>& d,
                        const std::vector<Edge>& local_edges,
                        std::size_t next_edge_index,
                        int constrained
                        ) const {
                    std::vector<int> face;
                    std::size_t num_faces = label_faces(d, face);
                    if (num_faces == 0) return McfResult::Unknown;

                    // remaining crossing capacity per drawn edge (edge labels are 0 ... edges.size()-1)
                    std::vector<int> cap(d.edges.size(), 0);
                    for (const auto& e : d.edges) cap[e.label] = static_cast<int>(klim) - static_cast<int>(e.ncr);

                    // dual graph: face -> (adjacent face, crossed edge)
                    std::vector<std::vector<std::pair<int, const HdsEdge*>>> dual(num_faces);
                    for (const auto& h : d.halfedges) {
                        if (h.edge->ncr >= klim) continue;
                        int f1 = face[h.label], f2 = face[h.twin->label];
                        if (f1 != f2) dual[f1].push_back({f2, h.edge});
                    }

                    typedef std::vector<std::size_t> McfPath; // sorted labels of crossed edges
                    std::vector<std::vector<McfPath>> paths;  // candidate paths per commodity
                    std::vector<char> is_target(num_faces), in_path(num_faces);
                    McfPath cur;

                    for (std::size_t rem_idx = next_edge_index; rem_idx < local_edges.size(); ++rem_idx) {
                        const Edge& le = local_edges[rem_idx];
                        std::size_t u = le[0], v = le[1];
                        std::vector<int> u_faces = incident_faces(d, face, u);
                        std::vector<int> v_faces = incident_faces(d, face, v);
                        // isolated endpoint: can be placed in any face, drop commodity
                        if (u_faces.empty() || v_faces.empty()) continue;

                        std::fill(is_target.begin(), is_target.end(), 0);
                        for (int f : v_faces) is_target[f] = 1;
                        bool free_edge = false;
                        for (int f : u_faces) if (is_target[f]) free_edge = true;
                        if (free_edge) continue; // uncrossed drawing possible, uses no capacity

                        // paths through another face of u are dominated by starting there
                        std::fill(in_path.begin(), in_path.end(), 0);
                        for (int f : u_faces) in_path[f] = 1;

                        std::vector<McfPath> cp;
                        bool overflow = false;
                        auto dfs = [&](auto& self, int f, std::size_t depth) -> void {
                            if (depth >= klim) return;
                            for (const auto& [nb, e] : dual[f]) {
                                if (overflow) return;
                                if (in_path[nb]) continue;
                                // same restrictions as Drawing::find_crossing
                                if (e->u == u || e->v == u || e->u == v || e->v == v) continue;
                                if (std::find(cur.begin(), cur.end(), e->label) != cur.end()) continue;
                                cur.push_back(e->label);
                                if (is_target[nb]) {
                                    // stop here: any extension is a dominated superset
                                    cp.push_back(cur);
                                    std::sort(cp.back().begin(), cp.back().end());
                                    if (cp.size() > config_.mcf_max_paths) overflow = true;
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

                    // backtracking with minimum-remaining-values commodity selection
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
                        if (++nodes > config_.mcf_node_budget) { aborted = true; return false; }

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

                bool check_remaining_edges_reachability(
                        const Drawing<klim>& d,
                        std::vector<Edge>& local_edges,
                        std::size_t next_edge_index
                        ) const {
                    std::vector<int> face;
                    std::size_t num_faces = label_faces(d, face);
                    if (num_faces == 0) return false;

                    // dual adjacency list filtered by available crossing capacity (capacity > 0):
                    // face -> (adjacent face, crossed edge)
                    std::vector<std::vector<std::pair<int, const HdsEdge*>>> dual_adj(num_faces);
                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                        int remaining_capacity = static_cast<int>(klim) - static_cast<int>(i->edge->ncr);
                        if (remaining_capacity > 0) {
                            int f1 = face[i->label];
                            int f2 = face[i->twin->label];
                            if (f1 >= 0 && f2 >= 0 && f1 != f2)
                                dual_adj[f1].push_back({f2, i->edge});
                        }
                    }
                    std::vector<std::size_t> crossed; // labels of edges crossed on the current dfs path

                    std::vector<std::pair<std::size_t, Edge>> rem_edges_with_counts;
                    rem_edges_with_counts.reserve(local_edges.size() - next_edge_index);

                    std::vector<bool> is_target(num_faces, false);
                    std::vector<bool> in_path(num_faces, false);

                    for (std::size_t rem_idx = next_edge_index; rem_idx < local_edges.size(); ++rem_idx) {
                        std::size_t u = local_edges[rem_idx][0];
                        std::size_t v = local_edges[rem_idx][1];
                        std::vector<int> u_faces = incident_faces(d, face, u);
                        std::vector<int> v_faces = incident_faces(d, face, v);

                        // an edge with an unplaced endpoint can always be drawn (the endpoint goes into a
                        // suitable face), so it is not checked and keeps its relative order among the
                        // unordered (>3 paths) edges, i.e. after every edge that places its endpoints.
                        // (This used to return false as a "fail safe"; with the original edge order all
                        // endpoints are placed at the checkpoints, but once edges are reordered before a
                        // checkpoint that pruned valid branches: on C6 only 9567 of 401026 leaves remained.)
                        if (u_faces.empty() || v_faces.empty()) {
                            rem_edges_with_counts.push_back({4, local_edges[rem_idx]});
                            continue;
                        }

                        std::fill(is_target.begin(), is_target.end(), false);
                        for (int f_v : v_faces) is_target[f_v] = true;

                        std::size_t path_count = 0;
                        std::fill(in_path.begin(), in_path.end(), false);

                        // perform dfs for each edge independently, determine number of ways to draw it
                        auto dfs_path_count = [&](auto& self, int curr_face, std::size_t depth) -> void {
                            if (path_count > 3) return;
                            if (is_target[curr_face]) {
                                path_count++;
                                return;
                            }
                            if (depth >= klim) return;
                            in_path[curr_face] = true;
                            for (const auto& [neighbor, ce] : dual_adj[curr_face]) {
                                if (in_path[neighbor]) continue;
                                // same restrictions as Drawing::find_crossing: no crossing of edges
                                // sharing an endpoint with uv, no crossing the same edge twice
                                if (ce->u == u || ce->v == u || ce->u == v || ce->v == v) continue;
                                if (std::find(crossed.begin(), crossed.end(), ce->label) != crossed.end()) continue;
                                crossed.push_back(ce->label);
                                self(self, neighbor, depth + 1);
                                crossed.pop_back();
                                if (path_count > 3) break;
                            }
                            in_path[curr_face] = false;
                        };

                        for (int f_u : u_faces) {
                            dfs_path_count(dfs_path_count, f_u, 0);
                            if (path_count > 3) break;
                        }

                        if (path_count == 0) return false;

                        rem_edges_with_counts.push_back({path_count, local_edges[rem_idx]});
                    }

                    // stably sort remaining edges: edges with fewer paths (1, 2, 3) are prioritized first.
                    // Edges with >3 paths maintain their original order at the end.
                    std::stable_sort(
                        rem_edges_with_counts.begin(),
                        rem_edges_with_counts.end(),
                        [](const auto& a, const auto& b) {
                            return a.first < b.first;
                        }
                    );

                    for (std::size_t i = 0; i < rem_edges_with_counts.size(); ++i) {
                        local_edges[next_edge_index + i] = rem_edges_with_counts[i].second;
                    }

                    return true;
                }


                // Early version of the leaf test "some face f has flow >= cycle_size from the active cycle
                // vertices" (potential final face), for a partial drawing whose new cycle vertices
                // nm, ..., nm+cycle_size-1 are all placed. Returns false only if no completion can have a
                // potential final face, so pruning on false is sound:
                //  - faces of a completion refine the current faces; contracting them maps every flow of the
                //    leaf network onto this network (arcs of new edges become loops and are dropped),
                //  - a current edge piece with x crossings (ncr, incl. prescribed ones) can be split into j
                //    pieces later, each with leaf capacity klim - ncr_final <= klim - x - (j-1), so all its
                //    descendants together carry at most bound(x) = max_j j*(klim-x-j+1) per direction
                //    (klim = 3: 4, 2, 1, 0 for x = 0, 1, 2, 3); the leaf capacity klim - x alone is not
                //    enough, since crossing a piece can increase the total capacity across it,
                //  - every later angle at a cycle vertex lies inside a current angle, and each cycle vertex
                //    sends at most 1 unit anyway.
                // Placement check for the remaining gadget edges (constraint propagation). Every unplaced
                // vertex x will lie inside some current face (faces are only split later). A remaining edge
                // with at most L crossings (L = klim - prescribed crossings) between x and a placed vertex p
                // forces x into a face within dual distance L of p, using only crossable edges (ncr < klim)
                // not incident to p, as in check_remaining_edges_reachability; an edge between two unplaced
                // vertices forces their faces within distance L of each other. The possible faces of all
                // unplaced vertices are narrowed until nothing changes; if some vertex has none left, no
                // completion exists. Only necessary conditions are used, so pruning on false is sound.
                // (Found on dead level-1 parents: the new cycle walks around and cannot be closed at the end.)
                // Optional outputs (for the placement filter): face labels of d, the number of faces, and the
                // narrowed possible faces of the unplaced vertices (only meaningful if true is returned).
                bool placements_still_possible(const Drawing<klim>& d, const std::vector<Edge>& local_edges,
                        std::size_t next_edge_index, int constrained,
                        std::vector<int>* face_out = nullptr, std::size_t* num_faces_out = nullptr,
                        std::map<std::size_t, std::vector<char>>* domains_out = nullptr) const {
                    std::vector<int> face;
                    const std::size_t num_faces = label_faces(d, face);
                    if (face_out) *face_out = face;
                    if (num_faces_out) *num_faces_out = num_faces;
                    // dual arcs over crossable edges: face -> (adjacent face, crossed edge)
                    std::vector<std::vector<std::pair<int, const HdsEdge*>>> dual(num_faces);
                    for (const auto& h : d.halfedges)
                        if (h.edge->ncr < klim && face[h.label] != face[h.twin->label])
                            dual[face[h.label]].push_back({face[h.twin->label], h.edge});

                    // faces within distance limit of the start set, not crossing edges incident to avoid
                    // (avoid = SIZE_MAX: no restriction)
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
                    for (std::size_t i = next_edge_index; i < local_edges.size(); ++i) {
                        const Edge& e = local_edges[i];
                        const std::size_t pcr = (constrained == 1 && e.size() == 3) ? e[2] : 0;
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
                    // narrow along the links until nothing changes
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

                // Crossing counts of the gadget edges in a complete (unreduced) drawing d whose new cycle is
                // nm, ..., nm+cycle_size-1: "cycle_crossings"[i] for the cycle edge {nm+i, nm+i+1} and
                // "matching_crossings"[i] for the edge from old vertex i (-1 if absent). Includes prescribed
                // crossings (Pass 2 cycle edges have klim). The expected drawing has all cycle edges at 0 and
                // all matching edges at klim.
                nlohmann::ordered_json gadget_crossings(const Drawing<klim>& d, std::size_t nm) const {
                    const std::size_t c = config_.cycle_size;
                    std::vector<int> cycle(c, -1), matching(c, -1);
                    for (const auto& e : d.edges) {
                        const std::size_t a = std::min(e.u, e.v), b = std::max(e.u, e.v);
                        if (a >= nm) {
                            for (std::size_t i = 0; i < c; ++i) {
                                const std::size_t x = nm + i, y = nm + (i + 1) % c;
                                if (a == std::min(x, y) && b == std::max(x, y)) cycle[i] = static_cast<int>(e.ncr);
                            }
                        } else if (b >= nm && a < c) {
                            matching[a] = static_cast<int>(e.ncr);
                        }
                    }
                    nlohmann::ordered_json m;
                    m["cycle_crossings"] = cycle;
                    m["matching_crossings"] = matching;
                    return m;
                }

                bool final_face_still_possible(const Drawing<klim>& d, std::size_t nm) const {
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

                SearchResult run() {
                    if (!config_.local_edges_builder)
                        throw std::runtime_error("NestedCycleSearcher Error: local_edges_builder callback is not set!");

                    if (config_.export_files)
                        std::filesystem::create_directories(config_.output_dir);


                    SearchResult result;
                    // split_mode: parents are only referenced here and loaded one at a time when processed
                    // (there can be far too many to keep in memory). An entry of input_parents is "base",
                    // a JSON file with one drawing, or a .jsonl file with one drawing per line.
                    struct ParentRef {
                        std::string path;        // "base" for the base drawing
                        std::streamoff offset;   // .jsonl: start of the line, otherwise -1
                        std::size_t line;        // .jsonl: line index (for messages)
                    };
                    std::vector<ParentRef> parent_refs;
                    if (config_.split_mode) {
                        std::filesystem::create_directories(config_.output_dir);
                        std::vector<std::string> inputs = config_.input_parents;
                        if (inputs.empty()) inputs.push_back("base");
                        for (const std::string& path : inputs) {
                            if (path.size() > 6 && path.compare(path.size() - 6, 6, ".jsonl") == 0) {
                                std::ifstream in(path, std::ios::binary);
                                if (!in) throw std::runtime_error("cannot open parents " + path);
                                std::size_t line_index = 0;
                                for (std::string line;;) {
                                    std::streamoff pos = in.tellg();
                                    if (!std::getline(in, line)) break;
                                    if (line.size() > 1) parent_refs.push_back({path, pos, line_index});
                                    ++line_index;
                                }
                            } else {
                                parent_refs.push_back({path, -1, 0});
                            }
                        }
                    } else {
                        result.solutions.push_back({create_base_drawing(), 0});
                    }
                    auto load_parent = [&](const ParentRef& r) -> Drawing<klim> {
                        if (r.path == "base") return create_base_drawing();
                        std::ifstream in(r.path, std::ios::binary);
                        if (!in) throw std::runtime_error("cannot open parent " + r.path);
                        if (r.offset < 0) return Drawing<klim>(nlohmann::json::parse(in));
                        in.seekg(r.offset);
                        std::string line;
                        std::getline(in, line);
                        return Drawing<klim>(nlohmann::json::parse(line));
                    };
                    auto describe_parent = [&](const ParentRef& r) {
                        return r.offset < 0 ? r.path : r.path + "#" + std::to_string(r.line);
                    };
                    // split_mode: only the input parents are processed, found children are not extended
                    const std::size_t num_inputs = config_.split_mode ? parent_refs.size() : result.solutions.size();
                    std::size_t unit = 0; // see Config::split_depth

                    // resuming: units completed by an earlier (killed) run of this task are listed in
                    // output_dir/done_units.txt and skipped; new children are appended to children.jsonl
                    // (a line cut off by a kill is skipped by the merge)
                    std::set<std::size_t> done_units;
                    if (config_.split_mode) {
                        std::ifstream in(config_.output_dir + "/done_units.txt");
                        for (std::size_t u; in >> u;) done_units.insert(u);
                        // a kill may have cut off the last child line: start new children on a fresh line
                        std::ifstream ch(config_.output_dir + "/children.jsonl", std::ios::ate | std::ios::binary);
                        if (ch && ch.tellg() > 0) {
                            ch.seekg(-1, std::ios::end);
                            if (ch.get() != '\n') std::ofstream(config_.output_dir + "/children.jsonl", std::ios::app) << "\n";
                        }
                    }
                    auto owns = [&](std::size_t u) {
                        return config_.task_count != 0 && u % config_.task_count == config_.task_index
                            && !done_units.count(u);
                    };
                    // the owned unit currently being searched; it is complete once the DFS reaches the
                    // next unit (owned or not) or finishes the pass
                    std::size_t open_unit = SIZE_MAX;
                    auto close_unit = [&]() {
                        if (open_unit == SIZE_MAX) return;
                        std::ofstream out(config_.output_dir + "/done_units.txt", std::ios::app);
                        out << open_unit << "\n";
                        ++result.units_done;
                        open_unit = SIZE_MAX;
                    };
                    std::set<std::size_t> reported_extensible; // parents already in extensible.txt
                    auto last_mcf_log = std::chrono::steady_clock::now();
                    const auto run_start = std::chrono::steady_clock::now();
                    auto last_progress = run_start;

                    // symmetries of the local edges on the active cycle, and canonical codes of all
                    // solutions found so far (code -> index into result.solutions)
                    const std::vector<CycleSym> syms = gadget_symmetries(config_.cycle_size, config_.local_edges_builder);
                    result.gadget_symmetry_count = syms.size();
                    std::unordered_map<IsoCode, std::size_t, IsoCodeHash> seen_codes;
                    for (std::size_t i = 0; i < result.solutions.size(); ++i)
                        seen_codes.emplace(canonical_code(result.solutions[i].drawing, syms), i);
                    if (config_.verbose)
                        std::cout << "Gadget symmetries on the active cycle: " << syms.size() << std::endl;

                    // we process the solutions in order, this is the current one we process
                    std::size_t solcount = 0;

                    while (solcount < (config_.split_mode ? num_inputs : result.solutions.size())) {
                        if (config_.split_mode && config_.split_depth == 0) {
                            // one unit per (parent, pass): skip parents without an owned unit before loading them
                            const std::size_t passes = (config_.run_pass1 ? 1 : 0) + (config_.run_pass2 ? 1 : 0);
                            bool any_owned = false;
                            for (std::size_t k = 0; k < passes; ++k) if (owns(unit + k)) any_owned = true;
                            if (!any_owned) {
                                unit += passes;
                                result.split_units += passes;
                                ++solcount;
                                continue;
                            }
                        }
                        Drawing<klim> d_parent = config_.split_mode ? load_parent(parent_refs[solcount]) : result.solutions[solcount].drawing;
                        const std::size_t current_depth = config_.split_mode ? 0 : result.solutions[solcount].depth;

                        if (config_.verbose) {
                            std::cout << "\n--- Processing Partial Drawing #" << solcount 
                                << " [Cycle Level " << current_depth << "]"
                                << " (Queue Size: " << result.solutions.size() << ") ---" << std::endl;
                        }

                        if (!config_.split_mode && current_depth >= config_.max_cycles) {
                            if (config_.verbose) {
                                std::cout << "--> Reached MAX_CYCLES limit (" << config_.max_cycles 
                                    << "). Skipping further cycle extension for Drawing #" << solcount << std::endl;
                            }
                            ++solcount;
                            continue;
                        }

                        bool is_extensible = false;
                        // Two-Pass Loop
                        // Pass 0 (constrained == 0): Unconstrained extension & subdrawing extraction
                        // Pass 1 (constrained == 1): Constrained (uncrossable) decagon completion
                        for (int constrained = 0; constrained < 2; ++constrained) {

                            // Pass 2 (constrained == 1) only runs if Pass 1 proved extensibility
                            if (config_.split_mode) {
                                if (!(constrained == 0 ? config_.run_pass1 : config_.run_pass2)) continue;
                                if (config_.split_depth == 0) {
                                    ++result.split_units;
                                    if (!owns(unit)) { ++unit; continue; }
                                    open_unit = unit++;
                                }
                            } else if (constrained == 1 && !is_extensible) {
                                if (config_.verbose) {
                                    std::cout << "--> Skipping Constrained Pass 2 (Drawing #" << solcount 
                                        << " is not extensible)" << std::endl;
                                }
                                continue;
                            }

                            if (config_.verbose) {
                                std::cout << "\n>>> Starting Pass " << (constrained == 0 ? "1 (Unconstrained)" : "2 (Constrained)")
                                    << " for Drawing #" << solcount << " <<<" << std::endl;
                            }

                            Drawing<klim> d = d_parent;
                            std::size_t nm = d.vertices.size();
                            d.add_vertices(config_.cycle_size);
                            std::vector<Edge> local_edges = config_.local_edges_builder(nm);
                            // placement filter cache per position: allowed landing faces of the new vertex placed by
                            // the edge at that position; valid while the edge before it was not re-placed
                            struct PlacementFilter { bool valid = false; std::uint64_t stamp = 0; std::vector<int> face; std::vector<char> allowed; };
                            std::vector<PlacementFilter> pfilter(local_edges.size());
                            std::vector<std::uint64_t> placed_stamp(local_edges.size(), 0);
                            std::uint64_t stamp_counter = 0;
                            if (config_.matching_edges_first) {
                                // edges to the old drawing (here: the matching) first, in their original
                                // relative order: they place the new vertices directly near their partners
                                std::stable_partition(local_edges.begin(), local_edges.end(),
                                        [nm](const Edge& e) { return e[0] < nm || e[1] < nm; });
                            }
                            auto e = local_edges.begin();
                            uint64_t total_local_iterations = 0;

                            // does path (for the edge at position pos, target vertex target) land in an allowed face?
                            auto landing_ok = [&](std::size_t pos, const HdsPath& path, std::size_t target) -> bool {
                                if (!config_.placement_filter || d.vertices[target].halfedge) return true;
                                PlacementFilter& f = pfilter[pos];
                                const std::uint64_t base = pos == 0 ? 0 : placed_stamp[pos - 1];
                                if (!f.valid || f.stamp != base) {
                                    std::size_t num_faces = 0;
                                    std::map<std::size_t, std::vector<char>> domains;
                                    const bool possible = placements_still_possible(d, local_edges, pos + 1, constrained, &f.face, &num_faces, &domains);
                                    auto it = domains.find(target);
                                    if (!possible) f.allowed.assign(num_faces, 0);
                                    else if (it == domains.end()) f.allowed.assign(num_faces, 1);
                                    else f.allowed = it->second;
                                    f.valid = true;
                                    f.stamp = base;
                                }
                                // the new vertex lands in the face of p[0] (no crossing) or behind the last crossing
                                const HdsHalfedge* h = path.size() == 2 ? path[0] : path[path.size() - 2]->twin;
                                if (f.allowed[f.face[h->label]]) return true;
                                ++result.filtered_paths_count;
                                return false;
                            };
                            auto next_path_filtered = [&](HdsPath& path, std::size_t target, int prescribed, std::size_t pos) {
                                while (d.next_path(path, target, prescribed))
                                    if (landing_ok(pos, path, target)) return true;
                                return false;
                            };

                            for (;;) {
                                ++total_local_iterations;
                                if (config_.progress_seconds > 0 && (total_local_iterations & 0xffff) == 0) {
                                    auto now = std::chrono::steady_clock::now();
                                    if (std::chrono::duration<double>(now - last_progress).count() >= config_.progress_seconds) {
                                        last_progress = now;
                                        std::cout << "[PROGRESS] t=" << std::chrono::duration<double>(now - run_start).count()
                                            << "s parent " << solcount << " pass " << (constrained + 1)
                                            << " dfs " << total_local_iterations
                                            << " leaves " << result.leaf_count
                                            << " children " << (config_.split_mode ? result.children_exported : result.solutions.size())
                                            << " dup " << result.discarded_count
                                            << " pruned " << result.pruned_early_count
                                            << " (placement " << result.pruned_placement_count
                                            << ", final face " << result.pruned_final_face_count
                                            << ", mcf " << result.pruned_mcf_count << ")"
                                            << " units " << result.split_units
                                            << " done " << done_units.size() + result.units_done << std::endl;
                                    }
                                }
                                std::size_t u = (*e)[0];
                                std::size_t v = (*e)[1];

                                // In Pass 2 (constrained == 1), newly added C12 cycle edges are uncrossable (pcr = klim). 
                                // In Pass 1 (constrained == 0), edges are unconstrained (pcr = 0).
                                int pcr = (constrained == 1 && e->size() == 3) ? (*e)[2] : 0;
                                HdsPath p = d.first_path(u, v, pcr);
                                {
                                    const std::size_t pos = static_cast<std::size_t>(e - local_edges.begin());
                                    if (!p.empty() && !landing_ok(pos, p, v) && !next_path_filtered(p, v, pcr, pos)) p.clear();
                                }

                                if (p.empty()) {
BACKUP:
                                    do {
                                        if (e == local_edges.begin()) {
                                            goto FINISH_PASS;
                                        }
                                        --e;
                                        u = (*e)[0];
                                        v = (*e)[1];
                                        pcr = (constrained == 1 && e->size() == 3) ? (*e)[2] : 0;
                                        p = d.edges.back().built;
                                        d.remove_edge();

                                        if (next_path_filtered(p, v, pcr, static_cast<std::size_t>(e - local_edges.begin()))) {
                                            break;
                                        }
                                    } while (true);
                                }
                                d.add_edge(p, v, pcr);

PLACED: // edge e has just been added (re-entered when the split check rejects a path)
                                // reachability check for remaining unplaced braid edges
                                std::size_t current_edge_idx = static_cast<std::size_t>(e - local_edges.begin());
                                placed_stamp[current_edge_idx] = ++stamp_counter;
                                if (config_.enable_early_pruning && (current_edge_idx + 1 < local_edges.size())) {
                                    bool is_checkpoint = false;
                                    is_checkpoint = (
                                        std::find(config_.early_prune_checkpoints.begin(),config_.early_prune_checkpoints.end(),current_edge_idx) 
                                        != config_.early_prune_checkpoints.end());

                                    // optional: the cheap reachability check after every edge from reach_every_edge_from on
                                    bool every_edge = current_edge_idx >= config_.reach_every_edge_from;
                                    bool final_face_check = constrained == 0 && current_edge_idx >= config_.final_face_every_edge_from;
                                    bool placement_check = current_edge_idx >= config_.placement_every_edge_from;

                                    if (is_checkpoint || every_edge || final_face_check || placement_check) {
                                        // joint routing check, only run if every edge is individually reachable
                                        auto mcf_prunes = [&]() -> bool {
                                            if (!config_.enable_mcf_pruning || !is_checkpoint) return false;
                                            auto t0 = std::chrono::steady_clock::now();
                                            McfResult r = check_remaining_edges_mcf(d, local_edges, current_edge_idx + 1, constrained);
                                            result.mcf_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                                            if (r == McfResult::Unknown) ++result.mcf_unknown_count;
                                            if (r == McfResult::Infeasible) ++result.pruned_mcf_count;

                                            // throttled progress log (at most every 0.5s)
                                            if (config_.mcf_log_progress && t0 - last_mcf_log >= std::chrono::milliseconds(500)) {
                                                last_mcf_log = t0;
                                                std::cout << "  [MCF] drawing #" << solcount
                                                    << " pass " << (constrained + 1)
                                                    << " checkpoint " << current_edge_idx
                                                    << ": pruned " << result.pruned_mcf_count
                                                    << " / early pruned " << result.pruned_early_count
                                                    << ", budget exhausted " << result.mcf_unknown_count
                                                    << ", mcf time " << result.mcf_seconds << "s" << std::endl;
                                            }
                                            return r == McfResult::Infeasible;
                                        };

                                        auto placement_prunes = [&]() -> bool {
                                            if (!placement_check || placements_still_possible(d, local_edges, current_edge_idx + 1, constrained)) return false;
                                            ++result.pruned_placement_count;
                                            return true;
                                        };

                                        auto final_face_prunes = [&]() -> bool {
                                            if (!final_face_check || final_face_still_possible(d, nm)) return false;
                                            ++result.pruned_final_face_count;
                                            return true;
                                        };

                                        while (!check_remaining_edges_reachability(d, local_edges, current_edge_idx + 1) || placement_prunes() || mcf_prunes() || final_face_prunes()) {
                                            // std::string filename2 = "../quasiDrawings/failExample.graphml";
                                            // std::ofstream of_graphml(filename2);
                                            // d.graphml_output(of_graphml);
                                            // of_graphml.close();
                                            // if (config_.verbose) {
                                            //     std::cout << "  [Early Prune] Edge index " << current_edge_idx 
                                            //         << ": Remaining edges cannot be routed. Early prune count: " << result.pruned_early_count << std::endl;
                                            // }
                                            ++result.pruned_early_count;
                                            p = d.edges.back().built; // edge e fails, remove and try alt paths
                                            d.remove_edge();

                                            // keep trying next path 
                                            if (next_path_filtered(p, v, pcr, current_edge_idx)) {
                                                d.add_edge(p, v, pcr);
                                                placed_stamp[current_edge_idx] = ++stamp_counter;
                                            }
                                            // until no alternatie paths
                                            else goto BACKUP;
                                        }
                                    }
                                }

                                // split_mode: only search below the prefixes owned by this process
                                if (config_.split_mode && config_.split_depth > 0 && current_edge_idx + 1 == config_.split_depth) {
                                    ++result.split_units;
                                    close_unit(); // reaching the next unit completes the previous one
                                    if (owns(unit)) {
                                        open_unit = unit++;
                                    } else {
                                        ++unit;
                                        p = d.edges.back().built;
                                        d.remove_edge();
                                        if (next_path_filtered(p, v, pcr, current_edge_idx)) {
                                            d.add_edge(p, v, pcr);
                                            goto PLACED;
                                        }
                                        goto BACKUP;
                                    }
                                }

                                if (++e == local_edges.end()) {
                                    if (constrained == 1)  {
                                        // =========================================================
                                        // PASS 2: CONSTRAINED COMPLETE DRAWING FOUND
                                        // =========================================================
                                        if (config_.verbose) {
                                            std::cout << ">>> [Pass 2] SUCCESS: Found Complete Drawing #" << result.full_solution_count 
                                                << " (from parent #" << solcount << ") <<<" << std::endl;
                                        }


                                        if (config_.split_mode) {
                                            // one file for all complete drawings of this task, one JSON per line
                                            std::ofstream of_json(config_.output_dir + "/solutions.jsonl", std::ios::app);
                                            nlohmann::ordered_json j = d.serialize_to_json();
                                            j["meta"] = gadget_crossings(d, nm);
                                            j["meta"]["parent"] = describe_parent(parent_refs[solcount]);
                                            of_json << j.dump() << "\n";
                                        } else if (config_.export_files) {
                                            std::string base_path = config_.output_dir + "/solution_" + std::to_string(result.full_solution_count);

                                            std::ofstream of_gml(base_path + ".graphml");
                                            d.graphml_output(of_gml);

                                            std::ofstream of_json(base_path + ".json");
                                            nlohmann::ordered_json output_json = d.serialize_to_json();
                                            of_json << output_json.dump(4);
                                        }

                                        ++result.full_solution_count;
                                        goto BACKUP;
                                    }

                                    ++result.leaf_count;
                                    // map halfedges to faces
                                    std::vector<int> face;
                                    std::size_t num_faces = label_faces(d, face);
                                    if (num_faces == 0) goto BACKUP;

                                    // Flow network vertex indexing:
                                    //  - Faces: 0 ... num_faces - 1
                                    //  - Active cycle vertices: num_faces ... num_faces + config_.cycle_size - 1
                                    //  - Source node: num_faces + config_.cycle_size
                                    std::size_t source = num_faces + config_.cycle_size;
                                    ThresholdFlow net(source + 1);
                                    // debug: the Boost push-relabel version of the same network, see Config::verify_threshold_flow
                                    std::unique_ptr<DualNetwork> check_net;
                                    if (config_.verify_threshold_flow) check_net = std::make_unique<DualNetwork>(source + 1);
                                    auto add_net_edge = [&](std::size_t from, std::size_t to, int capacity) {
                                        net.add_edge(static_cast<int>(from), static_cast<int>(to), capacity);
                                        if (check_net) check_net->add_edge(static_cast<int>(from), static_cast<int>(to), capacity);
                                    };
                                    // is the max flow from s to t at least limit?
                                    auto flow_at_least = [&](std::size_t s, std::size_t t, int limit) {
                                        bool r = net.flow(static_cast<int>(s), static_cast<int>(t), limit) >= limit;
                                        if (check_net && r != (check_net->flow(static_cast<int>(s), static_cast<int>(t)) >= limit))
                                            throw std::runtime_error("ThresholdFlow disagrees with push_relabel_max_flow");
                                        return r;
                                    };
                                    std::vector<std::vector<std::size_t>> dualg(num_faces);

                                    // source node to active cycle vertex
                                    for (std::size_t i = 0; i < config_.cycle_size; ++i)
                                        add_net_edge(source, num_faces + i, 1);

                                    // edges between faces and from (active cycle) vertices to incident faces
                                    std::size_t active_cycle_end = nm + config_.cycle_size;
                                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                                        // active cycle vertex -> incident face
                                        if (i->vertex->label >= nm && i->vertex->label < active_cycle_end)
                                            add_net_edge(num_faces + (i->vertex->label - nm), face[i->label], 1);

                                        // face -> adjacent face dual edge
                                        int remaining_capacity = static_cast<int>(klim) - static_cast<int>(i->edge->ncr);
                                        if (remaining_capacity > 0) {
                                            add_net_edge(face[i->label], face[i->twin->label], remaining_capacity);
                                            dualg[face[i->label]].push_back(face[i->twin->label]);
                                        }
                                    }

                                    // Potential final faces
                                    std::vector<int> pff;
                                    for (std::size_t f = 0; f < num_faces; ++f) {
                                        if (flow_at_least(source, f, static_cast<int>(config_.cycle_size))) {
                                            pff.push_back(static_cast<int>(f));
                                        }
                                    }
                                    if (pff.empty()) goto BACKUP;

                                    // Face classification:
                                    //   active = 3 -> possible final face (PFF)
                                    //   active = 2 -> active face
                                    //   active = 1 -> passive face
                                    //   active = 0 -> transit face
                                    //   active = -1 -> irrelevant face
                                    std::vector<int> active(num_faces, -1);
                                    for (int target_f : pff) active[target_f] = 3;

                                    // Active faces: flow >= 3 to a PFF
                                    for (std::size_t i = 0; i < num_faces; ++i) {
                                        if (active[i] != -1) continue;
                                        for (int target_f : pff) {
                                            if (flow_at_least(i, target_f, 3)) {
                                                active[i] = 2; break;
                                            }
                                        }
                                    }

                                    // Compute map: face -> representative boundary halfedge
                                    std::vector<const HdsHalfedge*> fhedge(num_faces, nullptr);
                                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                                        if (fhedge[face[i->label]] == nullptr) fhedge[face[i->label]] = &*i;
                                    }

                                    // dual distances (over crossable edges, i.e. dualg) from a set of faces
                                    auto dual_distance = [&](auto in_set) {
                                        std::vector<int> dist(num_faces, -1);
                                        std::vector<std::size_t> queue;
                                        for (std::size_t f = 0; f < num_faces; ++f) if (in_set(f)) { dist[f] = 0; queue.push_back(f); }
                                        for (std::size_t qi = 0; qi < queue.size(); ++qi)
                                            for (std::size_t g : dualg[queue[qi]])
                                                if (dist[g] < 0) { dist[g] = dist[queue[qi]] + 1; queue.push_back(g); }
                                        return dist;
                                    };
                                    const int kmax = static_cast<int>(klim); // max #crossings per edge (3)
                                    const std::vector<int> dist_active = dual_distance([&](std::size_t f) { return active[f] >= 2; });

                                    // Passive faces: not active, incident to an active cycle vertex v, and a dual path of length
                                    // <= klim (= 3) to an active face/PFF whose first edge is not incident to v: an edge from v to
                                    // the next cycle has at most klim crossings and ends in an active face (the 2-planar paper
                                    // uses length <= 2)
                                    for (std::size_t i = 0; i < num_faces; ++i) {
                                        if (fhedge[i] == nullptr) throw std::runtime_error("no edge for face");
                                        if (active[i] != -1) continue;
                                        const HdsHalfedge* e_curr = fhedge[i];
                                        do {
                                            std::size_t v = e_curr->vertex->label;
                                            if (v >= nm && v < active_cycle_end) // active vertex
                                                for (const HdsHalfedge* f = e_curr->next->next; f != e_curr; f = f->next) {
                                                    if (f->edge->ncr < klim &&
                                                            f->vertex->label != v &&
                                                            f->twin->vertex->label != v) {
                                                        // first edge crosses f, the rest (<= klim-1 steps) leads to an active face
                                                        int dfn = dist_active[face[f->twin->label]];
                                                        if (dfn >= 0 && dfn <= kmax - 1) { active[i] = 1; break; }
                                                    }
                                                }
                                            e_curr = e_curr->next;
                                        } while (active[i] == -1 && e_curr != fhedge[i]);
                                    }

                                    // Transit faces: an edge with <= klim (= 3) crossings that starts in an active or passive face
                                    // and ends in an active face passes through up to klim-1 (= 2) faces in between, so a face is
                                    // transit if (dual distance from an active/passive face) + (dual distance to an active face)
                                    // <= klim (the 2-planar paper: adjacent to an active and an active/passive face)
                                    const std::vector<int> dist_ap = dual_distance([&](std::size_t f) { return active[f] >= 1; });
                                    for (std::size_t i = 0; i < num_faces; ++i)
                                        if (active[i] == -1 && dist_ap[i] >= 1 && dist_active[i] >= 1 && dist_ap[i] + dist_active[i] <= kmax)
                                            active[i] = 0;

                                    // determine relevant vertices
                                    // ensure that the active cycle vertices offset, ..., offset+12-1, are mapped to 0...13 in relv
                                    std::vector<const HdsHalfedge*> relv(config_.cycle_size, nullptr); // relevant halfedges
                                    std::vector<const HdsHalfedge*> prelv; // possibly relevant halfedges

                                    for (auto i = d.vertices.begin(); i != d.vertices.end(); ++i) {
                                        auto e = i->halfedge;
                                        std::size_t nirf = 0; // #incident relevant faces
                                        auto a = e;
                                        do {
                                            if (active[face[e->label]] >= 0) { ++nirf; a = e; }
                                            e = e->next->twin;
                                        } while (e != i->halfedge);

                                        if (i->label >= nm && i->label < active_cycle_end) {
                                            // if no incident face is relevant, then 
                                            // this isnt a valid solution
                                            if (nirf == 0) goto BACKUP;
                                            relv[(i->label) - nm] = a;
                                        } else if (nirf >= 2)
                                            relv.push_back(a);
                                        else if (nirf == 1)
                                            prelv.push_back(a);
                                    }

                                    // ... and crossings
                                    for (auto i = d.crossings.begin(); i != d.crossings.end(); ++i) {
                                        auto e = i->halfedge;
                                        std::size_t nirf = 0;
                                        auto a = e;
                                        do {
                                            if (active[face[e->label]] >= 0) { ++nirf; a = e; }
                                            e = e->next->twin;
                                        } while (e != i->halfedge);

                                        if (nirf >= 2)
                                            relv.push_back(a);
                                        else if (nirf == 1)
                                            prelv.push_back(a);
                                    }

                                    // make sure triangles are not contracted to (non-simple) lenses
                                    std::vector<std::size_t> fsize(num_faces, 0); // size of relevant faces
                                    std::size_t total_nv = d.vertices.size() + d.crossings.size();
                                    std::vector<std::size_t> vact(total_nv, total_nv); // vertex mapping: old label -> new label
                                    std::size_t vic = 0;

                                    for (auto i = relv.begin(); i != relv.end(); ++i) {
                                        auto j = *i;
                                        do {
                                            ++fsize[face[j->label]];
                                            j = j->next->twin;
                                        } while (j != *i);
                                        vact[(*i)->vertex->label] = vic++;
                                    }

                                    for (auto i = prelv.begin(); i != prelv.end(); ++i)
                                        if (fsize[face[(*i)->label]] < 3) { // if lense, add prelv halfedge
                                            ++fsize[face[(*i)->label]];
                                            relv.push_back(*i);
                                            vact[(*i)->vertex->label] = vic++;
                                        }

                                    Drawing<klim> nd(relv.size() + 1); // extract relevant parts into pruned drawing
                                    std::deque<const HdsHalfedge*> todo(1, relv[0]);
                                    // vertex status: -1 == not considered yet, 0 == in queue and in nd, 
                                    // 1 == handled and incident edges in nd
                                    std::vector<int> done(d.vertices.size() + d.crossings.size(), -1);
                                    // We merge connected regions of irrelevant faces into "black" regions. Record them here. 
                                    std::vector<HdsHalfedge*> black;

                                    while (!todo.empty()) {
                                        auto i = todo.front();
                                        todo.pop_front();
                                        std::size_t u = i->vertex->label;

                                        if (active[face[i->label]] < 0)
                                            throw std::runtime_error("irrelevant face");

                                        HdsHalfedge* last = nd.vertices[vact[u]].halfedge;
                                        if (nd.edges.empty()) {
                                            // first edge
                                            auto p = prev_active<klim>(i, face, active, vact);
                                            std::size_t v = p.first->vertex->label;
                                            last = nd.add_first_edge(vact[u], vact[v], p.second)->twin;
                                            if (last->vertex->label != 0)
                                                throw std::runtime_error("wrong first edge in new drawing");
                                            todo.push_back(p.first);
                                            done[v] = 0;
                                        } else if (last == nullptr) {
                                            todo.push_back(i);
                                            continue;
                                        }

                                        // act[u] is now connected in nd and last points to act[u]; 
                                        // we build neighborhood of act[u] in nd
                                        done[u] = 1;
                                        std::size_t wnd = last->twin->vertex->label; // neighbor of u in nd
                                        std::size_t w = 0; // label of wnd in d
                                        auto wi = i;
                                        do {
                                            auto e_prev = prev_active<klim>(wi, face, active, vact).first;
                                            w = e_prev->vertex->label;
                                            if (vact[w] == wnd) break;
                                            wi = wi->next->twin;
                                        } while (wi != i);
                                        if (vact[w] != wnd) throw std::runtime_error("no w");

                                        // edges incident to u in d
                                        for (auto j = wi->next->twin; j != wi; j = j->next->twin) {
                                            if (active[face[j->label]] < 0 && active[face[j->twin->label]] < 0)
                                                continue;

                                            auto p = prev_active<klim>(j, face, active, vact);
                                            std::size_t v = p.first->vertex->label;
                                            if (last->next->vertex->label == vact[v]) {
                                                // edge uv already present
                                                last = last->next->twin;
                                                continue;
                                            }

                                            HdsPath path(2, last);
                                            if (!nd.find_target(path, vact[v]))
                                                throw std::runtime_error("cannot find v");

                                            auto ne = nd.add_edge(path, vact[v], p.second);
                                            if (active[face[j->label]] < 0)
                                                black.push_back(ne->twin);
                                            else if (active[face[j->twin->label]] < 0)
                                                black.push_back(ne);

                                            last = last->next->twin;
                                            if (done[v] < 0) {
                                                done[v] = 0;
                                                todo.push_back(p.first);
                                            }
                                        }
                                    }

                                    // process black faces
                                    std::vector<bool> blackdone(nd.halfedges.size(), false);
                                    std::vector<HdsHalfedge*> large_black_faces;

                                    for (auto i = black.begin(); i != black.end(); ++i) {
                                        if (blackdone[(*i)->label]) continue;
                                        auto j = *i;
                                        std::size_t bc = 0;
                                        do {
                                            blackdone[j->label] = true;
                                            ++bc;
                                            j = j->next;
                                        } while (j != *i);

                                        if (bc >= 4) large_black_faces.push_back(*i);
                                    }
                                    // allocate extra star vertices if more than one black face
                                    if (large_black_faces.size() > 1)
                                        nd.add_vertices(large_black_faces.size() - 1);

                                    // add an uncrossable star in each (large) black face
                                    std::size_t star_v = relv.size();
                                    for (auto start_edge : large_black_faces) {
                                        auto j = start_edge;
                                        // add spoke to first boundary vertex
                                        auto x = nd.add_edge(HdsPath({j, nullptr}), star_v, klim);
                                        for (;;) {
                                            j = j->next->twin->next; // jump over spoke
                                            if (j == start_edge) break;
                                            nd.add_edge(HdsPath({j, x}), star_v, klim);
                                        }
                                        ++star_v;
                                    }


                                    if (config_.split_mode && reported_extensible.insert(solcount).second) {
                                        // written right away, so it survives a job killed by its time limit
                                        std::ofstream ext(config_.output_dir + "/extensible.txt", std::ios::app);
                                        ext << describe_parent(parent_refs[solcount]) << "\n";
                                    }
                                    auto t_iso = std::chrono::steady_clock::now();
                                    IsoCode code = canonical_code(nd, syms);
                                    auto known = seen_codes.find(code);

                                    if (config_.verify_iso) {
                                        bool pairwise = false;
                                        for (const auto& x : result.solutions)
                                            if (are_isomorphic_gadget(x.drawing, nd, syms)) { pairwise = true; break; }
                                        if (pairwise != (known != seen_codes.end()))
                                            throw std::runtime_error("canonical code disagrees with are_isomorphic_gadget");
                                        if (!pairwise)
                                            for (const auto& x : result.solutions)
                                                if (are_isomorphic(x.drawing, nd)) { ++result.label_agnostic_only_count; break; }
                                    }
                                    result.iso_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_iso).count();

                                    if (known != seen_codes.end()) {
                                        if (config_.verbose) {
                                            std::cout << "--- Discard drawing #" << result.discarded_count
                                                << ", isomorphic to solution #" << known->second << std::endl;
                                        }
                                        ++result.discarded_count;
                                        goto BACKUP;
                                    }

                                    // we have a new, valid solution -> record it
                                    if (!nd.is_valid()) throw std::runtime_error("nd is invalid");
                                    if (!config_.split_mode || config_.split_keep_children) {
                                        seen_codes.emplace(std::move(code), result.solutions.size());
                                        result.solutions.push_back({nd, current_depth+1});
                                    } else {
                                        if (seen_codes.size() >= config_.split_dedup_limit) seen_codes.clear();
                                        seen_codes.emplace(std::move(code), num_inputs + result.children_exported);
                                    }
                                    is_extensible = true; // Mark as extensible to enable Pass 2
                                    if (config_.split_mode) {
                                        if (config_.split_export_children) {
                                            std::ofstream of_json(config_.output_dir + "/children.jsonl", std::ios::app);
                                            nlohmann::ordered_json j = nd.serialize_to_json();
                                            j["meta"] = gadget_crossings(d, nm);
                                            j["meta"]["parent"] = describe_parent(parent_refs[solcount]);
                                            of_json << j.dump() << "\n";
                                        }
                                        ++result.children_exported;
                                    } else {
                                        std::cout << "Drawing #" << result.solutions.size()-1 << ":\n"
                                            << d << std::endl;
                                    }

                                    if (pff.size() > 1 && (!config_.split_mode || config_.verbose))
                                        std::cout << "!!! Final face is not unique for this drawing ---"
                                            << std::endl;

                                    if (config_.verbose) {
                                        std::cout << "[Pass 1] Added Intermediate Drawing #" << result.solutions.size() - 1 
                                            << " to Queue." << std::endl;
                                    }

                                    goto BACKUP; 
                                } // if (++e == local_edges.end())
                            } // for (auto e = edges.begin();;)

FINISH_PASS:
                            close_unit(); // the last owned unit of this pass is complete
                            std::cout << "Finished Pass " << (constrained == 0 ? "1" : "2") 
                                << " for Drawing #" << solcount 
                                << " (Total DFS Iterations: " << total_local_iterations
                                << ", MCF pruned so far: " << result.pruned_mcf_count << ")" << std::endl;
                        } // for (int constrained = 0; constrained < 2; constrained++)
                        ++solcount;
                    } // while (solcount < solutions.size())

                    result.total_processed = solcount;
                    return result;
                }

            private:
                Config config_;
        };
} // namespace hds

#endif // NESTED_CYCLE_SEARCH_HPP


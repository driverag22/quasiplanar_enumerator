// Shared helpers for the validation drivers.
// Include exactly one hds_*.h header and iso.h before this file.
#ifndef VALIDATION_COMMON_H
#define VALIDATION_COMMON_H
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

inline double now_s() {
    static const auto t0 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// peak resident memory of this process in MB (VmHWM)
inline double peak_mem_mb() {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("VmHWM:", 0) == 0) return std::stod(line.substr(6)) / 1024.0;
    return -1;
}

// ---------------------------------------------------------------------------
// Canonical codes of the planarization.
//
// From a start halfedge h pointing to a real vertex, a BFS numbers the nodes (vertices and crossings) in order of
// discovery and walks the rotation of every node, starting at the halfedge through which it was discovered. The
// code lists, for every node, its neighbours in rotation order (number and type real/crossing). Two drawings
// have a common code iff their planarizations are isomorphic by a map that sends vertices to vertices, crossings to
// crossings and preserves rotations (mirror = false), i.e. iff they are orientation-preservingly isomorphic up to a
// relabelling of the vertices (which is then automatically an automorphism of the graph). With mirror = true all
// rotations are read reversed, which gives the codes of the mirror image.
// The canonical code is the lexicographic minimum over all start halfedges.
// PRE: the drawn part is connected (isolated vertices are ignored).
// ---------------------------------------------------------------------------
typedef std::vector<uint32_t> Code;

template <class D>
struct Canonizer {
    std::vector<uint32_t> num;
    std::vector<const HdsHalfedge*> queue;
    Code cur;

    // code from start halfedge h; aborts early (returns false) as soon as the code is known to be > best
    // (best may be empty = +infinity). Returns true if cur < best (cur is then complete).
    bool code_from(const D& d, const HdsHalfedge* h, bool mirror, const Code& best) {
        const std::size_t np = d.vertices.size();
        const std::size_t nn = np + d.crossings.size();
        num.assign(nn, UINT32_MAX);
        queue.clear();
        cur.clear();
        int cmp = best.empty() ? -1 : 0;   // -1: already smaller, 0: equal so far
        auto emit = [&](uint32_t x) -> bool {
            if (cmp == 0) {
                const std::size_t k = cur.size();
                if (k >= best.size()) return false;   // longer with equal prefix: greater
                if (x < best[k]) cmp = -1;
                else if (x > best[k]) return false;
            }
            cur.push_back(x);
            return true;
        };
        uint32_t next = 0;
        num[h->vertex->label] = next++;
        queue.push_back(h);
        for (std::size_t qi = 0; qi < queue.size(); ++qi) {
            const HdsHalfedge* s = queue[qi];
            const HdsHalfedge* t = s;
            do {
                const std::size_t o = t->twin->vertex->label;
                if (num[o] == UINT32_MAX) {
                    num[o] = next++;
                    queue.push_back(t->twin);
                }
                if (!emit(2 * num[o] + (o < np ? 2 : 1))) return false;
                t = mirror ? t->twin->prev : t->next->twin;
            } while (t != s);
            if (!emit(0)) return false;   // end of rotation
        }
        if (cmp == 0) return cur.size() < best.size();
        return true;
    }

    // canonical code over all starts at real vertices, for one orientation
    Code canon(const D& d, bool mirror) {
        Code best;
        for (const auto& v : d.vertices) {
            if (v.halfedge == nullptr) continue;
            const HdsHalfedge* h = v.halfedge;
            do {
                if (code_from(d, h, mirror, best)) best = cur;
                h = h->next->twin;
            } while (h != v.halfedge);
        }
        return best;
    }

    // label-preserving, orientation-preserving code: start at the lowest drawn vertex, at its edge to the
    // lowest neighbour; real vertices are written with their labels
    Code labelled(const D& d) {
        const std::size_t np = d.vertices.size();
        const HdsHalfedge* start = nullptr;
        std::size_t w0 = SIZE_MAX;
        std::size_t v0 = SIZE_MAX;
        for (std::size_t v = 0; v < np && start == nullptr; ++v) {
            if (d.vertices[v].halfedge == nullptr) continue;
            v0 = v;
            const HdsHalfedge* h = d.vertices[v].halfedge;
            do {
                const std::size_t w = (h->edge->u == v ? h->edge->v : h->edge->u);
                if (w < w0) { w0 = w; start = h; }
                h = h->next->twin;
            } while (h != d.vertices[v].halfedge);
        }
        Code res;
        if (start == nullptr) return res;
        res.push_back(v0);
        const std::size_t nn = np + d.crossings.size();
        num.assign(nn, UINT32_MAX);
        queue.clear();
        uint32_t next = 0;
        num[start->vertex->label] = next++;
        queue.push_back(start);
        for (std::size_t qi = 0; qi < queue.size(); ++qi) {
            const HdsHalfedge* s = queue[qi];
            const HdsHalfedge* t = s;
            do {
                const std::size_t o = t->twin->vertex->label;
                if (num[o] == UINT32_MAX) {
                    num[o] = next++;
                    queue.push_back(t->twin);
                }
                res.push_back(o < np ? 2 * o + 2 : 2 * num[o] + 1);
                t = t->next->twin;
            } while (t != s);
            res.push_back(0);
        }
        return res;
    }
};

inline uint64_t hash_code(const Code& c) {
    uint64_t h = 1469598103934665603ull;
    for (uint32_t x : c) {
        h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        h *= 1099511628211ull;
    }
    return h;
}
// 128-bit fingerprint of a code (two independent 64-bit hashes); used as key to keep memory small
struct Key {
    uint64_t a, b;
    bool operator==(const Key& o) const { return a == o.a && b == o.b; }
    bool operator<(const Key& o) const { return a < o.a || (a == o.a && b < o.b); }
};
struct KeyHash {
    std::size_t operator()(const Key& k) const { return k.a ^ (k.b * 0x9e3779b97f4a7c15ull); }
};
inline Key code_str(const Code& c) {
    uint64_t h2 = 0x84222325cbf29ce4ull;
    for (uint32_t x : c) {
        h2 += x;
        h2 *= 0xff51afd7ed558ccdull;
        h2 ^= h2 >> 33;
    }
    return Key{hash_code(c), h2 ^ c.size()};
}

// cheap invariant: number of crossings and sorted crossing numbers of the edges
template <class D>
std::string invariant(const D& d) {
    std::vector<std::size_t> ncr;
    for (const auto& e : d.edges) ncr.push_back(e.ncr);
    std::sort(ncr.begin(), ncr.end());
    std::string s = std::to_string(d.crossings.size()) + ":";
    for (auto x : ncr) s += std::to_string(x) + ",";
    return s;
}

// ---------------------------------------------------------------------------
// Exhaustive drawing of edges[i..j-1] in all ways (no pruning), as in tools/enumeration.cpp / k10.cpp.
// edges[k][0] must be drawn when edge k is drawn.
// ---------------------------------------------------------------------------
template <class D>
void draw_all(D& d, const Edges& edges, std::size_t i, std::size_t j, uint64_t& nodes,
              const std::function<void(D&)>& done) {
    if (i == j) {
        done(d);
        return;
    }
    const std::size_t a = edges[i][0], b = edges[i][1];
    for (HdsPath p = d.first_path(a, b); !p.empty();) {
        d.add_edge(p, b);
        ++nodes;
        draw_all(d, edges, i + 1, j, nodes, done);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, b)) break;
    }
}

// ---------------------------------------------------------------------------
// Class bookkeeping with optional cross-checks against iso.h
// ---------------------------------------------------------------------------
template <class D>
struct ClassTable {
    struct Info {
        Key oriented;           // oriented canonical code of the representative
        bool achiral;
        std::size_t rep;        // index into reps (if stored)
        uint64_t count = 0;     // labelled drawings in this class
    };
    bool check_iso = false;
    std::unordered_map<Key, Info, KeyHash> classes;   // strong canonical code -> info
    std::unordered_set<Key, KeyHash> oriented;          // oriented canonical codes seen
    std::vector<D> reps;
    std::vector<uint64_t> labelled_hashes;
    bool keep_hashes = false;
    uint64_t labelled = 0;
    uint64_t iso_mismatch = 0;   // disagreements between iso.h and the canonical codes
    uint64_t iso_calls = 0;
    Canonizer<D> cz;

    void add(const D& d) {
        ++labelled;
        if (keep_hashes) labelled_hashes.push_back(hash_code(cz.labelled(d)));
        Code co = cz.canon(d, false);
        Code cm = cz.canon(d, true);
        const bool achiral = (co == cm);
        const Key so = code_str(co);
        const Key ss = code_str(std::min(co, cm));
        oriented.insert(so);
        auto it = classes.find(ss);
        if (it == classes.end()) {
            Info inf{so, achiral, reps.size(), 1};
            if (check_iso) {
                reps.push_back(d);
                ++iso_calls;
                if (are_isomorphic_oriented(d, d, true) != achiral) ++iso_mismatch;
            }
            classes.emplace(ss, inf);
        } else {
            ++it->second.count;
            if (check_iso) {
                const D& r = reps[it->second.rep];
                iso_calls += 2;
                if (!are_isomorphic(r, d)) ++iso_mismatch;
                if (are_isomorphic_oriented(r, d, false) != (so == it->second.oriented)) ++iso_mismatch;
            }
        }
    }
};

#endif

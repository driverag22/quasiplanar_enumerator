// Canonical codes of drawings, used by enum.cpp to count drawings up to strong isomorphism.
//
// A drawing is stored as the planarization (vertices and crossings, with rotations). From a start halfedge h at a
// vertex, a BFS numbers the nodes in order of discovery and walks the rotation of every node, starting at the halfedge
// through which it was discovered. The code lists, for every node, its neighbours in rotation order (BFS number and
// type, vertex or crossing). Two connected drawings have the same code for some pair of start halfedges iff their
// planarizations are isomorphic by a map that sends vertices to vertices and crossings to crossings and preserves all
// rotations. With mirror = true, all rotations are read in reverse, which gives the codes of the mirror image.
// The canonical code is the lexicographic minimum over all start halfedges, and the strong code is the minimum of the
// canonical codes of the drawing and of its mirror image. So two drawings have the same strong code iff they are
// strongly isomorphic; vertex labels are not used, so the relabelling is automatically an automorphism of the graph.
// Include it after hds_kplanar.h or hds_quasiplanar.h.
#ifndef VALIDATION_CANON_H
#define VALIDATION_CANON_H
#include <cstdint>
#include <vector>

typedef std::vector<uint32_t> Code;

template <class D>
struct Canonizer {
    std::vector<uint32_t> num;               // BFS number of every node (vertices first, then crossings)
    std::vector<const HdsHalfedge*> queue;
    Code cur;

    // the code from start halfedge h; stops early and returns false as soon as it is known to be >= best
    // (an empty best means +infinity); returns true if cur < best, and then cur is complete
    bool code_from(const D& d, const HdsHalfedge* h, bool mirror, const Code& best) {
        const std::size_t nv = d.vertices.size();
        num.assign(nv + d.crossings.size(), UINT32_MAX);
        queue.clear();
        cur.clear();
        bool smaller = best.empty();   // cur is already smaller than best (otherwise: equal so far)
        auto emit = [&](uint32_t x) {
            if (!smaller) {
                const std::size_t k = cur.size();
                if (k >= best.size() || x > best[k]) return false;
                if (x < best[k]) smaller = true;
            }
            cur.push_back(x);
            return true;
        };
        uint32_t next = 0;
        num[h->vertex->label] = next++;
        queue.push_back(h);
        for (std::size_t qi = 0; qi < queue.size(); ++qi) {
            const HdsHalfedge* t = queue[qi];
            do {
                const std::size_t o = t->twin->vertex->label;
                if (num[o] == UINT32_MAX) {
                    num[o] = next++;
                    queue.push_back(t->twin);
                }
                if (!emit(2 * num[o] + (o < nv ? 2 : 1))) return false;   // even: vertex, odd: crossing
                t = mirror ? t->twin->prev : t->next->twin;
            } while (t != queue[qi]);
            if (!emit(0)) return false;   // end of the rotation
        }
        return smaller || cur.size() < best.size();
    }

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

    Code strong(const D& d) { return std::min(canon(d, false), canon(d, true)); }
};

// 128-bit fingerprint of a code (two independent 64-bit hashes), to keep memory small for millions of drawings
struct Key {
    uint64_t a, b;
    bool operator==(const Key& o) const { return a == o.a && b == o.b; }
};
struct KeyHash {
    std::size_t operator()(const Key& k) const { return k.a ^ (k.b * 0x9e3779b97f4a7c15ull); }
};
inline Key fingerprint(const Code& c) {
    uint64_t h1 = 1469598103934665603ull, h2 = 0x84222325cbf29ce4ull;
    for (uint32_t x : c) {
        h1 ^= x + 0x9e3779b97f4a7c15ull + (h1 << 6) + (h1 >> 2);
        h1 *= 1099511628211ull;
        h2 += x;
        h2 *= 0xff51afd7ed558ccdull;
        h2 ^= h2 >> 33;
    }
    return Key{h1, h2 ^ c.size()};
}

#endif

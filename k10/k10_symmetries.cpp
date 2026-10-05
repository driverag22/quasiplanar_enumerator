#include "hds_quasiplanar.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <tuple>
#include <string>
#include <utility>
#include <nlohmann/json.hpp>

// Computes the strong automorphism group of each K10 drawing, i.e., all
// vertex relabelings that yield a strongly isomorphic drawing.
//
// Instead of relabeling and calling are_isomorphic (iso.h) for each of the
// 10! permutations, we run the same simultaneous-BFS as in iso.h with d1 == d2
// and collect every successful map. An automorphism of the planarization is
// fixed by the image of a single halfedge (and orientation), so there are at
// most 2 * #halfedges candidates and each success is a distinct automorphism.

constexpr int KPLANE_LIMIT = 15;
const std::string DIR = "../quasiDrawings/K10_all_quasi/";
constexpr int MAX_FILES = 9;

struct Automorphism {
    std::vector<std::size_t> perm; // restricted to the original vertices
    bool mirror;                   // true iff orientation reversing
};

// Try to extend h1 -> h2 to an automorphism of the planarization of d.
// On success, phi holds the map on all vertices and crossings.
template < int kplane >
bool extend_automorphism(const Drawing<kplane>& d, const HdsHalfedge* h1, const HdsHalfedge* h2,
                         bool mirror, std::vector<std::size_t>& phi) {
    std::size_t np = d.vertices.size();
    std::size_t n = np + d.crossings.size();
    phi.assign(n, n);
    std::vector<std::size_t> psi(n, n);
    std::size_t a = h1->vertex->label, b = h2->vertex->label;
    if ((a < np) != (b < np)) return false;
    phi[a] = b;
    psi[b] = a;

    std::queue<std::tuple<const HdsHalfedge*, const HdsHalfedge*>> queue;
    queue.emplace(h1, h2);
    while (!queue.empty()) {
        const HdsHalfedge* s1 = std::get<0>(queue.front());
        const HdsHalfedge* s2 = std::get<1>(queue.front());
        queue.pop();
        const HdsHalfedge* t1 = s1;
        const HdsHalfedge* t2 = s2;
        do {
            if (t1->edge->ncr != t2->edge->ncr) return false;
            std::size_t i1 = t1->twin->vertex->label;
            std::size_t i2 = t2->twin->vertex->label;
            if ((i1 < np) != (i2 < np)) return false;
            if (phi[i1] == n) {
                if (psi[i2] < n) return false;
                phi[i1] = i2;
                psi[i2] = i1;
                queue.emplace(t1->twin, t2->twin);
            } else if (phi[i1] != i2 || psi[i2] != i1) {
                return false;
            }
            t1 = (mirror ? t1->twin->prev : t1->next->twin);
            t2 = t2->next->twin;
        } while (t1 != s1 && t2 != s2);
        if (t1 != s1 || t2 != s2) return false;
    }
    for (std::size_t z = 0; z < n; ++z)
        if (phi[z] == n || psi[phi[z]] != z)
            throw std::runtime_error("internal error in automorphism search");
    return true;
}

template < int kplane >
std::vector<Automorphism> automorphisms(const Drawing<kplane>& d) {
    std::vector<Automorphism> result;
    std::size_t np = d.vertices.size();
    const HdsHalfedge* h1 = d.vertices[0].halfedge;
    std::vector<std::size_t> phi;
    for (bool mirror : {false, true}) {
        for (std::size_t j = 0; j < np; ++j) {
            const HdsHalfedge* h2 = d.vertices[j].halfedge;
            const HdsHalfedge* s = h2;
            do {
                if (extend_automorphism(d, h1, s, mirror, phi))
                    result.push_back({std::vector<std::size_t>(phi.begin(), phi.begin() + np), mirror});
                s = s->next->twin;
            } while (s != h2);
        }
    }
    return result;
}

// Independent sanity check: the vertex permutation must preserve which
// pairs of (abstract) edges cross.
template < int kplane >
std::set<std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t>>>
crossing_pairs(const Drawing<kplane>& d, const std::vector<std::size_t>& perm) {
    std::size_t np = d.vertices.size();
    auto key = [&](const HdsEdge* e) {
        std::size_t u = perm[e->u], v = perm[e->v];
        return std::make_pair(std::min(u, v), std::max(u, v));
    };
    std::set<std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t>>> res;
    for (const auto& h : d.halfedges) {
        if (h.vertex->label < np) continue;
        auto e1 = key(h.edge), e2 = key(h.next->twin->edge);
        res.insert(std::minmax(e1, e2));
    }
    return res;
}

std::string cycle_notation(const std::vector<std::size_t>& perm) {
    std::vector<bool> seen(perm.size(), false);
    std::string out;
    for (std::size_t i = 0; i < perm.size(); ++i) {
        if (seen[i] || perm[i] == i) { seen[i] = true; continue; }
        out += "(";
        std::size_t j = i;
        bool first = true;
        while (!seen[j]) {
            seen[j] = true;
            if (!first) out += " ";
            out += std::to_string(j);
            first = false;
            j = perm[j];
        }
        out += ")";
    }
    return out.empty() ? "id" : out;
}

int main() {
    for (int i = 0; i < MAX_FILES; ++i) {
        std::string filename = DIR + std::to_string(i) + ".json";
        std::ifstream input_file(filename);
        if (!input_file) {
            std::cout << "=== " << filename << ": not found, skipping\n\n";
            continue;
        }
        nlohmann::json import_data;
        input_file >> import_data;
        input_file.close();

        Drawing<KPLANE_LIMIT> d(import_data);
        std::size_t np = d.vertices.size();
        std::vector<Automorphism> autos = automorphisms(d);

        std::vector<std::size_t> id(np);
        for (std::size_t v = 0; v < np; ++v) id[v] = v;
        auto base = crossing_pairs(d, id);

        std::size_t n_rot = 0, n_mir = 0;
        for (const auto& a : autos) (a.mirror ? n_mir : n_rot)++;

        std::cout << "=== " << filename << ": " << d.crossings.size() << " crossings, |Aut| = "
                  << autos.size() << " (" << n_rot << " orientation preserving, "
                  << n_mir << " orientation reversing)\n";

        for (const auto& a : autos) {
            bool ok = (crossing_pairs(d, a.perm) == base);
            std::cout << "  " << (a.mirror ? "[mirror] " : "[rotate] ") << cycle_notation(a.perm)
                      << (ok ? "" : "   !! crossing pairs NOT preserved") << "\n";
        }

        // vertex orbits under the automorphism group
        std::vector<std::size_t> orbit(np, np);
        std::size_t n_orbits = 0;
        for (std::size_t v = 0; v < np; ++v) {
            if (orbit[v] != np) continue;
            std::vector<std::size_t> stack = {v};
            orbit[v] = n_orbits;
            while (!stack.empty()) {
                std::size_t x = stack.back();
                stack.pop_back();
                for (const auto& a : autos)
                    if (orbit[a.perm[x]] == np) {
                        orbit[a.perm[x]] = n_orbits;
                        stack.push_back(a.perm[x]);
                    }
            }
            ++n_orbits;
        }
        std::cout << "  vertex orbits:";
        for (std::size_t o = 0; o < n_orbits; ++o) {
            std::cout << " {";
            bool first = true;
            for (std::size_t v = 0; v < np; ++v)
                if (orbit[v] == o) {
                    std::cout << (first ? "" : ",") << v;
                    first = false;
                }
            std::cout << "}";
        }
        std::cout << "\n\n";
    }
    return 0;
}

// Counts the simple k-plane drawings (or, with -DQUASI, the simple quasiplane drawings with at most k crossings per edge)
// of K_n, of K_{a,b}, or of K_n minus a set F of edges, up to strong isomorphism.
//
// Build: g++ -O3 -std=c++17 -Iinclude -pthread -DKLIM=<k> [-DQUASI] validation/enum.cpp -o validation/.bin/enum_k<k>
// Run:   validation/.bin/enum_k3 K 7 [threads] [iso]                      K_7
//        validation/.bin/enum_k3 B 3 4 [threads] [iso]                    K_{3,4}
//        validation/.bin/enum_k3 K 9 minus 1-2,2-3,3-4 [threads] [iso]    K_9 minus a path with three edges
//   iso: also check the classes against the isomorphism test of iso.h, which the other programs use
//
// Method. The graph is G = K_n - F on the vertices 0, ..., n-1. The star at vertex 0, to all its neighbours N(0) in
// increasing order, is drawn first, with a fixed rotation. In a simple drawing of G, the edges at 0 do not cross, so a
// permutation p of N(0) makes the rotation at 0 equal to the fixed one, and turns G into K_n - p(F). So every drawing
// of G is, up to relabelling, a drawing of K_n - p(F) with the fixed star, for one of the distinct images p(F) of F,
// and we enumerate these for all of them. The other edges are drawn in all possible ways, in colex order (by larger
// endpoint, then by smaller endpoint), each starting at an endpoint that is already drawn.
// The drawings are compared by their strong canonical codes (canon.h), which do not use vertex labels.
// K_{a,b} (a <= b) is K_n minus the edges inside its two sides, where 0, b+1, ..., n-1 form the side of size a and
// 1, ..., b the side of size b; then F is invariant under the permutations of N(0) = {1, ..., b}.
#ifdef QUASI
#include "hds_quasiplanar.h"
#else
#include "hds_kplanar.h"
#endif
#include "iso.h"
#include "canon.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>

#ifndef KLIM
#define KLIM 3
#endif
typedef Drawing<KLIM> D;
typedef std::pair<std::size_t, std::size_t> E;   // an edge (u, v); when drawn, u is the endpoint that is already drawn

// draws edges[i..j-1] in all possible ways on d, and calls done for every resulting drawing
template <class F>
void draw_all(D& d, const std::vector<E>& edges, std::size_t i, std::size_t j, const F& done) {
    if (i == j) return done(d);
    const auto [a, b] = edges[i];
    for (HdsPath p = d.first_path(a, b); !p.empty();) {
        d.add_edge(p, b);
        draw_all(d, edges, i + 1, j, done);
        p = d.edges.back().built;
        d.remove_edge();
        if (!d.next_path(p, b)) break;
    }
}

// K_n - F with the star at 0 drawn, and the remaining edges in drawing order
std::pair<D, std::vector<E>> start(std::size_t n, const std::set<E>& F) {
    std::vector<bool> drawn(n, false);
    D d(n);
    drawn[0] = true;
    for (std::size_t v = 1; v < n; ++v) {
        if (F.count({0, v})) continue;
        if (d.edges.empty()) d.add_first_edge(0, v);
        else d.add_edge(d.first_path(0, v), v);
        drawn[v] = true;
    }
    std::vector<E> rest;
    for (std::size_t v = 2; v < n; ++v)
        for (std::size_t u = 1; u < v; ++u) {
            if (F.count({u, v})) continue;
            if (!drawn[u] && !drawn[v]) throw std::runtime_error("edge with no drawn endpoint");
            rest.push_back(drawn[u] ? E{u, v} : E{v, u});
            drawn[u] = drawn[v] = true;
        }
    return {d, rest};
}

int main(int argc, char** argv) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<std::string> args(argv + 1, argv + argc);
    if (args.size() < 2) {
        std::cerr << "usage: " << argv[0] << " K n [minus u-v,...] | B a b   [threads] [iso]\n";
        return 2;
    }
    std::size_t n, i;   // i: index of the first optional argument
    std::set<E> F;
    std::string name;
    if (args[0] == "B") {
        std::size_t a = std::stoul(args[1]), b = std::stoul(args[2]);
        if (a > b) std::swap(a, b);
        n = a + b;
        name = "K" + std::to_string(a) + "," + std::to_string(b);
        auto big_side = [&](std::size_t v) { return 1 <= v && v <= b; };
        for (std::size_t v = 1; v < n; ++v)
            for (std::size_t u = 0; u < v; ++u)
                if (big_side(u) == big_side(v)) F.insert({u, v});
        i = 3;
    } else {
        n = std::stoul(args[1]);
        name = "K" + std::to_string(n);
        i = 2;
        if (args.size() > 3 && args[2] == "minus") {   // a comma-separated list u-v,...
            const std::string s = args[3] + ",";
            for (std::size_t pos = 0, c; (c = s.find(',', pos)) != std::string::npos; pos = c + 1) {
                const std::string t = s.substr(pos, c - pos);
                F.insert(std::minmax(std::stoul(t.substr(0, t.find('-'))), std::stoul(t.substr(t.find('-') + 1))));
            }
            name += " minus " + args[3];
            i = 4;
        }
    }
    const std::size_t threads = args.size() > i ? std::stoul(args[i]) : 1;
    const bool check_iso = args.size() > i + 1 && args[i + 1] == "iso";
#ifdef QUASI
    const std::string cls = "quasiplane (at most " + std::to_string(KLIM) + " crossings per edge)";
#else
    const std::string cls = std::to_string(KLIM) + "-plane";
#endif

    // the distinct images of F under the permutations of N(0)
    std::vector<std::size_t> nb;
    for (std::size_t v = 1; v < n; ++v)
        if (!F.count({0, v})) nb.push_back(v);
    std::set<std::set<E>> images;
    std::vector<std::size_t> p = nb;   // nb[k] is mapped to p[k]
    do {
        std::vector<std::size_t> img(n);
        for (std::size_t v = 0; v < n; ++v) img[v] = v;
        for (std::size_t k = 0; k < nb.size(); ++k) img[nb[k]] = p[k];
        std::set<E> f;
        for (const auto& [u, v] : F) f.insert(std::minmax(img[u], img[v]));
        images.insert(f);
    } while (std::next_permutation(p.begin(), p.end()));

    // tasks: the partial drawings after the star and the next `split` edges; many small tasks keep all threads busy
    std::vector<D> tasks;
    std::vector<std::vector<E>> task_rest;
    for (std::size_t split = 0;; ++split) {
        tasks.clear();
        task_rest.clear();
        bool all = true;   // split covers all remaining edges, for every image
        for (const auto& f : images) {
            auto [d, rest] = start(n, f);
            const std::size_t s = std::min(split, rest.size());
            all = all && s == rest.size();
            draw_all(d, rest, 0, s, [&](D& x) {
                tasks.push_back(x);
                task_rest.emplace_back(rest.begin() + s, rest.end());
            });
        }
        // take first split with enough tasks
        if (threads == 1 || all || tasks.empty() || tasks.size() >= 1024 * threads) break;
    }

    // every thread keeps the strong codes it has seen, with a representative of each class if iso is checked
    struct Table {
        std::unordered_map<Key, std::size_t, KeyHash> classes;   // strong code -> index of the representative
        std::vector<D> reps;
        uint64_t labelled = 0, iso_mismatch = 0;
    };
    std::vector<Table> tables(threads);
    std::atomic<std::size_t> next(0);
    auto work = [&](Table& T) {
        Canonizer<D> cz;
        for (std::size_t t = next++; t < tasks.size(); t = next++)
            draw_all(tasks[t], task_rest[t], 0, task_rest[t].size(), [&](D& d) {
                ++T.labelled;
                const Key k = fingerprint(cz.strong(d));
                const auto it = T.classes.find(k);
                if (it == T.classes.end()) {
                    T.classes.emplace(k, T.reps.size());
                    if (check_iso) T.reps.push_back(d);
                } else if (check_iso && !are_isomorphic(T.reps[it->second], d)) {
                    ++T.iso_mismatch;   // same code, but not isomorphic by iso.h
                }
            });
    };
    std::vector<std::thread> pool;
    for (auto& T : tables) pool.emplace_back(work, std::ref(T));
    for (auto& th : pool) th.join();

    // merge the tables; with iso, also check that representatives of different classes are not isomorphic by iso.h
    std::unordered_map<Key, const D*, KeyHash> classes;
    uint64_t labelled = 0, iso_mismatch = 0;
    for (const auto& T : tables) {
        labelled += T.labelled;
        iso_mismatch += T.iso_mismatch;
        for (const auto& [k, r] : T.classes) {
            const auto [it, is_new] = classes.emplace(k, check_iso ? &T.reps[r] : nullptr);
            if (!is_new && check_iso && !are_isomorphic(*it->second, T.reps[r])) ++iso_mismatch;
        }
    }
    if (check_iso) {
        std::map<std::vector<std::size_t>, std::vector<const D*>> by_crossings;   // sorted crossing numbers of the edges
        for (const auto& [k, r] : classes) {
            std::vector<std::size_t> c;
            for (const auto& e : r->edges) c.push_back(e.ncr);
            std::sort(c.begin(), c.end());
            by_crossings[c].push_back(r);
        }
        for (const auto& [c, rs] : by_crossings)
            for (std::size_t x = 0; x < rs.size(); ++x)
                for (std::size_t y = x + 1; y < rs.size(); ++y) iso_mismatch += are_isomorphic(*rs[x], *rs[y]);
    }

    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::cout << name << ", simple " << cls << " drawings: " << classes.size() << " up to strong isomorphism ("
              << labelled << " labelled drawings with the star at 0 fixed, " << images.size() << " image"
              << (images.size() == 1 ? "" : "s") << " of the missing edges)";
    if (check_iso) std::cout << ", iso.h mismatches: " << iso_mismatch;
    std::cout << ", " << threads << " thread" << (threads == 1 ? "" : "s") << ", " << secs << " s" << std::endl;
    return iso_mismatch == 0 ? 0 : 1;
}

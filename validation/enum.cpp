// Count all simple k-plane (or quasiplane, -DQUASI) drawings of K_n or K_{a,b} up to strong isomorphism and up to
// orientation-preserving isomorphism (both with relabelling of the vertices by graph automorphisms).
//
// Build: g++ -O3 -Wall -std=c++17 -Iinclude -DKLIM=<k> [-DQUASI] validation/enum.cpp -o validation/.bin/enum_k<k>
// Run:   validation/.bin/enum_k3 K 7 [threads] [split] [flags]       (complete graph K_7)
//        validation/.bin/enum_k3 B 3 4 [threads] [split] [flags]     (complete bipartite K_{3,4})
//   threads: worker threads (default 1); split: edges drawn after the star before the work is distributed (default 2)
//   flags: "iso" = cross-check every drawing against iso.h, "lab" = check that no labelled drawing is found twice
//
// Symmetry reduction (star fixing):
//   K_n: the star at vertex 0 is drawn first, with the fixed rotation produced by first_path. In any simple drawing
//     the star edges at 0 do not cross each other, and relabelling 1..n-1 makes the rotation at 0 equal to the fixed
//     one; any permutation of 1..n-1 is an automorphism of K_n. So every drawing is, up to a label-permuting
//     automorphism and an orientation-preserving homeomorphism, one of those enumerated.
//   K_{a,b} (a <= b), sides A = {0..a-1}, B = {a..a+b-1}: the star from 0 to B is drawn first with a fixed rotation;
//     any permutation of B is an automorphism of K_{a,b}, so the same argument applies. The remaining edges are drawn
//     vertex by vertex for x = 1..a-1: first (a, x), which places x, then (a+j, x).
//   No mirror image is used in the reduction, so the enumerated set contains a representative of every
//   orientation-preserving class; counting distinct oriented canonical codes therefore counts these classes directly.
#ifdef QUASI
#include "hds_quasiplanar.h"
#else
#include "hds_kplanar.h"
#endif
#include "iso.h"
#include "common.h"

#ifndef KLIM
#define KLIM 3
#endif
const std::size_t klim = KLIM;
typedef Drawing<klim> D;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " K n | B a b  [threads] [split] [iso] [lab]\n";
        return 2;
    }
    const bool bip = std::string(argv[1]) == "B";
    int ai = 2;
    std::size_t a = std::stoul(argv[ai++]);
    std::size_t b = bip ? std::stoul(argv[ai++]) : 0;
    if (bip && a > b) std::swap(a, b);
    std::size_t threads = argc > ai ? std::stoul(argv[ai]) : 1;
    ++ai;
    std::size_t split = argc > ai ? std::stoul(argv[ai]) : 2;
    ++ai;
    bool check_iso = false, check_lab = false;
    for (int k = ai; k < argc; ++k) {
        if (std::string(argv[k]) == "iso") check_iso = true;
        if (std::string(argv[k]) == "lab") check_lab = true;
    }
    const std::size_t n = bip ? a + b : a;
    const std::string name = bip ? "K" + std::to_string(a) + "," + std::to_string(b) : "K" + std::to_string(n);
#ifdef QUASI
    const std::string cls = "quasi(klim=" + std::to_string(klim) + ")";
#else
    const std::string cls = std::to_string(klim) + "-planar";
#endif

    // edge list: star first, then the rest
    Edges star, rest;
    if (!bip) {
        for (std::size_t v = 1; v < n; ++v) star.push_back({0, v});
        for (std::size_t v = 2; v < n; ++v)
            for (std::size_t u = 1; u < v; ++u) rest.push_back({u, v});
    } else {
        for (std::size_t j = 0; j < b; ++j) star.push_back({0, a + j});
        for (std::size_t x = 1; x < a; ++x)
            for (std::size_t j = 0; j < b; ++j) rest.push_back({a + j, x});
    }
    if (split > rest.size()) split = rest.size();

    D d(n);
    d.add_first_edge(star[0][0], star[0][1]);
    for (std::size_t i = 1; i < star.size(); ++i) {
        HdsPath p = d.first_path(star[i][0], star[i][1]);
        if (p.empty()) throw std::runtime_error("star");
        d.add_edge(p, star[i][1]);
    }
    std::cout << name << " " << cls << ": star of " << star.size() << " edges, " << rest.size()
              << " further edges, threads " << threads << ", split " << split << (check_iso ? ", iso check" : "")
              << (check_lab ? ", label check" : "") << std::endl;

    uint64_t nodes0 = 0;
    std::vector<D> partial;
    draw_all<D>(d, rest, 0, split, nodes0, [&](D& p) { partial.push_back(p); });
    std::cout << partial.size() << " partial drawings after " << split << " edges (" << now_s() << " s)" << std::endl;

    std::vector<ClassTable<D>> tables(threads);
    std::vector<uint64_t> nodes(threads, 0);
    std::atomic<std::size_t> next(0), finished(0);
    std::mutex out_mtx;
    auto work = [&](std::size_t t) {
        ClassTable<D>& T = tables[t];
        T.check_iso = check_iso;
        T.keep_hashes = check_lab;
        for (std::size_t i = next++; i < partial.size(); i = next++) {
            draw_all<D>(partial[i], rest, split, rest.size(), nodes[t], [&](D& f) { T.add(f); });
            const std::size_t fin = ++finished;
            if (fin % std::max<std::size_t>(1, partial.size() / 20) == 0) {
                std::lock_guard<std::mutex> lk(out_mtx);
                std::cout << "  progress " << fin << "/" << partial.size() << " (" << now_s() << " s)" << std::endl;
            }
        }
    };
    std::vector<std::thread> th;
    for (std::size_t t = 0; t < threads; ++t) th.emplace_back(work, t);
    for (auto& x : th) x.join();
    const double t_enum = now_s();

    // merge
    uint64_t labelled = 0, iso_mismatch = 0, iso_calls = 0, total_nodes = nodes0;
    std::unordered_map<Key, std::pair<bool, uint64_t>, KeyHash> classes;   // strong code -> (achiral, count)
    std::unordered_map<Key, std::pair<std::size_t, std::size_t>, KeyHash> class_rep;   // strong code -> (table, rep)
    std::unordered_set<Key, KeyHash> oriented;
    std::vector<uint64_t> hashes;
    for (std::size_t t = 0; t < threads; ++t) {
        ClassTable<D>& T = tables[t];
        labelled += T.labelled;
        iso_mismatch += T.iso_mismatch;
        iso_calls += T.iso_calls;
        total_nodes += nodes[t];
        oriented.insert(T.oriented.begin(), T.oriented.end());
        hashes.insert(hashes.end(), T.labelled_hashes.begin(), T.labelled_hashes.end());
        for (auto& [code, inf] : T.classes) {
            auto it = classes.find(code);
            if (it == classes.end()) {
                classes.emplace(code, std::make_pair(inf.achiral, inf.count));
                if (check_iso) class_rep.emplace(code, std::make_pair(t, inf.rep));
            } else {
                it->second.second += inf.count;
                if (check_iso) {   // same class found by two threads: iso.h must agree
                    auto [t2, r2] = class_rep[code];
                    ++iso_calls;
                    if (!are_isomorphic(tables[t2].reps[r2], T.reps[inf.rep])) ++iso_mismatch;
                }
            }
        }
    }
    std::size_t achiral = 0;
    for (auto& [code, x] : classes) achiral += x.first;
    uint64_t dup_labelled = 0;
    if (check_lab) {
        std::sort(hashes.begin(), hashes.end());
        for (std::size_t i = 1; i < hashes.size(); ++i) dup_labelled += (hashes[i] == hashes[i - 1]);
    }
    // negative cross-check: representatives of different classes must not be isomorphic by iso.h
    uint64_t iso_neg_pairs = 0;
    if (check_iso) {
        std::map<std::string, std::vector<const D*>> buckets;
        for (auto& [code, tr] : class_rep) {
            const D& r = tables[tr.first].reps[tr.second];
            buckets[invariant(r)].push_back(&r);
        }
        for (auto& [inv, v] : buckets)
            for (std::size_t i = 0; i < v.size(); ++i)
                for (std::size_t j = i + 1; j < v.size(); ++j) {
                    ++iso_neg_pairs;
                    ++iso_calls;
                    if (are_isomorphic(*v[i], *v[j])) ++iso_mismatch;
                }
    }
    const double t_all = now_s();
    const std::size_t C = classes.size();
    std::cout << "RESULT " << name << " " << cls << " strong=" << C << " achiral=" << achiral
              << " oriented=" << oriented.size() << " (2C-A=" << 2 * C - achiral << ")"
              << " labelled=" << labelled << " nodes=" << total_nodes;
    if (check_lab) std::cout << " dup_labelled=" << dup_labelled;
    if (check_iso)
        std::cout << " iso_calls=" << iso_calls << " iso_neg_pairs=" << iso_neg_pairs << " iso_mismatch=" << iso_mismatch;
    std::cout << " t_enum=" << t_enum << "s t_total=" << t_all << "s peak_mem=" << peak_mem_mb() << "MB"
              << std::endl;
    return 0;
}

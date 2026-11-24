// Hoffmann-Reddy-Wang checks for simple 3-plane drawings of K_9 minus 2 or 3 edges and of K_10 minus a perfect
// matching.
//
// Build: g++ -O3 -Wall -std=c++17 -Iinclude validation/hrw.cpp -o validation/.bin/hrw
// Run:   validation/.bin/hrw <case> [threads]     case = adj2 | dis2 | tri | path3 | star3 | p3k2 | 3k2 | k10pm | all
//        validation/.bin/hrw direct <case> [threads]     (independent second method, see run_direct)
//        validation/.bin/hrw direct <case> [threads]     (independent second method, see run_direct)
//
// Method. Let G = K_N - F and let S be a vertex cover of F, so that G - S = K_m (m = N - |S|). Every simple 3-plane
// drawing of G contains a simple 3-plane drawing of K_m on V - S, which (relabelled) is one of the representatives of
// the classes of simple 3-plane drawings of K_m, enumerated here with enum.cpp's method (star fixed, all labelled
// drawings, deduplicated by canonical code). The vertices of S get the labels m..N-1; the endpoints of F outside S
// ("placeholders") are mapped to old labels 0..m-1 in all injective ways, which covers every relabelling of V - S.
// Placements that differ only by a permutation of the new labels m..N-1 are equivalent (the permutation maps the
// drawings extending the base to each other), so only one of them is kept. For every representative and placement,
// all ways of drawing the edges of G at S are enumerated (no pruning), and the drawings of G found are counted up to
// strong isomorphism (canonical codes; the codes do not use labels, so drawings from different placements compare).
#include "hds_kplanar.h"
#include "iso.h"
#include "common.h"
#include <memory>
#include <set>

const std::size_t klim = 3;
typedef Drawing<klim> D;

// pattern: missing edges, vertices >= 100 are placeholders (100 + k), others are new vertices (index j -> label m+j)
struct Case {
    std::string name, desc;
    std::size_t N, m;
    std::vector<std::pair<std::size_t, std::size_t>> missing;   // (new index j, placeholder or new index)
};
const std::size_t P = 100;

std::vector<Case> cases() {
    return {
        {"adj2", "K9 minus two adjacent edges", 9, 8, {{0, P + 0}, {0, P + 1}}},
        {"dis2", "K9 minus two disjoint edges", 9, 7, {{0, P + 0}, {1, P + 1}}},
        {"tri", "K9 minus a triangle", 9, 7, {{0, P + 0}, {1, P + 0}, {0, 1}}},
        {"path3", "K9 minus a path with 3 edges", 9, 7, {{0, P + 0}, {0, 1}, {1, P + 1}}},
        {"star3", "K9 minus a star K_{1,3}", 9, 8, {{0, P + 0}, {0, P + 1}, {0, P + 2}}},
        {"p3k2", "K9 minus P3+K2", 9, 7, {{0, P + 0}, {0, P + 1}, {1, P + 2}}},
        {"3k2", "K9 minus 3K2", 9, 6, {{0, P + 0}, {1, P + 1}, {2, P + 2}}},
        {"k10pm", "K10 minus a perfect matching", 10, 5, {{0, P + 0}, {1, P + 1}, {2, P + 2}, {3, P + 3}, {4, P + 4}}},
    };
}

typedef std::vector<std::pair<std::size_t, std::size_t>> EdgeSet;   // sorted, (min, max)

// all placements of the placeholders, up to permutations of the new vertices
std::vector<EdgeSet> placements(const Case& c) {
    std::size_t nph = 0, nnew = c.N - c.m;
    for (auto& e : c.missing) if (e.second >= P) nph = std::max(nph, e.second - P + 1);
    std::vector<std::size_t> img(nph);
    std::set<EdgeSet> seen;
    std::vector<EdgeSet> res;
    std::function<void(std::size_t, std::vector<bool>&)> rec = [&](std::size_t k, std::vector<bool>& used) {
        if (k == nph) {
            std::vector<std::size_t> perm(nnew);
            for (std::size_t j = 0; j < nnew; ++j) perm[j] = j;
            EdgeSet best;
            bool first = true;
            do {
                EdgeSet f;
                for (auto& e : c.missing) {
                    std::size_t x = c.m + perm[e.first];
                    std::size_t y = e.second >= P ? img[e.second - P] : c.m + perm[e.second];
                    f.push_back({std::min(x, y), std::max(x, y)});
                }
                std::sort(f.begin(), f.end());
                if (first || f < best) best = f;
                first = false;
            } while (std::next_permutation(perm.begin(), perm.end()));
            if (seen.insert(best).second) res.push_back(best);
            return;
        }
        for (std::size_t o = 0; o < c.m; ++o) {
            if (used[o]) continue;
            used[o] = true;
            img[k] = o;
            rec(k + 1, used);
            used[o] = false;
        }
    };
    std::vector<bool> used(c.m, false);
    rec(0, used);
    return res;
}

// representatives of the classes of simple 3-plane drawings of K_m, drawn on vertices 0..m-1 of a drawing with N
// vertices (the others isolated)
std::vector<D> base_classes(std::size_t m, std::size_t N, uint64_t& labelled) {
    Edges star, rest;
    for (std::size_t v = 1; v < m; ++v) star.push_back({0, v});
    for (std::size_t v = 2; v < m; ++v)
        for (std::size_t u = 1; u < v; ++u) rest.push_back({u, v});
    D d(N);
    d.add_first_edge(0, 1);
    for (std::size_t i = 1; i < star.size(); ++i) d.add_edge(d.first_path(star[i][0], star[i][1]), star[i][1]);
    std::set<Key> seen;
    std::vector<D> reps;
    Canonizer<D> cz;
    uint64_t nodes = 0;
    labelled = 0;
    draw_all<D>(d, rest, 0, rest.size(), nodes, [&](D& f) {
        ++labelled;
        Key k = code_str(std::min(cz.canon(f, false), cz.canon(f, true)));
        if (seen.insert(k).second) reps.push_back(f);
    });
    return reps;
}

void run(const Case& c, std::size_t threads) {
    const double t0 = now_s();
    uint64_t base_labelled = 0;
    std::vector<D> reps = base_classes(c.m, c.N, base_labelled);
    std::vector<EdgeSet> pl = placements(c);
    std::cout << c.name << " (" << c.desc << "): base K" << c.m << " has " << reps.size() << " classes ("
              << base_labelled << " labelled), " << pl.size() << " placements, " << reps.size() * pl.size()
              << " extension searches" << std::endl;

    struct Out {
        std::map<Key, std::pair<bool, uint64_t>> classes;   // strong code -> (achiral, count)
        std::map<Key, D*> rep;
        uint64_t labelled = 0, nodes = 0, iso_mismatch = 0;
        std::vector<std::unique_ptr<D>> store;
    };
    std::vector<Out> out(threads);
    std::atomic<std::size_t> next(0);
    const std::size_t jobs = reps.size() * pl.size();
    auto work = [&](std::size_t t) {
        Canonizer<D> cz;
        for (std::size_t w = next++; w < jobs; w = next++) {
            const D& base = reps[w / pl.size()];
            const EdgeSet& f = pl[w % pl.size()];
            std::set<std::pair<std::size_t, std::size_t>> miss(f.begin(), f.end());
            Edges edges;   // edges at the new vertices, from the earlier endpoint
            for (std::size_t x = c.m; x < c.N; ++x)
                for (std::size_t u = 0; u < x; ++u)
                    if (!miss.count({u, x})) edges.push_back({u, x});
            if (edges.empty() || edges[0][0] >= c.m) throw std::runtime_error("bad edge order");
            D d(base);
            draw_all<D>(d, edges, 0, edges.size(), out[t].nodes, [&](D& g) {
                ++out[t].labelled;
                Code co = cz.canon(g, false), cm = cz.canon(g, true);
                Key k = code_str(std::min(co, cm));
                auto it = out[t].classes.find(k);
                if (it == out[t].classes.end()) {
                    out[t].classes[k] = {co == cm, 1};
                    out[t].store.emplace_back(new D(g));
                    out[t].rep[k] = out[t].store.back().get();
                    if (are_isomorphic_oriented(g, g, true) != (co == cm)) ++out[t].iso_mismatch;
                } else {
                    ++it->second.second;
                    if (!are_isomorphic(*out[t].rep[k], g)) ++out[t].iso_mismatch;
                }
            });
        }
    };
    std::vector<std::thread> th;
    for (std::size_t t = 0; t < threads; ++t) th.emplace_back(work, t);
    for (auto& x : th) x.join();

    std::map<Key, std::pair<bool, uint64_t>> classes;
    std::map<Key, D*> rep;
    uint64_t labelled = 0, nodes = 0, iso_mismatch = 0;
    for (auto& o : out) {
        labelled += o.labelled;
        nodes += o.nodes;
        iso_mismatch += o.iso_mismatch;
        for (auto& [k, v] : o.classes) {
            auto it = classes.find(k);
            if (it == classes.end()) {
                classes[k] = v;
                rep[k] = o.rep[k];
            } else {
                it->second.second += v.second;
                if (!are_isomorphic(*rep[k], *o.rep[k])) ++iso_mismatch;
            }
        }
    }
    std::size_t achiral = 0;
    for (auto& [k, v] : classes) achiral += v.first;
    // different classes must not be isomorphic by iso.h
    for (auto i = rep.begin(); i != rep.end(); ++i)
        for (auto j = std::next(i); j != rep.end(); ++j)
            if (are_isomorphic(*i->second, *j->second)) ++iso_mismatch;
    std::size_t idx = 0;
    for (auto& [k, r] : rep) {
        std::ofstream js("logs/hrw_" + c.name + "_" + std::to_string(idx++) + ".json");
        js << r->serialize_to_json().dump(2);
    }
    std::cout << "RESULT " << c.name << " (" << c.desc << ") strong=" << classes.size() << " achiral=" << achiral
              << " oriented=" << 2 * classes.size() - achiral << " labelled_found=" << labelled
              << " nodes=" << nodes << " iso_mismatch=" << iso_mismatch;
    if (!rep.empty()) {
        std::cout << " crossings=";
        for (auto& [k, r] : rep) std::cout << r->crossings.size() << ",";
    }
    std::cout << " time=" << now_s() - t0 << "s peak_mem=" << peak_mem_mb() << "MB" << std::endl;
}

// Independent second method ("direct"): vertex 0 is not covered by F except, for k10pm, by its partner 9. The star
// at 0 (to all its neighbours) is drawn first with a fixed rotation. Relabelling the leaves 1..8 makes the rotation at
// 0 the fixed one and turns F into another edge set F' (still containing 0-9 for k10pm), so the union over all
// distinct images F' of F under permutations of 1..8 covers all drawings. Each G = K_N - F' is enumerated
// exhaustively (no pruning), the edges after the star in colex order.
void run_direct(const std::string& name, std::size_t threads) {
    const double t0 = now_s();
    std::size_t N = 9;
    EdgeSet F;
    if (name == "adj2") F = {{1, 2}, {1, 3}};
    else if (name == "dis2") F = {{1, 2}, {3, 4}};
    else if (name == "tri") F = {{1, 2}, {1, 3}, {2, 3}};
    else if (name == "path3") F = {{1, 2}, {2, 3}, {3, 4}};
    else if (name == "star3") F = {{1, 2}, {1, 3}, {1, 4}};
    else if (name == "p3k2") F = {{1, 2}, {1, 3}, {4, 5}};
    else if (name == "3k2") F = {{1, 2}, {3, 4}, {5, 6}};
    else if (name == "k10pm") { N = 10; F = {{0, 9}, {1, 2}, {3, 4}, {5, 6}, {7, 8}}; }
    else throw std::runtime_error("unknown case " + name);
    std::set<EdgeSet> images;
    std::vector<std::size_t> perm(N);
    for (std::size_t i = 0; i < N; ++i) perm[i] = i;
    do {
        EdgeSet f;
        for (auto& e : F) {
            const std::size_t x = perm[e.first], y = perm[e.second];
            f.push_back({std::min(x, y), std::max(x, y)});
        }
        std::sort(f.begin(), f.end());
        images.insert(f);
    } while (std::next_permutation(perm.begin() + 1, perm.begin() + 9));   // permute 1..8 only
    const std::vector<EdgeSet> imgs(images.begin(), images.end());
    std::cout << "direct " << name << ": " << imgs.size() << " placements of the missing edges" << std::endl;

    struct Out {
        std::map<Key, std::pair<bool, uint64_t>> classes;
        uint64_t labelled = 0, nodes = 0;
    };
    std::vector<Out> out(threads);
    std::atomic<std::size_t> next(0);
    auto work = [&](std::size_t t) {
        Canonizer<D> cz;
        for (std::size_t w = next++; w < imgs.size(); w = next++) {
            const std::set<std::pair<std::size_t, std::size_t>> miss(imgs[w].begin(), imgs[w].end());
            Edges star, rest;
            for (std::size_t v = 1; v < N; ++v)
                if (!miss.count({0, v})) star.push_back({0, v});
            for (std::size_t v = 2; v < N; ++v)
                for (std::size_t u = 1; u < v; ++u)
                    if (!miss.count({u, v})) rest.push_back({u, v});
            D d(N);
            d.add_first_edge(star[0][0], star[0][1]);
            for (std::size_t i = 1; i < star.size(); ++i)
                d.add_edge(d.first_path(star[i][0], star[i][1]), star[i][1]);
            draw_all<D>(d, rest, 0, rest.size(), out[t].nodes, [&](D& g) {
                ++out[t].labelled;
                const Code co = cz.canon(g, false), cm = cz.canon(g, true);
                auto& c = out[t].classes[code_str(std::min(co, cm))];
                c.first = (co == cm);
                ++c.second;
            });
        }
    };
    std::vector<std::thread> th;
    for (std::size_t t = 0; t < threads; ++t) th.emplace_back(work, t);
    for (auto& x : th) x.join();
    std::map<Key, std::pair<bool, uint64_t>> classes;
    uint64_t labelled = 0, nodes = 0;
    for (auto& o : out) {
        labelled += o.labelled;
        nodes += o.nodes;
        for (auto& [k, v] : o.classes) {
            classes[k].first = v.first;
            classes[k].second += v.second;
        }
    }
    std::size_t achiral = 0;
    for (auto& [k, v] : classes) achiral += v.first;
    std::cout << "RESULT direct " << name << " strong=" << classes.size() << " achiral=" << achiral
              << " oriented=" << 2 * classes.size() - achiral << " labelled_found=" << labelled << " nodes=" << nodes
              << " time=" << now_s() - t0 << "s peak_mem=" << peak_mem_mb() << "MB" << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " adj2|dis2|tri|path3|star3|p3k2|3k2|k10pm|all [threads]\n";
        return 2;
    }
    const std::string which = argv[1];
    if (which == "direct") {   // hrw direct <case> [threads]
        run_direct(argv[2], argc > 3 ? std::stoul(argv[3]) : 1);
        return 0;
    }
    const std::size_t threads = argc > 2 ? std::stoul(argv[2]) : 1;
    for (const Case& c : cases())
        if (which == "all" || which == c.name) run(c, threads);
    return 0;
}

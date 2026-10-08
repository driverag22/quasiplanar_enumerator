#ifndef SPLIT_DRIVER_H
#define SPLIT_DRIVER_H

#include "iterative_nested_split.h"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// klim and make_config() (the gadget) are defined by the including iterative_nested_c<k>_split.cpp
using namespace nested_cycle_build;

// Split driver for cluster job arrays (shared by iterative_nested_c<k>_split.cpp, which define klim and
// make_config() with the gadget before including this file).
//
//   search <task_index> <task_count> <split_depth> <passes> <output_dir> [parents_file]
//       search the prefixes (paths of the first split_depth gadget edges) with
//       index % task_count == task_index; passes = 1, 2 or 12. Without parents_file the base
//       drawing is the only parent. New children are appended to output_dir/children.jsonl,
//       Pass 2 complete drawings to output_dir/solutions.jsonl.
//   count <split_depth> <passes> [parents_file]
//       only count the prefixes at split_depth (to choose split_depth and task_count)
//   merge <output_dir> <input_dir>...
//       dedup all children (children.jsonl) of the input dirs by canonical code, write the
//       distinct ones to output_dir/parents.jsonl (one drawing per line), output_dir/parents.txt
//       (usable as parents_file, it just names parents.jsonl) and the union of the extensible.txt files.
//
//   merge_solutions <output_dir> <input_dir>...
//       dedup all Pass 2 complete drawings (solutions.jsonl) up to strong isomorphism (no fixed
//       labels), write output_dir/final.jsonl and final_<k>.json/.graphml (first 1000).
//
// parents_file: one entry per line: "base" (the base drawing), a JSON file with one drawing, or a
// .jsonl file with one drawing per line. Parents are loaded one at a time.
// Environment: REACH_EVERY_EDGE_FROM=<index> enables the reachability check after every edge from
// that index on, FINAL_FACE_EVERY_EDGE_FROM=<index> the early final-face test. All tasks of a run (and any resume) must use the same setting, since pruning
// above split_depth changes the prefix numbering.

// Environment options shared by all drivers (see Config): REACH_EVERY_EDGE_FROM,
// FINAL_FACE_EVERY_EDGE_FROM, PLACEMENT_EVERY_EDGE_FROM, PLACEMENT_FILTER, MATCHING_FIRST,
// PROGRESS_SECONDS, EXPORT_CHILDREN (0 = count children only, for benchmarks).
template <typename C>
void apply_env_options(C& config) {
    // optional: reachability check after every edge from this index on (see Config::reach_every_edge_from)
    if (const char* r = std::getenv("REACH_EVERY_EDGE_FROM")) config.reach_every_edge_from = std::stoul(r);
    if (const char* r = std::getenv("FINAL_FACE_EVERY_EDGE_FROM")) config.final_face_every_edge_from = std::stoul(r);
    if (const char* r = std::getenv("PLACEMENT_EVERY_EDGE_FROM")) config.placement_every_edge_from = std::stoul(r);
    if (const char* r = std::getenv("MATCHING_FIRST")) config.matching_edges_first = std::string(r) == "1";
    if (const char* r = std::getenv("PLACEMENT_FILTER")) config.placement_filter = std::string(r) == "1";
    if (const char* r = std::getenv("EXPORT_CHILDREN")) config.split_export_children = std::string(r) != "0";
    if (const char* r = std::getenv("PROGRESS_SECONDS")) config.progress_seconds = std::stod(r);
    if (const char* r = std::getenv("MCF_PRUNING")) config.enable_mcf_pruning = std::string(r) != "0";
}

std::vector<std::string> read_lines(const std::string& file) {
    std::ifstream in(file);
    if (!in) throw std::runtime_error("cannot open " + file);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) lines.push_back(line);
    }
    return lines;
}

int search(int argc, char** argv, bool count_only) {
    auto config = make_config();
    apply_env_options(config);
    int a = 2;
    if (!count_only) {
        config.task_index = std::stoul(argv[a++]);
        config.task_count = std::stoul(argv[a++]);
    } else {
        config.task_count = 0;
    }
    config.split_depth = std::stoul(argv[a++]);
    std::string passes = argv[a++];
    config.run_pass1 = passes.find('1') != std::string::npos;
    config.run_pass2 = passes.find('2') != std::string::npos;
    config.output_dir = count_only ? std::string("./count_tmp") : std::string(argv[a++]);
    if (a < argc) config.input_parents = read_lines(argv[a++]);
    config.split_mode = true;
    config.export_files = true; // Pass 2 full solutions

    auto t0 = std::chrono::steady_clock::now();
    NestedCycleSearcher<klim> searcher(config);
    auto result = searcher.run();
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::cout << "\n==================================================" << std::endl;
    std::cout << (count_only ? "COUNT" : "SEARCH") << " task " << config.task_index << "/" << config.task_count
              << ", split_depth " << config.split_depth << ", passes " << passes << std::endl;
    std::cout << "Units (prefixes) enumerated: " << result.split_units << std::endl;
    std::cout << "Children exported: " << result.children_exported << std::endl;
    std::cout << "Full solutions: " << result.full_solution_count << std::endl;
    std::cout << "Discarded isomorphisms: " << result.discarded_count << std::endl;
    std::cout << "Early pruned: " << result.pruned_early_count << " (placement " << result.pruned_placement_count << ", final face " << result.pruned_final_face_count << ", MCF " << result.pruned_mcf_count
              << ", budget exhausted " << result.mcf_unknown_count << ", mcf time " << result.mcf_seconds << "s)" << std::endl;
    std::cout << "Paths skipped by the placement filter: " << result.filtered_paths_count << std::endl;
    std::cout << "Iso time: " << result.iso_seconds << "s, total time: " << secs << "s" << std::endl;
    std::cout << "==================================================" << std::endl;
    return 0;
}

int merge(int argc, char** argv) {
    namespace fs = std::filesystem;
    auto config = make_config();
    const auto syms = gadget_symmetries(config.cycle_size, config.gadget_edges_builder());
    const fs::path out = fs::absolute(argv[2]);
    fs::create_directories(out);

    std::unordered_map<IsoCode, std::size_t, IsoCodeHash> seen;
    std::set<std::string> extensible;
    std::ofstream parents(out / "parents.jsonl");
    std::size_t total = 0, broken = 0;

    for (int i = 3; i < argc; ++i) {
        const fs::path dir = argv[i];
        std::ifstream in(dir / "children.jsonl");
        for (std::string line; std::getline(in, line);) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            nlohmann::json j;
            try { j = nlohmann::json::parse(line); }
            catch (const nlohmann::json::parse_error&) { ++broken; continue; } // cut off by a kill
            Drawing<klim> d(j);
            ++total;
            if (!seen.emplace(canonical_code(d, syms), seen.size()).second) continue;
            parents << d.serialize_to_json().dump() << "\n";
        }
        if (fs::exists(dir / "extensible.txt"))
            for (const auto& line : read_lines((dir / "extensible.txt").string())) extensible.insert(line);
    }

    std::ofstream(out / "parents.txt") << (out / "parents.jsonl").string() << "\n";
    std::ofstream ext(out / "extensible.txt");
    for (const auto& line : extensible) ext << line << "\n";
    std::cout << "Merged " << total << " children (" << broken << " broken lines skipped) into " << seen.size()
              << " non-isomorphic drawings, " << extensible.size() << " extensible parents" << std::endl;
    return 0;
}

int merge_solutions(int argc, char** argv) {
    // dedup the Pass 2 complete drawings (solutions.jsonl) of the input dirs up to strong isomorphism
    namespace fs = std::filesystem;
    const fs::path out = fs::absolute(argv[2]);
    fs::create_directories(out);
    std::unordered_map<IsoCode, std::size_t, IsoCodeHash> seen;
    std::ofstream all(out / "final.jsonl");
    std::size_t total = 0, broken = 0;
    for (int i = 3; i < argc; ++i) {
        std::ifstream in(fs::path(argv[i]) / "solutions.jsonl");
        for (std::string line; std::getline(in, line);) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            nlohmann::json j;
            try { j = nlohmann::json::parse(line); }
            catch (const nlohmann::json::parse_error&) { ++broken; continue; }
            Drawing<klim> d(j);
            ++total;
            if (!seen.emplace(canonical_code_free(d), seen.size()).second) continue;
            const std::size_t k = seen.size() - 1;
            all << line << "\n";
            if (k < 1000) { // individual files for viewing, as long as there are not too many
                const std::string base = (out / ("final_" + std::to_string(k))).string();
                std::ofstream(base + ".json") << d.serialize_to_json().dump(4);
                std::ofstream gml(base + ".graphml");
                d.graphml_output(gml);
            }
        }
    }
    std::cout << "Merged " << total << " complete drawings (" << broken << " broken lines skipped) into "
              << seen.size() << " non-isomorphic final drawings" << std::endl;
    return 0;
}

// Consistency check of a drawing for any klim (Drawing::is_valid only follows edges through at most 2 crossings):
// Euler's formula for the planarization (connected, isolated vertices ignored), every edge followed through its
// crossings reaches its other endpoint with at most klim crossings and ncr in [#crossings, klim], every crossing is
// between two different, non-adjacent edges, and no two edges cross twice. Returns "" if fine, else a description.
std::string drawing_problems(const Drawing<klim>& D) {
    const std::size_t nv = D.vertices.size();
    std::size_t isolated = 0;
    for (const auto& v : D.vertices) if (!v.halfedge) ++isolated;
    std::vector<int> face;
    const std::size_t nf = NestedCycleSearcher<klim>::label_faces(D, face);
    const long V = static_cast<long>(nv - isolated + D.crossings.size()), E = static_cast<long>(D.halfedges.size() / 2);
    if (V - E + static_cast<long>(nf) != 2) return "Euler formula fails";
    std::set<const HdsEdge*> checked;
    for (const auto& h : D.halfedges) {
        const HdsEdge* e = h.edge;
        if (h.twin->vertex->label != e->u || checked.count(e)) continue; // start at the halfedge leaving u
        checked.insert(e);
        const HdsHalfedge* cur = &h;
        std::size_t crossings = 0;
        while (cur->vertex->label >= nv) { // a crossing: continue straight through it
            cur = cur->next->twin->next;
            if (cur->edge != e) return "edge does not continue straight through a crossing";
            if (++crossings > klim) return "edge with more than klim crossings";
        }
        if (cur->vertex->label != e->v) return "edge does not end at its endpoint";
        if (e->ncr < crossings || e->ncr > klim) return "inconsistent ncr";
    }
    if (checked.size() != D.edges.size()) return "edge without a start halfedge";
    std::set<std::pair<const HdsEdge*, const HdsEdge*>> crossing_pairs;
    for (const auto& c : D.crossings) {
        const HdsHalfedge* a = c.halfedge;
        const HdsEdge* e1 = a->edge;
        const HdsEdge* e2 = a->next->edge;
        if (e1 == e2) return "edge crosses itself";
        if (e1->u == e2->u || e1->u == e2->v || e1->v == e2->u || e1->v == e2->v) return "adjacent edges cross";
        if (!crossing_pairs.insert({std::min(e1, e2), std::max(e1, e2)}).second) return "two edges cross twice";
    }
    return "";
}

// debug: crossing sequence of every edge, keyed by its sorted endpoints, read from the smaller endpoint
std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<std::size_t, std::size_t>>> crossing_sequences(const Drawing<klim>& D) {
    const std::size_t nv = D.vertices.size();
    auto key = [](const HdsEdge* e) { return std::make_pair(std::min(e->u, e->v), std::max(e->u, e->v)); };
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::pair<std::size_t, std::size_t>>> r;
    for (const auto& h : D.halfedges) {
        const HdsEdge* e = h.edge;
        const std::size_t from = std::min(e->u, e->v);
        if (h.twin->vertex->label != from) continue; // start at the halfedge leaving the smaller endpoint
        auto& seq = r[key(e)];
        if (!seq.empty()) continue;
        const HdsHalfedge* cur = &h;
        while (cur->vertex->label >= nv) {
            seq.push_back(key(cur->next->edge));
            cur = cur->next->twin->next;
        }
        if (seq.empty()) seq.clear();
    }
    for (const auto& e : D.edges) r[key(&e)]; // uncrossed edges too
    return r;
}

// debug: cyclic order of the incident edges (sorted endpoints) at every real vertex
std::map<std::size_t, std::vector<std::pair<std::size_t, std::size_t>>> rotations(const Drawing<klim>& D) {
    std::map<std::size_t, std::vector<std::pair<std::size_t, std::size_t>>> r;
    for (const auto& v : D.vertices) {
        if (!v.halfedge) continue;
        const HdsHalfedge* h = v.halfedge;
        do {
            r[v.label].push_back({std::min(h->edge->u, h->edge->v), std::max(h->edge->u, h->edge->v)});
            h = h->next->twin;
        } while (h != v.halfedge);
    }
    return r;
}

// is the cyclic sequence a equal to b restricted to the elements of a (as cyclic sequences)?
template <typename T>
bool cyclic_restriction_equal(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> rb;
    for (const T& x : b) if (std::find(a.begin(), a.end(), x) != a.end()) rb.push_back(x);
    if (rb.size() != a.size()) return false;
    if (a.empty()) return true;
    for (std::size_t s = 0; s < rb.size(); ++s) {
        bool eq = true;
        for (std::size_t i = 0; i < a.size() && eq; ++i) eq = a[i] == rb[(s + i) % rb.size()];
        if (eq) return true;
    }
    return false;
}

// Saturation test of the complete 3-cycle drawings (D1 = base cycle 0..c-1, D2 = first new cycle, D3 = second new cycle,
// uncrossed). Pass 1 on the base drawing (split like search); at every leaf whose reduced drawing is one of the finishable
// intermediate drawings, Pass 2 runs on the complete (unreduced) drawing. Every completion gets an uncrossable star in the
// outer face of D3, and every missing edge between two different cycles or within D2 (pairs within D1 and within D3
// count as present) is tested with first_path (all routes with <= klim crossings). Output per distinct complete drawing
// (strong isomorphism, per task) in <outdir>/saturation.jsonl, summary in <outdir>/summary.txt.
int saturate(int argc, char** argv) {
    namespace fs = std::filesystem;
    auto config = make_config();
    apply_env_options(config);
    const std::size_t c = config.cycle_size;
    const auto syms = gadget_symmetries(c, config.gadget_edges_builder());
    config.task_index = std::stoul(argv[2]);
    config.task_count = std::stoul(argv[3]);
    config.split_depth = std::stoul(argv[4]);
    config.output_dir = argv[5];
    config.split_mode = true;
    config.run_pass1 = true;
    config.run_pass2 = false;
    config.split_export_children = false;
    // SAT_CYCLES=4: four cycles D1..D4 (at every leaf a second Pass 1 adds the crossable D3, then Pass 2 the uncrossed
    // D4), so that D2 and D3 are both middle cycles; default 3 (D1, D2, D3)
    const std::size_t ncycles = std::getenv("SAT_CYCLES") ? std::stoul(std::getenv("SAT_CYCLES")) : 3;
    if (ncycles != 3 && ncycles != 4) throw std::runtime_error("SAT_CYCLES must be 3 or 4");
    // the later cycles and the outer star, so that no copy needs a JSON round trip
    config.reserve_vertices = (ncycles - 2) * c + 1;
    fs::create_directories(config.output_dir);
    // STOP_ON_UNSATURATED=<marker file> (sweeps): stop at the first complete drawing that is not saturated and
    // append to the marker; a task whose marker already exists (another task of the same gadget stopped) skips
    const char* stop_marker = std::getenv("STOP_ON_UNSATURATED");
    if (stop_marker && fs::exists(stop_marker)) {
        std::cout << "SKIP: " << stop_marker << " exists, the gadget is not maximal" << std::endl;
        return 0;
    }

    // finishable intermediate drawings = parents (meta.parent "<path>#<line>") of the 3-cycle final drawings
    // intermediate file "-": no filter, Pass 2 runs on the complete drawing at every leaf
    const bool filter = std::string(argv[6]) != "-";
    std::set<std::size_t> finishable_lines;
    if (filter) {
        std::ifstream in(argv[7]);
        for (std::string line; std::getline(in, line);) {
            if (line.empty()) continue;
            std::string parent = nlohmann::json::parse(line)["meta"]["parent"];
            finishable_lines.insert(std::stoul(parent.substr(parent.rfind('#') + 1)));
        }
    }
    std::unordered_map<IsoCode, std::size_t, IsoCodeHash> finishable; // code -> line in the intermediate drawings
    if (filter) {
        std::ifstream in(argv[6]);
        std::size_t k = 0;
        for (std::string line; std::getline(in, line); ++k)
            if (!line.empty() && finishable_lines.count(k))
                finishable.emplace(canonical_code(Drawing<klim>(nlohmann::json::parse(line)), syms), k);
    }
    // per finishable class: reached by a complete drawing at all, and by one that can itself be completed
    std::set<std::size_t> classes_reached, classes_completed;
    std::cout << "finishable intermediate drawings: " << finishable_lines.size() << " lines, " << finishable.size() << " codes" << std::endl;

    std::size_t leaves = 0, matching = 0, completions = 0, distinct = 0, saturated = 0, errors = 0;
    std::size_t max_braid = 0, max_total = 0; // greedy augmentation maxima
    std::map<std::string, std::size_t> type_count; // insertable edge types over the distinct drawings
    std::unordered_set<IsoCode, IsoCodeHash> seen;
    std::ofstream out(fs::path(config.output_dir) / "saturation.jsonl", std::ios::app);

    // debug: DEBUG_FIND_LEAF=<drawing.json> reports when a leaf isomorphic to that drawing is reached
    std::optional<IsoCode> debug_target;
    if (const char* f = std::getenv("DEBUG_FIND_LEAF"))
        debug_target = canonical_code(Drawing<klim>(nlohmann::json::parse(std::ifstream(f))), syms);
    if (debug_target && std::getenv("DEBUG_PRUNE")) { // report placement prunes of restrictions of the target leaf
        const Drawing<klim> td(nlohmann::json::parse(std::ifstream(std::getenv("DEBUG_FIND_LEAF"))));
        const auto target = crossing_sequences(td);
        const auto target_rot = rotations(td);
        config.placement_prune_callback = [target, target_rot, &config](const Drawing<klim>& d, const char* reason) {
            for (const auto& [v, rot] : rotations(d)) {
                auto it = target_rot.find(v);
                if (it == target_rot.end() || !cyclic_restriction_equal(rot, it->second)) return;
            }
            const auto cur = crossing_sequences(d);
            for (const auto& [k, seq] : cur) {
                auto it = target.find(k);
                if (it == target.end()) return;
                std::vector<std::pair<std::size_t, std::size_t>> restricted;
                for (const auto& x : it->second) if (cur.count(x)) restricted.push_back(x);
                if (restricted != seq) return;
            }
            static int n = 0;
            std::cout << "DEBUG: " << reason << " prune of a restriction of the target, " << d.edges.size() << " edges" << std::endl;
            std::ofstream(fs::path(config.output_dir) / ("pruned_" + std::string(reason) + "_" + std::to_string(n++) + ".json")) << d.serialize_to_json().dump();
        };
    }
    // debug: CHECK_ROUNDTRIP=1 serializes every leaf and complete drawing, reads it back and compares the rotation
    // systems (with crossings) with the original
    const bool check_roundtrip = std::getenv("CHECK_ROUNDTRIP") != nullptr;
    std::size_t roundtrip_checked = 0, roundtrip_failed = 0;
    auto roundtrip = [&](const Drawing<klim>& x) {
        ++roundtrip_checked;
        std::ostringstream a, b;
        a << x;
        try { b << Drawing<klim>(nlohmann::json::parse(x.serialize_to_json().dump())); }
        catch (const std::exception& ex) { b << "exception: " << ex.what(); }
        if (a.str() != b.str()) { ++roundtrip_failed; std::cout << "ROUNDTRIP MISMATCH" << std::endl; }
    };
    // a complete drawing with the cycles starting at starts[0] = 0 (D1), starts[1], ..., the last one uncrossed:
    // outer star, saturation test, record
    auto complete = [&](const Drawing<klim>& D, const std::vector<std::size_t>& starts) {
        ++completions;
        // exact copy; the last (reserved, still isolated) vertex becomes the outer star
        Drawing<klim> E(D);
        const std::size_t star = E.vertices.size() - 1;
        if (E.vertices[star].halfedge) { ++errors; return; }
        // outer face of the last cycle: bounded by exactly its c edges (all other vertices are on the other side)
        std::vector<int> face;
        const std::size_t nf = NestedCycleSearcher<klim>::label_faces(E, face);
        std::vector<std::vector<HdsHalfedge*>> fh(nf);
        for (auto& h : E.halfedges) fh[face[h.label]].push_back(&h);
        const std::size_t k = starts.size(), nmB = starts.back();
        auto inD3 = [&](std::size_t x) { return x >= nmB && x < nmB + c; };
        HdsHalfedge* start = nullptr;
        std::size_t outer = 0;
        for (std::size_t f = 0; f < nf; ++f) {
            bool only_d3 = fh[f].size() == c;
            for (auto* h : fh[f]) only_d3 = only_d3 && inD3(h->edge->u) && inD3(h->edge->v);
            if (only_d3) { ++outer; start = fh[f][0]; }
        }
        if (outer != 1) { ++errors; std::cerr << "outer face of the last cycle not unique: " << outer << std::endl; return; }
        // uncrossable star in that face (as for the black regions of the reduction)
        HdsHalfedge* j = start;
        HdsHalfedge* x = E.add_edge(HdsPath({j, nullptr}), star, klim);
        for (;;) {
            j = j->next->twin->next;
            if (j == start) break;
            E.add_edge(HdsPath({j, x}), star, klim);
        }
        // missing edges: different cycles, or both in a middle cycle; pairs within the first / last cycle count as present
        auto cycle_of = [&](std::size_t v) { int r = 0; for (std::size_t q = 0; q < k; ++q) if (v >= starts[q] && v < starts[q] + c) r = static_cast<int>(q) + 1; return r; };
        auto pos_of = [&](std::size_t v) { return v - starts[cycle_of(v) - 1]; };
        auto middle = [&](int cyc) { return cyc > 1 && cyc < static_cast<int>(k); };
        std::set<std::pair<std::size_t, std::size_t>> adj;
        for (const auto& e : D.edges) adj.insert({std::min(e.u, e.v), std::max(e.u, e.v)});
        std::vector<std::size_t> verts;
        for (std::size_t q = 0; q < k; ++q) for (std::size_t i = 0; i < c; ++i) verts.push_back(starts[q] + i);
        std::sort(verts.begin(), verts.end());
        nlohmann::json insertable = nlohmann::json::array();
        for (std::size_t a = 0; a < verts.size(); ++a)
            for (std::size_t b = a + 1; b < verts.size(); ++b) {
                const std::size_t u = verts[a], v = verts[b];
                const int cu = cycle_of(u), cv = cycle_of(v);
                if (cu == cv && !middle(cu)) continue;
                if (adj.count({u, v})) continue;
                if (!E.first_path(u, v, 0).empty())
                    insertable.push_back({cu, pos_of(u), cv, pos_of(v)});
            }
        if (check_roundtrip) roundtrip(E);
        const IsoCode ecode = canonical_code_free(E);
        if (!seen.insert(ecode).second) return;
        ++distinct;
        if (const std::string pr = drawing_problems(E); !pr.empty()) { ++errors; std::cerr << "complete drawing: " << pr << std::endl; }
        if (!insertable.empty()) { // insert the first insertable edge into a copy and check the result as well
            const auto& t = insertable[0];
            auto label = [&](int cyc, std::size_t pos) { return starts[cyc - 1] + pos; };
            const std::size_t u = label(t[0], t[1]), v = label(t[2], t[3]);
            Drawing<klim> F(E);
            F.add_edge(F.first_path(u, v, 0), v, 0);
            if (const std::string pr = drawing_problems(F); !pr.empty()) { ++errors; std::cerr << "after insertion: " << pr << std::endl; }
        }
        if (insertable.empty()) ++saturated;
        std::set<std::string> types;
        for (const auto& t : insertable) types.insert("D" + std::to_string(t[0].get<int>()) + "-D" + std::to_string(t[2].get<int>()));
        for (const auto& t : types) ++type_count[t];
        // greedy augmentation: add insertable edges one at a time until none fits, braid edges (between consecutive
        // cycles) first, then the others; a lower bound for the number of edges that fit simultaneously
        Drawing<klim> G(E);
        std::set<std::pair<std::size_t, std::size_t>> gadj = adj;
        auto add_greedy = [&](bool braid_only) {
            nlohmann::json added = nlohmann::json::array();
            for (bool changed = true; changed;) {
                changed = false;
                for (std::size_t a = 0; a < verts.size(); ++a)
                    for (std::size_t b = a + 1; b < verts.size(); ++b) {
                        const std::size_t u = verts[a], v = verts[b];
                        const int cu = cycle_of(u), cv = cycle_of(v);
                        if (cu == cv && !middle(cu)) continue;
                        const bool braid = cu - cv == 1 || cv - cu == 1;
                        if (braid_only != braid) continue;
                        if (gadj.count({u, v})) continue;
                        HdsPath p = G.first_path(u, v, 0);
                        if (p.empty()) continue;
                        G.add_edge(p, v, 0);
                        gadj.insert({u, v});
                        added.push_back({cu, pos_of(u), cv, pos_of(v)});
                        changed = true;
                    }
            }
            return added;
        };
        const nlohmann::json added_braid = add_greedy(true);
        const nlohmann::json added_other = add_greedy(false);
        if (const std::string pr = drawing_problems(G); !pr.empty()) { ++errors; std::cerr << "after augmentation: " << pr << std::endl; }
        max_braid = std::max<std::size_t>(max_braid, added_braid.size());
        max_total = std::max<std::size_t>(max_total, added_braid.size() + added_other.size());

        nlohmann::ordered_json rec;
        rec["code_hash"] = std::to_string(IsoCodeHash()(ecode)); // for deduplication across tasks
        rec["n_insertable"] = insertable.size();
        rec["greedy_braid"] = added_braid;  // braid edges added together (greedy, in this order)
        rec["greedy_other"] = added_other;  // then further edges (within middle cycles, between non-consecutive cycles)
        rec["insertable"] = insertable; // [cycle, position, cycle, position], cycles 1 = D1, 2 = D2, ...
        nlohmann::ordered_json labels;
        for (std::size_t q = 0; q < k; ++q) labels["D" + std::to_string(q + 1)] = starts[q];
        labels["outer_star"] = star;
        rec["labels"] = labels;
        rec["drawing"] = E.serialize_to_json();
        out << rec.dump() << "\n";
        out.flush();
        // fail fast (sweeps): one drawing that is not saturated decides that the gadget is not maximal
        if (stop_marker && !insertable.empty()) {
            std::ofstream(stop_marker, std::ios::app) << "task " << config.task_index << "\n";
            std::cout << "STOP: drawing with an insertable edge found (marker " << stop_marker << ")" << std::endl;
            std::exit(0);
        }
    };
    // Pass 2 on d: the uncrossed last cycle is attached to the cycle starting at starts.back()
    auto run_pass2 = [&](const Drawing<klim>& d, std::vector<std::size_t> starts) {
        auto cfg = make_config();
        cfg.split_mode = true;
        cfg.run_pass1 = false;
        cfg.run_pass2 = true;
        cfg.split_depth = 0;
        cfg.task_index = 0;
        cfg.task_count = 1;
        cfg.input_drawing = &d;
        if (const char* r = std::getenv("MCF_PRUNING")) cfg.enable_mcf_pruning = std::string(r) != "0";
        cfg.split_files = false;
        cfg.export_files = false;
        for (std::size_t i = 0; i < c; ++i) cfg.active_cycle.push_back(starts.back() + i);
        cfg.solution_callback = [&](const Drawing<klim>& D, std::size_t nmB) {
            ++completions;
            std::vector<std::size_t> all = starts;
            all.push_back(nmB);
            complete(D, all);
        };
        NestedCycleSearcher<klim>(cfg).run();
    };
    std::size_t inner_leaves = 0; // SAT_CYCLES=4: leaves of the second Pass 1
    config.leaf_callback = [&](const Drawing<klim>& d, std::size_t nmA, const IsoCode& code) {
        ++leaves;
        if (check_roundtrip) roundtrip(d);
        if (debug_target && canonical_code(d, syms) == *debug_target) { // (code is of the reduced drawing)
            std::cout << "DEBUG: target leaf reached" << std::endl;
            std::ofstream(fs::path(config.output_dir) / "target_leaf.json") << d.serialize_to_json().dump();
        }
        auto fc = finishable.find(code);
        if (filter && fc == finishable.end()) return;
        ++matching;
        if (filter) classes_reached.insert(fc->second);
        const std::size_t completions_before = completions;
        if (ncycles == 3) {
            run_pass2(d, {0, nmA});
        } else {
            // second Pass 1: crossable D3 attached to D2 (from the reserved vertices), Pass 2 at each of its leaves
            auto cfg = make_config();
            apply_env_options(cfg);
            cfg.split_mode = true;
            cfg.run_pass1 = true;
            cfg.run_pass2 = false;
            cfg.split_depth = 0;
            cfg.task_index = 0;
            cfg.task_count = 1;
            cfg.input_drawing = &d;
            cfg.split_files = false;
            cfg.export_files = false;
            cfg.split_export_children = false;
            for (std::size_t i = 0; i < c; ++i) cfg.active_cycle.push_back(nmA + i);
            cfg.leaf_callback = [&](const Drawing<klim>& d2, std::size_t nmB, const IsoCode&) {
                ++inner_leaves;
                run_pass2(d2, {0, nmA, nmB});
            };
            NestedCycleSearcher<klim>(cfg).run();
        }
        if (filter && completions > completions_before) classes_completed.insert(fc->second);
    };

    // debug: DEBUG_PLACEMENT=<drawing.json> checks placements_still_possible on every prefix of the drawing's recipe
    // (remaining recipe edges as the local edges still to place); every prefix can be completed, so false is a bug
    if (const char* dp = std::getenv("DEBUG_PLACEMENT")) {
        const nlohmann::json j = nlohmann::json::parse(std::ifstream(dp));
        const auto& rec = j["drawing_recipe"];
        NestedCycleSearcher<klim> s(config);
        for (std::size_t k = 1; k < rec.size(); ++k) {
            nlohmann::json pj = j;
            pj["drawing_recipe"] = nlohmann::json(std::vector<nlohmann::json>(rec.begin(), rec.begin() + k));
            Drawing<klim> d(pj);
            std::vector<Edge> rest;
            for (std::size_t i = k; i < rec.size(); ++i) rest.push_back({rec[i]["u"].get<std::size_t>(), rec[i]["v"].get<std::size_t>()});
            if (std::getenv("DEBUG_SELFTEST")) { // is every prefix recognized as a restriction of the full drawing?
                static const Drawing<klim> full(j);
                static const auto fseq = crossing_sequences(full);
                static const auto frot = rotations(full);
                bool ok = true;
                for (const auto& [v, rot] : rotations(d)) ok = ok && cyclic_restriction_equal(rot, frot.at(v));
                const auto cur = crossing_sequences(d);
                for (const auto& [kk, seq] : cur) {
                    std::vector<std::pair<std::size_t, std::size_t>> r;
                    for (const auto& x : fseq.at(kk)) if (cur.count(x)) r.push_back(x);
                    ok = ok && r == seq;
                }
                if (!ok) std::cout << "selftest: prefix " << k << " not recognized" << std::endl;
            }
            if (!s.placements_still_possible(d, rest, 0, 0))
                std::cout << "placements_still_possible false after " << k << " edges (next " << rec[k]["u"] << "-" << rec[k]["v"] << ")" << std::endl;
        }
        std::cout << "debug placement done" << std::endl;
        return 0;
    }
    // debug: DEBUG_PLACEMENT_STATE=<drawing.json> runs placements_still_possible on that Pass 1 state with all gadget
    // edges not yet in the drawing as the remaining edges
    if (const char* ds = std::getenv("DEBUG_PLACEMENT_STATE")) {
        Drawing<klim> d(nlohmann::json::parse(std::ifstream(ds)));
        std::set<std::pair<std::size_t, std::size_t>> have;
        for (const auto& e : d.edges) have.insert({std::min(e.u, e.v), std::max(e.u, e.v)});
        std::vector<Edge> rest;
        for (const Edge& e : config.gadget_edges_builder()(c + 1))
            if (!have.count({std::min(e[0], e[1]), std::max(e[0], e[1])})) { rest.push_back(e); std::cout << e[0] << "-" << e[1] << " "; }
        std::cout << "\n" << rest.size() << " remaining edges; possible: "
                  << NestedCycleSearcher<klim>(config).placements_still_possible(d, rest, 0, 0) << std::endl;
        return 0;
    }
    // debug: DEBUG_LEAF=<drawing.json> runs only Pass 2 on that Pass 1 leaf (D1 at 0..c-1, star c, D2 at c+1..2c,
    // reserved vertices for D3 and the outer star at the end) instead of the search
    if (const char* dl = std::getenv("DEBUG_LEAF")) {
        Drawing<klim> d(nlohmann::json::parse(std::ifstream(dl)));
        std::cout << "debug leaf " << dl << ": problems '" << drawing_problems(d) << "'" << std::endl;
        config.leaf_callback(d, c + 1, canonical_code(d, syms));
        std::cout << "complete 3-cycle drawings: " << completions << ", distinct: " << distinct << ", saturated: " << saturated
                  << ", errors: " << errors << std::endl;
        return 0;
    }
    auto t0 = std::chrono::steady_clock::now();
    NestedCycleSearcher<klim> searcher(config);
    auto result = searcher.run();
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::ofstream summary(fs::path(config.output_dir) / "summary.txt");
    for (std::ostream* o : {static_cast<std::ostream*>(&std::cout), static_cast<std::ostream*>(&summary)}) {
        *o << "task " << config.task_index << "/" << config.task_count << ", split_depth " << config.split_depth
           << ", units " << result.split_units << ", time " << secs << "s\n"
           << "leaves (with potential final face): " << leaves << ", with finishable reduction: " << matching << "\n"
           << "cycles: " << ncycles << (ncycles == 4 ? ", leaves of the second Pass 1: " + std::to_string(inner_leaves) : std::string()) << "\n"
           << "complete 3-cycle drawings: " << completions << ", distinct (this task): " << distinct
           << ", saturated: " << saturated << ", not saturated: " << distinct - saturated << ", errors: " << errors << "\n"
           << "finishable classes reached: " << classes_reached.size() << ", of which completed from the complete drawing: "
           << classes_completed.size() << "\n"
           << "greedy augmentation: max braid edges added together " << max_braid << ", max edges in total " << max_total << "\n";
        if (check_roundtrip) *o << "roundtrip checked: " << roundtrip_checked << ", failed: " << roundtrip_failed << "\n";
        for (const auto& [t, n] : type_count) *o << "  drawings with an insertable " << t << " edge: " << n << "\n";
    }
    // class lists for the global evaluation (line numbers in the intermediate drawings)
    std::ofstream cr(fs::path(config.output_dir) / "classes_reached.txt");
    for (std::size_t k : classes_reached) cr << k << (classes_completed.count(k) ? " completed" : " not_completed") << "\n";
    return 0;
}

int main(int argc, char** argv) {
    std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "search" && argc >= 7) return search(argc, argv, false);
    if (mode == "count" && argc >= 4) return search(argc, argv, true);
    if (mode == "merge" && argc >= 4) return merge(argc, argv);
    if (mode == "merge_solutions" && argc >= 4) return merge_solutions(argc, argv);
    if (mode == "saturate" && argc >= 8) return saturate(argc, argv);
    std::cerr << "usage:\n"
              << "  " << argv[0] << " search <task_index> <task_count> <split_depth> <passes> <output_dir> [parents_file]\n"
              << "  " << argv[0] << " count <split_depth> <passes> [parents_file]\n"
              << "  " << argv[0] << " merge <output_dir> <input_dir>...\n"
              << "  " << argv[0] << " merge_solutions <output_dir> <input_dir>...\n"
              << "  " << argv[0] << " saturate <task_index> <task_count> <split_depth> <output_dir> <intermediate.jsonl> <final_3cycles.jsonl>\n";
    return 1;
}

#endif // SPLIT_DRIVER_H

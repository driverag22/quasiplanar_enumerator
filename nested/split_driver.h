#ifndef SPLIT_DRIVER_H
#define SPLIT_DRIVER_H

#include "iterative_nested_split.h"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
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
    const auto syms = gadget_symmetries(config.cycle_size, config.local_edges_builder);
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

int main(int argc, char** argv) {
    std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "search" && argc >= 7) return search(argc, argv, false);
    if (mode == "count" && argc >= 4) return search(argc, argv, true);
    if (mode == "merge" && argc >= 4) return merge(argc, argv);
    if (mode == "merge_solutions" && argc >= 4) return merge_solutions(argc, argv);
    std::cerr << "usage:\n"
              << "  " << argv[0] << " search <task_index> <task_count> <split_depth> <passes> <output_dir> [parents_file]\n"
              << "  " << argv[0] << " count <split_depth> <passes> [parents_file]\n"
              << "  " << argv[0] << " merge <output_dir> <input_dir>...\n"
              << "  " << argv[0] << " merge_solutions <output_dir> <input_dir>...\n";
    return 1;
}

#endif // SPLIT_DRIVER_H

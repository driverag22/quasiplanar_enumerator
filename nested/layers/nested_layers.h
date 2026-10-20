// Nested cycles with k layers: all 3-plane drawings, shared by nested_c12_layers.cpp and nested_c14_layers.cpp
// (which fix the cycle length, the matching, the number of cycles and the split settings).
//
// The graph G_{c,k}: cycles D1, ..., Dk of length c, consecutive cycles connected by the braided matching.
// D1 and Dk are uncrossed (in the full graph the gadgets force this), D2, ..., D(k-1) may be crossed.
//
// The search, for one task of a split run (the tasks can run in parallel, e.g. as a Slurm array):
//   1. Pass 1 on the base drawing (the uncrossed D1) adds the crossable D2. The prefixes (paths of the first
//      split_depth edges) are numbered identically in every task; a task searches those with
//      index % num_tasks == task_index.
//   2. Every intermediate drawing found by Pass 1 (j cycles, the last one crossable, with a potential final face) is
//      reduced to its relevant part (see nested_search.h) and extended: by Pass 1 again for the next cycle as long
//      as j + 1 < k, by Pass 2 for the uncrossed Dk otherwise.
// With k = 2 there is no Pass 1: Pass 2 adds the uncrossed D2 to D1 directly (task 0 only).
//
// Counted: the distinct intermediate drawings with j cycles (up to the symmetries of the gadget on the last
// cycle), those with k-1 cycles that can be finished, and the distinct complete drawings (strong isomorphism).
// The complete drawings are completions of reduced intermediate drawings.
//
// Files of a task directory <run_dir>/t<task_index>:
//   <name>_complete.jsonl   the distinct complete drawings of the task, one per line (<name> = c12_3cycles etc.)
//   intermediate.txt        hashes of the distinct intermediate drawings ("<cycles> <hash>"), of those with k-1
//                           cycles that can be finished ("finishable <hash>"), and "finished" at the end;
//                           needed to count distinct drawings over all tasks
//   done_units.txt          finished prefixes (written by the searcher): a killed task restarted with the same
//                           arguments skips them
// The report combines the tasks and writes <run_dir>/<name>_complete.jsonl (all distinct complete drawings).

#include "nested_search.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>
#include <nlohmann/json.hpp>

namespace nested_layers {

using namespace nested_search;
namespace fs = std::filesystem;

const std::size_t klim = 3;  // at most 3 crossings per edge
typedef LayerSearcher<klim> Searcher;
typedef Searcher::Config Config;

struct Setup {
    std::size_t cycle_len;                                      // c
    std::size_t num_cycles;                                     // k >= 2
    std::string run_dir;                                        // output directory of the run
    std::size_t num_tasks;                                      // tasks of a run (e.g. the Slurm array size)
    std::size_t split_depth;                                    // edges of a prefix
    std::function<std::vector<Edge>(std::size_t)> layer_edges;  // new cycle (nm, ..., nm+c-1) + matching
    std::vector<std::size_t> early_prune_checkpoints;           // multicommodity flow check after these edges
};

// name of a run, e.g. c12_3cycles
inline std::string run_name(const Setup& s) {
    return "c" + std::to_string(s.cycle_len) + "_" + std::to_string(s.num_cycles) + "cycles";
}

// searcher settings for adding one layer (Pass 1: crossable new cycle, Pass 2: uncrossed new cycle)
inline Config make_config(const Setup& s, bool pass1) {
    Config config;
    config.cycle_size = s.cycle_len;
    config.layer_edges = s.layer_edges;
    config.early_prune_checkpoints = s.early_prune_checkpoints;
    config.pass1 = pass1;
    return config;
}

// --------------------------------------------------------------------------------------------------------------
// run: one task

class LayerSearch {
public:
    LayerSearch(const Setup& s, const fs::path& dir) : s_(s), dir_(dir), seen_(s.num_cycles) {
        fs::create_directories(dir_);
        hash_file_ = open_append("intermediate.txt");
        complete_file_ = open_append(run_name(s_) + "_complete.jsonl");
    }

    // this task's part of the search
    void run(std::size_t task_index) {
        if (s_.num_cycles == 2) {
            if (task_index == 0) add_next_cycle(nullptr, 1, 0);  // nothing to split: task 0 does it all
        } else {
            Config config = make_config(s_, true);
            config.task_index = task_index;
            config.task_count = s_.num_tasks;
            config.split_depth = s_.split_depth;
            config.output_dir = dir_.string();
            config.intermediate_callback = [this](const Drawing<klim>& d, const IsoCode& code) { intermediate(d, 2, code); };
            Searcher(config).run();
        }
        hash_file_ << "finished\n" << std::flush;
    }

    void print_summary() const {
        for (std::size_t j = 2; j < s_.num_cycles; ++j)
            std::cout << "distinct intermediate drawings with " << j << " cycles: " << seen_[j].size() << "\n";
        if (s_.num_cycles > 2) std::cout << "finishable: " << finishable_.size() << "\n";
        std::cout << "distinct complete drawings: " << complete_seen_.size() << "\n";
    }

private:
    // a restarted (killed) task appends: remove a line cut off by the kill first
    std::ofstream open_append(const std::string& name) {
        const fs::path f = dir_ / name;
        drop_partial_line(f.string());
        return std::ofstream(f, std::ios::app);
    }

    // a reduced intermediate drawing with j cycles; code: canonical code up to the gadget symmetries
    void intermediate(const Drawing<klim>& d, std::size_t j, const IsoCode& code) {
        const std::size_t h = IsoCodeHash()(code);
        if (!seen_[j].insert(h).second) return;  // an isomorphic drawing was already extended in this task
        hash_file_ << j << " " << h << "\n" << std::flush;
        add_next_cycle(&d, j, h);
    }

    // add the next cycle to d (j cycles, nullptr: the base drawing): crossable (Pass 1) if it is a middle cycle,
    // else the uncrossed last cycle (Pass 2)
    void add_next_cycle(const Drawing<klim>* d, std::size_t j, std::size_t hash) {
        const bool last = j + 1 == s_.num_cycles;
        Config config = make_config(s_, !last);
        config.input_drawing = d;
        if (last) {
            config.solution_callback = [this, j, hash](const Drawing<klim>& next) {
                if (j > 1 && finishable_.insert(hash).second) hash_file_ << "finishable " << hash << "\n" << std::flush;
                complete(next);
            };
        } else {
            config.intermediate_callback = [this, j](const Drawing<klim>& next, const IsoCode& code) {
                intermediate(next, j + 1, code);
            };
        }
        Searcher(config).run();
    }

    // a complete drawing: record it unless it is isomorphic to an earlier one of this task
    void complete(const Drawing<klim>& d) {
        if (complete_seen_.insert(canonical_code_free(d)).second)
            complete_file_ << d.serialize_to_json().dump() << "\n" << std::flush;
    }

    const Setup& s_;
    const fs::path dir_;
    std::vector<std::unordered_set<std::size_t>> seen_;  // [j]: hashes of the intermediate drawings with j cycles
    std::unordered_set<std::size_t> finishable_;
    std::unordered_set<IsoCode, IsoCodeHash> complete_seen_;
    std::ofstream hash_file_, complete_file_;
};

inline int run(const Setup& s, std::size_t task_index) {
    if (task_index >= s.num_tasks) throw std::runtime_error("task index must be below the number of tasks");
    LayerSearch search(s, fs::path(s.run_dir) / ("t" + std::to_string(task_index)));
    search.run(task_index);
    std::cout << run_name(s) << ", task " << task_index << "/" << s.num_tasks << "\n";
    search.print_summary();
    return 0;
}

// --------------------------------------------------------------------------------------------------------------
// report: combine the tasks of a run

inline int report(const Setup& s) {
    const fs::path run_dir = s.run_dir;
    const std::string name = run_name(s);
    std::map<std::size_t, std::unordered_set<std::string>> intermediate;  // cycles -> distinct hashes
    std::unordered_set<std::string> finishable;
    std::unordered_set<IsoCode, IsoCodeHash> complete_seen;
    std::ofstream all(run_dir / (name + "_complete.jsonl"));
    std::size_t unfinished = 0;
    const std::size_t tasks = s.num_cycles == 2 ? 1 : s.num_tasks;       // with 2 cycles task 0 does it all
    for (std::size_t t = 0; t < tasks; ++t) {
        const fs::path dir = run_dir / ("t" + std::to_string(t));
        bool finished = false;
        std::ifstream hashes(dir / "intermediate.txt");
        for (std::string what; hashes >> what;) {
            if (what == "finished") { finished = true; continue; }
            std::string h;
            hashes >> h;
            if (what == "finishable") finishable.insert(h);
            else intermediate[std::stoul(what)].insert(h);
        }
        if (!finished) ++unfinished;
        std::ifstream in(dir / (name + "_complete.jsonl"));
        for (std::string line; std::getline(in, line);)
            if (!line.empty() && complete_seen.insert(canonical_code_free(Drawing<klim>(nlohmann::json::parse(line)))).second)
                all << line << "\n";
    }

    std::cout << name << ", " << tasks << " tasks, unfinished: " << unfinished << "\n";
    for (const auto& [j, hashes] : intermediate)
        std::cout << "distinct intermediate drawings with " << j << " cycles: " << hashes.size() << "\n";
    if (s.num_cycles > 2) std::cout << "finishable intermediate drawings with " << s.num_cycles - 1 << " cycles: " << finishable.size() << "\n";
    std::cout << "distinct complete drawings: " << complete_seen.size() << "\n"
              << "complete drawings written to " << (run_dir / (name + "_complete.jsonl")).string() << "\n";
    return 0;
}

inline int main(const Setup& s, int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "";
    if (mode == "run" && argc == 3) return run(s, std::stoul(argv[2]));
    if (mode == "report" && argc == 2) return report(s);
    std::cerr << "usage:\n"
              << "  " << argv[0] << " run <task_index>    (task_index < " << s.num_tasks << ", output in " << s.run_dir << "/t<task_index>)\n"
              << "  " << argv[0] << " report\n";
    return 2;
}

}  // namespace nested_layers

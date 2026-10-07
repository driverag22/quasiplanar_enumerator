#include <array>
#include "iterative_nested_split.h"

// C12 gadget for the split driver (see split_driver.h for usage).
using namespace nested_cycle_build;
const std::size_t klim = 3;

NestedCycleSearcher<klim>::Config make_config() {
    NestedCycleSearcher<klim>::Config config;
    config.cycle_size = 12;
    config.max_cycles = 2;
    config.local_edges_builder = [](std::size_t nm) -> std::vector<Edge> {
        return {
             // one matching
             {0,nm+9},
             // outer uncrossed cycle
             {nm+9 ,nm+10,klim},
             {nm+10,nm+11,klim},
             {nm+11,nm   ,klim},
             {nm   ,nm+1 ,klim},
             {nm+1 ,nm+2 ,klim},
             {nm+2 ,nm+3 ,klim},
             {nm+3 ,nm+4 ,klim},
             {nm+4 ,nm+5 ,klim},
             {nm+5 ,nm+6 ,klim},
             {nm+6 ,nm+7 ,klim},
             {nm+7 ,nm+8 ,klim},
             {nm+8 ,nm+9 ,klim},
             // rest of matching
             {2,nm+11},{4,nm+1},{6,nm+3},{8,nm+5},{10,nm+7},
             {1,nm+4},{3,nm+6},{5,nm+8},{7,nm+10},{9,nm},{11,nm+2},
        };
    };
    // optional extra edge classes (environment, several classes separated by '/', each "d:p[:s]" or "d,p[,s]"):
    //   EXTRA_BRAID: old vertex i -> new vertex i+d (mod 12) for all i with i mod s == p (in every layer)
    //   EXTRA_CHORD: chord new vertex i -> new vertex i+d (mod 12) for all i with i mod s == p, Pass 1 only (pairs
    //                within D1 / D3 count as present in the saturation test)
    // step s defaults to 2
    auto parse_classes = [](const char* name, const std::string& spec) {
        std::vector<std::array<std::size_t, 3>> classes; // {d, p, s}
        std::size_t begin = 0;
        for (;;) {
            const std::size_t end = spec.find('/', begin);
            const std::string item = spec.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
            std::vector<std::size_t> v;
            std::size_t b = 0;
            for (;;) {
                const std::size_t sep = item.find_first_of(",:", b);
                const std::string tok = item.substr(b, sep == std::string::npos ? std::string::npos : sep - b);
                if (tok.empty() || tok.find_first_not_of("0123456789") != std::string::npos)
                    throw std::runtime_error(std::string(name) + ": classes must be d:p[:s], got " + spec);
                v.push_back(std::stoul(tok));
                if (sep == std::string::npos) break;
                b = sep + 1;
            }
            if (v.size() < 2 || v.size() > 3 || (v.size() == 3 && (v[2] == 0 || 12 % v[2] != 0)))
                throw std::runtime_error(std::string(name) + ": classes must be d:p[:s] with s dividing 12, got " + spec);
            classes.push_back({v[0], v[1], v.size() == 3 ? v[2] : 2});
            if (end == std::string::npos) break;
            begin = end + 1;
        }
        return classes;
    };
    if (const char* x = std::getenv("EXTRA_BRAID")) {
        const auto classes = parse_classes("EXTRA_BRAID", x);
        auto base_builder = config.local_edges_builder;
        config.local_edges_builder = [base_builder, classes](std::size_t nm) {
            auto es = base_builder(nm);
            std::set<std::pair<std::size_t, std::size_t>> have;
            for (const auto& e : es) have.insert({e[0], e[1]});
            for (const auto& [d, par, step] : classes)
                for (std::size_t i = par % step; i < 12; i += step)
                    if (have.insert({i, nm + (i + d) % 12}).second) es.push_back({i, nm + (i + d) % 12});
            return es;
        };
    }
    if (const char* x = std::getenv("EXTRA_CHORD")) {
        const auto classes = parse_classes("EXTRA_CHORD", x);
        config.pass1_extra_edges_builder = [classes](std::size_t nm) {
            std::set<std::pair<std::size_t, std::size_t>> chords;
            for (const auto& [k, par, step] : classes)
                for (std::size_t i = par % step; i < 12; i += step)
                    if (k % 12 > 1 && k % 12 < 11) // not a cycle edge or loop
                        chords.insert({nm + std::min(i, (i + k) % 12), nm + std::max(i, (i + k) % 12)});
            std::vector<Edge> es;
            for (const auto& [a, b] : chords) es.push_back({a, b});
            return es;
        };
    }
    config.early_prune_checkpoints = {16,18,21};
    config.enable_mcf_pruning = true;
    config.mcf_log_progress = false;
    config.verbose = false;
    return config;
}

#include "split_driver.h"

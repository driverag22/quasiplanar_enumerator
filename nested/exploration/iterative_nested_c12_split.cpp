#include <array>
#include "iterative_nested_split.h"

// C12 gadget for the split driver (see split_driver.h for usage).
using namespace nested_cycle_build;
const std::size_t klim = 3;
using Config = NestedCycleSearcher<klim>::Config;

// extra edge classes, several classes separated by '/', each "d:p[:s]" or "d,p[,s]" (step s defaults to 2):
//   braid: old vertex i -> new vertex i+d (mod 12), chord: new vertex i -> new vertex i+d (mod 12),
//   for all i with i mod s == p
static std::vector<std::array<std::size_t, 3>> parse_classes(const std::string& name, const std::string& spec) {
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
                throw std::runtime_error(name + ": classes must be d:p[:s], got " + spec);
            v.push_back(std::stoul(tok));
            if (sep == std::string::npos) break;
            b = sep + 1;
        }
        if (v.size() < 2 || v.size() > 3 || (v.size() == 3 && (v[2] == 0 || 12 % v[2] != 0)))
            throw std::runtime_error(name + ": classes must be d:p[:s] with s dividing 12, got " + spec);
        classes.push_back({v[0], v[1], v.size() == 3 ? v[2] : 2});
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return classes;
}

// braid edges of the classes that are not among the edges of base (yet)
static std::vector<Edge> braid_edges(const std::vector<std::array<std::size_t, 3>>& classes, std::size_t nm,
        const std::vector<Edge>& base) {
    std::set<std::pair<std::size_t, std::size_t>> have;
    for (const auto& e : base) have.insert({e[0], e[1]});
    std::vector<Edge> es;
    for (const auto& [d, par, step] : classes)
        for (std::size_t i = par % step; i < 12; i += step)
            if (have.insert({i, nm + (i + d) % 12}).second) es.push_back({i, nm + (i + d) % 12});
    return es;
}

// chords of the new cycle
static std::vector<Edge> chord_edges(const std::vector<std::array<std::size_t, 3>>& classes, std::size_t nm) {
    std::set<std::pair<std::size_t, std::size_t>> chords;
    for (const auto& [k, par, step] : classes)
        for (std::size_t i = par % step; i < 12; i += step)
            if (k % 12 > 1 && k % 12 < 11) // not a cycle edge or loop
                chords.insert({nm + std::min(i, (i + k) % 12), nm + std::max(i, (i + k) % 12)});
    std::vector<Edge> es;
    for (const auto& [a, b] : chords) es.push_back({a, b});
    return es;
}

Config make_config() {
    Config config;
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
    // optional extra edge classes (environment, syntax as in parse_classes):
    //   EXTRA_BRAID: braid edges in every layer
    //   EXTRA_CHORD: chords of the new cycle, Pass 1 only (pairs within D1 / D3 count as present in the saturation test)
    //   EXTRA_BRAID_P1 / EXTRA_BRAID_P2: like EXTRA_BRAID, but only in Pass 1 (D1-D2) / only in Pass 2 (D2-D3),
    //                for braid edges in every other layer only
    // (per layer and per cycle with more cycles: SAT_LAYERS / SAT_CHORDS / SAT_JUMPS in split_driver.h, via configure_layer)
    if (const char* x = std::getenv("EXTRA_BRAID")) {
        const auto classes = parse_classes("EXTRA_BRAID", x);
        auto base_builder = config.local_edges_builder;
        config.local_edges_builder = [base_builder, classes](std::size_t nm) {
            auto es = base_builder(nm);
            const auto extra = braid_edges(classes, nm, es);
            es.insert(es.end(), extra.begin(), extra.end());
            return es;
        };
    }
    if (const char* x = std::getenv("EXTRA_BRAID_P2")) {
        const auto classes = parse_classes("EXTRA_BRAID_P2", x);
        config.pass2_extra_edges_builder = [base_builder = config.local_edges_builder, classes](std::size_t nm) {
            return braid_edges(classes, nm, base_builder(nm));
        };
    }
    if (const char* x = std::getenv("EXTRA_CHORD")) {
        const auto classes = parse_classes("EXTRA_CHORD", x);
        config.pass1_extra_edges_builder = [classes](std::size_t nm) { return chord_edges(classes, nm); };
    }
    if (const char* x = std::getenv("EXTRA_BRAID_P1")) { // braid edges first, then the chords (if any)
        const auto classes = parse_classes("EXTRA_BRAID_P1", x);
        config.pass1_extra_edges_builder = [base_builder = config.local_edges_builder, chords = config.pass1_extra_edges_builder,
                classes](std::size_t nm) {
            auto es = braid_edges(classes, nm, base_builder(nm));
            if (chords) { const auto c = chords(nm); es.insert(es.end(), c.begin(), c.end()); }
            return es;
        };
    }
    config.early_prune_checkpoints = {16,18,21};
    config.enable_mcf_pruning = true;
    config.mcf_log_progress = false;
    config.verbose = false;
    return config;
}

// jump edges: vertex i of the cycle below the old one -> new vertex i+d (mod 12), written as {12 + i, nm + (i+d) % 12};
// the searcher's active_cycle then lists that cycle after the old one (Config::jump_edges_builder, SAT_JUMPS)
static std::vector<Edge> jump_edges(const std::vector<std::array<std::size_t, 3>>& classes, std::size_t nm) {
    std::set<std::pair<std::size_t, std::size_t>> jumps;
    for (const auto& [d, par, step] : classes)
        for (std::size_t i = par % step; i < 12; i += step) jumps.insert({12 + i, nm + (i + d) % 12});
    std::vector<Edge> es;
    for (const auto& [a, b] : jumps) es.push_back({a, b});
    return es;
}

// local edges with other base braids than 9:0/3:1 (classes in parse_classes syntax): the base edge from the smallest
// old vertex, the new cycle starting at its end, then the other base edges in class order (as make_config for 9:0/3:1)
static std::function<std::vector<Edge>(std::size_t)> base_builder_for(const std::vector<std::array<std::size_t, 3>>& classes) {
    return [classes](std::size_t nm) {
        std::vector<Edge> braids = braid_edges(classes, nm, {});
        if (braids.empty()) throw std::runtime_error("base braid classes without edges");
        std::stable_sort(braids.begin(), braids.end(), [](const Edge& a, const Edge& b) { return a[0] == 0 && b[0] != 0; });
        const std::size_t t0 = braids[0][1] - nm;
        std::vector<Edge> es{braids[0]};
        for (std::size_t j = 0; j < 12; ++j) es.push_back({nm + (t0 + j) % 12, nm + (t0 + j + 1) % 12, klim});
        es.insert(es.end(), braids.begin() + 1, braids.end());
        return es;
    };
}

// extra edges of one layer (for the saturation test with SAT_LAYERS / SAT_CHORDS / SAT_JUMPS, see split_driver.h):
// braid classes between the old and the new cycle, chord classes of the new cycle and jump classes from the cycle
// below the old one to the new cycle ("-" or empty = none); pass2: the new cycle is the uncrossed last one (chords not
// allowed). Replaces the extras and the jump edges of config. A braid entry "<extras>|<base>" (e.g. "1:1|3:0/3:1" or
// "-|3:0/3:1") also replaces the base braids 9:0/3:1 of the layer.
void configure_layer(Config& config, const std::string& braid_spec, const std::string& chord, const std::string& jump, bool pass2) {
    const std::size_t bar = braid_spec.find('|');
    const std::string braid = braid_spec.substr(0, bar);
    if (bar != std::string::npos && braid_spec.substr(bar + 1) != "-")
        config.local_edges_builder = base_builder_for(parse_classes("SAT_LAYERS base", braid_spec.substr(bar + 1)));
    const bool has_braid = !braid.empty() && braid != "-", has_chord = !chord.empty() && chord != "-";
    const bool has_jump = !jump.empty() && jump != "-";
    const auto bc = has_braid ? parse_classes("SAT_LAYERS", braid) : std::vector<std::array<std::size_t, 3>>();
    const auto cc = has_chord ? parse_classes("SAT_CHORDS", chord) : std::vector<std::array<std::size_t, 3>>();
    const auto jc = has_jump ? parse_classes("SAT_JUMPS", jump) : std::vector<std::array<std::size_t, 3>>();
    config.jump_edges_builder = nullptr;
    if (has_jump) config.jump_edges_builder = [jc](std::size_t nm) { return jump_edges(jc, nm); };
    if (pass2) {
        if (has_chord) throw std::runtime_error("configure_layer: no chords on the uncrossed last cycle");
        config.pass2_extra_edges_builder = nullptr;
        if (has_braid)
            config.pass2_extra_edges_builder = [base_builder = config.local_edges_builder, bc](std::size_t nm) {
                return braid_edges(bc, nm, base_builder(nm));
            };
        return;
    }
    config.pass1_extra_edges_builder = nullptr;
    if (has_braid || has_chord)
        config.pass1_extra_edges_builder = [base_builder = config.local_edges_builder, bc, cc](std::size_t nm) {
            auto es = braid_edges(bc, nm, base_builder(nm));
            const auto c = chord_edges(cc, nm);
            es.insert(es.end(), c.begin(), c.end());
            return es;
        };
}

#include "split_driver.h"

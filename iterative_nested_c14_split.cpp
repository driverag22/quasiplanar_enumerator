#include "iterative_nested_split.cpp"

// C14 gadget for the split driver (see split_driver.h for usage).
using namespace nested_cycle_build;
const std::size_t klim = 3;

NestedCycleSearcher<klim>::Config make_config() {
    NestedCycleSearcher<klim>::Config config;
    config.cycle_size = 14;
    config.max_cycles = 2;
    config.local_edges_builder = [](std::size_t nm) -> std::vector<Edge> {
        return {
            // one matching edge to prev layer
            {0,nm+11}, 
            // outer cycle
            {nm+11,nm+12,klim},
            {nm+12,nm+13,klim},
            {nm+13,nm   ,klim},
            {nm   ,nm+1 ,klim},
            {nm+1 ,nm+2 ,klim},
            {nm+2 ,nm+3 ,klim},
            {nm+3 ,nm+4 ,klim},
            {nm+4 ,nm+5 ,klim},
            {nm+5 ,nm+6 ,klim},
            {nm+6 ,nm+7 ,klim},
            {nm+7 ,nm+8 ,klim},
            {nm+8 ,nm+9 ,klim},
            {nm+9 ,nm+10,klim},
            {nm+10,nm+11,klim},
            // rest of matching
            // EVEN: (i, nm + ( (i+11) % 14)) == (i, nm + ( (i-3) % 14))
            // ODD: (i, nm + ( (i+3) % 14)) == (i, nm + ( (i-11) % 14))
            {2,nm+13},{4,nm+1},{6,nm+3},{8,nm+5},{10,nm+7},{12,nm+9},
            {1,nm+4},{3,nm+6},{5,nm+8},{7,nm+10},{9,nm+12},{11,nm+0},{13,nm+2}
        };
    };
    // at 25%, 50% and 75%
    config.early_prune_checkpoints = {19,21,24};
    config.enable_mcf_pruning = true;
    config.mcf_log_progress = false;
    config.verbose = false;
    return config;
}

#include "split_driver.h"

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
    config.early_prune_checkpoints = {16,18,21};
    config.enable_mcf_pruning = true;
    config.mcf_log_progress = false;
    config.verbose = false;
    return config;
}

#include "split_driver.h"

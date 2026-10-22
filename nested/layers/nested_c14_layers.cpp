// All 3-plane drawings of the nested C14 graph G_{14,k} (D1 and Dk uncrossed), see nested_layers.h.
//
// Build:  g++ -O3 -std=c++17 -Iinclude nested/layers/nested_c14_layers.cpp -o .bin/nested_c14_layers
// Run:    .bin/nested_c14_layers run <task_index>     for task_index = 0, ..., num_tasks-1 (in parallel)
//         .bin/nested_c14_layers report               after all tasks finished

#include "nested_layers.h"

using namespace nested_layers;

int main(int argc, char** argv) {
    Setup setup;
    setup.cycle_len = 14;
    setup.num_cycles = 3;     // k: 2 = G_{14,2}, 3 = G_{14,3}, ...
    setup.run_dir = "../quasiDrawings/nestedCycles/nested_c14/";  // tasks in run_dir/t<task_index>
    setup.num_tasks = 100;
    setup.split_depth = 12;
    setup.layer_edges = [](std::size_t nm) -> std::vector<Edge> {
        return {
            // one matching edge
            {0, nm + 11},
            // new cycle (uncrossable in Pass 2)
            {nm + 11, nm + 12, klim},
            {nm + 12, nm + 13, klim},
            {nm + 13, nm     , klim},
            {nm     , nm + 1 , klim},
            {nm + 1 , nm + 2 , klim},
            {nm + 2 , nm + 3 , klim},
            {nm + 3 , nm + 4 , klim},
            {nm + 4 , nm + 5 , klim},
            {nm + 5 , nm + 6 , klim},
            {nm + 6 , nm + 7 , klim},
            {nm + 7 , nm + 8 , klim},
            {nm + 8 , nm + 9 , klim},
            {nm + 9 , nm + 10, klim},
            {nm + 10, nm + 11, klim},
            // rest of the matching
            // EVEN: (i, nm + (i-3) % 14), ODD: (i, nm + (i+3) % 14)
            {2, nm + 13}, {4, nm + 1}, {6, nm + 3}, {8, nm + 5}, {10, nm + 7}, {12, nm + 9},
            {1, nm + 4}, {3, nm + 6}, {5, nm + 8}, {7, nm + 10}, {9, nm + 12}, {11, nm}, {13, nm + 2},
        };
    };
    setup.early_prune_checkpoints = {19, 21, 24};  // flow check after the edges at these positions (0-based, of 28)
    return nested_layers::main(setup, argc, argv);
}

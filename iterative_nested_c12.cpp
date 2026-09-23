#include "iterative_nested.cpp"

int main() {
    const std::size_t klim = 3;
    using namespace nested_cycle_build;

    NestedCycleSearcher<klim>::Config config;
    config.cycle_size = 12;
    config.max_cycles = 2;
    config.output_dir = "../quasiDrawings/nested_c12";
    config.local_edges_builder = [](std::size_t nm) -> std::vector<nested_cycle_build::Edge> {
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
             // EVEN: (i, nm + ( (i+9) % 12)) == (i, nm + ( (i-3) % 12))
             // ODD: (i, nm + ( (i+3) % 12)) == (i, nm + ( (i-9) % 12))
             {2,nm+11},{4,nm+1},{6,nm+3},{8,nm+5},{10,nm+7},
             {1,nm+4},{3,nm+6},{5,nm+8},{7,nm+10},{9,nm},{11,nm+2},
        };
    };
    // at 25%, 50% and 75%
    config.early_prune_checkpoints = {15,18,21}; // indices of local_edges_builder to check early prune

    // Instantiate and execute search
    NestedCycleSearcher<3> searcher(config);

    std::cout << "Starting 3-Planar Nested Cycle Search with C12..." << std::endl;
    auto result = searcher.run();

    std::cout << "\n==================================================" << std::endl;
    std::cout << "C12 Search Complete!" << std::endl;
    std::cout << "Total Subdrawings Processed: " << result.total_processed << std::endl;
    std::cout << "Total Intermediate Subdrawings: " << result.solutions.size() << std::endl;
    std::cout << "Total Terminal Solutions Found: " << result.full_solution_count << std::endl;
    std::cout << "Total Early Pruned Backtracks: "  << result.pruned_early_count << std::endl;
    std::cout << "Total Discarded Isomorphisms: " << result.discarded_count << std::endl;
    std::cout << "==================================================" << std::endl;

    return 0;
}


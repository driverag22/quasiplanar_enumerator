#include "iterative_nested.cpp"

int main() {
    using namespace nested_cycle_build;

    // 1. Configure for C12
    NestedCycleSearcher<3>::Config config14;
    config14.cycle_size = 12;
    config14.max_cycles = 2;
    config14.output_dir = "../quasiDrawings/nested_c12";
    config14.local_edges_builder = [](std::size_t nm) -> std::vector<nested_cycle_build::Edge> {
        return {
            {0,nm+11}, // one matching
            // outer cycle
            {nm+11,nm+12,3},
            {nm+12,nm+13,3},
            {nm+13,nm   ,3},
            {nm   ,nm+1 ,3},
            {nm+1 ,nm+2 ,3},
            {nm+2 ,nm+3 ,3},
            {nm+3 ,nm+4 ,3},
            {nm+4 ,nm+5 ,3},
            {nm+5 ,nm+6 ,3},
            {nm+6 ,nm+7 ,3},
            {nm+7 ,nm+8 ,3},
            {nm+8 ,nm+9 ,3},
            {nm+9 ,nm+10,3},
            {nm+10,nm+11,3},
            // rest of matching
            // EVEN: (i, nm + ( (i+11) % 14)) == (i, nm + ( (i-3) % 14))
            // ODD: (i, nm + ( (i+3) % 14)) == (i, nm + ( (i-11) % 14))
            {2,nm+13},{4,nm+1},{6,nm+3},{8,nm+5},{10,nm+7},{12,nm+9},
            {1,nm+4},{3,nm+6},{5,nm+8},{7,nm+10},{9,nm+11},{11,nm+0},{13,nm+2}
        };
    };

    // Instantiate and execute search
    NestedCycleSearcher<3> searcher(config14);

    std::cout << "Starting 3-Planar Nested Cycle Search with C14..." << std::endl;
    auto result = searcher.run();

    std::cout << "\n==================================================" << std::endl;
    std::cout << "C14 Search Complete!" << std::endl;
    std::cout << "Total Subdrawings Processed: " << result.total_processed << std::endl;
    std::cout << "Total Intermediate Subdrawings: " << result.solutions.size() << std::endl;
    std::cout << "Total Terminal Solutions Found: " << result.full_solution_count << std::endl;
    std::cout << "Total Discarded Isomorphisms: " << result.discarded_count << std::endl;
    std::cout << "==================================================" << std::endl;

    return 0;
}


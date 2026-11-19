#include <toggle/toggle.h>

#ifndef GRAPH_SYMMETRIZED
#define GRAPH_SYMMETRIZED true
#endif

#ifndef GRAPH_WEIGHT_TYPE
#define GRAPH_WEIGHT_TYPE void
#endif

int main(int argc, char** argv) {
    std::string file_path = argv[1];
    uint32_t num_rounds = (argc <= 2) ? 3 : std::atoi(argv[2]);
    std::string source_path = (argc <= 3) ? "disabled" : argv[3];
    std::string dump_path = (argc <= 4) ? "disabled" : argv[4];

    toggle::Timer t; double tt = 0, ttt = 0;
    t.start();
    auto G = toggle::read_graph<GRAPH_SYMMETRIZED, GRAPH_WEIGHT_TYPE>(file_path);
    tt = t.stop();
    std::string graph_name = toggle::extract_graph_name(file_path);
    std::cout << "==================================================================\n";
    std::cout << std::right << std::setw(66) << ("Graph: " + graph_name) << "\n";
    std::cout << "Load: " << tt << "  Dump: " << dump_path << "\n";
    std::cout << "Threads: " << toggle::num_workers() << "  Rounds: " << num_rounds << "\n";

    num_rounds = std::min<uint32_t>(G.n(), num_rounds);
    auto sources = toggle::get_sources(source_path, graph_name, num_rounds);
    for (uint32_t round = 0; round < num_rounds; round++) {
        auto source = sources[round];
        {
            t.start();
            auto result = toggle::dijkstra(G, source);
            tt = t.stop();
            std::cout << "    Round " << round + 1 << "  source: " << source << "  Warmup: "  << std::setprecision(2) << tt << std::setprecision(6);
        }
        {
            t.start();
            auto result = toggle::dijkstra(G, source);
            tt = t.stop();
            std::cout << "  time = " << tt << " sec\n";
            ttt += tt;
            if (round == num_rounds - 1) {
                toggle::process_result(ttt / num_rounds, result, ".", graph_name, "ToT", dump_path);
            }
        }
    }

    return 0;
}

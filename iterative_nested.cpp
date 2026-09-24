#ifndef NESTED_CYCLE_SEARCH_HPP
#define NESTED_CYCLE_SEARCH_HPP

#include "hds_kplanar.h"
#include "iso.h"
#include <iostream>
#include <vector>
#include <deque>
#include <fstream>
#include <cassert>
#include <nlohmann/json.hpp>

// BGL Flow
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

namespace nested_cycle_build {

    // Max-Flow Dual Network Structure
    struct DualNetwork {
        typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;

        typedef boost::adjacency_list<
            boost::vecS, boost::vecS, boost::directedS,
            boost::no_property,
            boost::property<boost::edge_capacity_t, int,
            boost::property<boost::edge_residual_capacity_t, int,
            boost::property<boost::edge_reverse_t, Traits::edge_descriptor>>>
                > Graph;

        Graph G;
        DualNetwork(std::size_t num_nodes) : G(num_nodes) {}

        void add_edge(int from, int to, int capacity) {
            auto c_map = boost::get(boost::edge_capacity, G);
            auto r_map = boost::get(boost::edge_reverse, G);
            const auto e = boost::add_edge(from, to, G).first;
            const auto rev_e = boost::add_edge(to, from, G).first;
            c_map[e] = capacity;
            c_map[rev_e] = 0; // residual reverse capacity
            r_map[e] = rev_e;
            r_map[rev_e] = e;
        }

        int flow(int source, int sink) {return boost::push_relabel_max_flow(G, source, sink);}
    };

    typedef std::vector<std::size_t> Edge;
    typedef std::vector<Edge> Edges;

    template <std::size_t klim>
    inline void check_and_terminate_if_invalid(const Drawing<klim>& d, const std::string& filename = "../quasiDrawings/failExample.graphml") {
        // 1. Check for self-loops (u == v)
        for (const auto& edge : d.edges) {
            if (edge.u == edge.v) {
                std::cout << "\n[INVALID DRAWING] Self-loop detected at vertex " << edge.u << "!" << std::endl;
                std::ofstream of_graphml(filename);
                d.graphml_output(of_graphml);
                of_graphml.close();
                std::cerr << "Saved invalid drawing to " << filename << ". Terminating execution." << std::endl;
                std::exit(1);
            }
        }

        // 2. Check for parallel edges (duplicate {u, v} pairs)
        std::set<std::pair<std::size_t, std::size_t>> seen_edges;
        for (const auto& edge : d.edges) {
            std::size_t u = std::min(edge.u, edge.v);
            std::size_t v = std::max(edge.u, edge.v);
            if (seen_edges.count({u, v}) > 0) {
                std::cout << "\n[INVALID DRAWING] Parallel edge detected between vertices " << u << " and " << v << "!" << std::endl;
                std::ofstream of_graphml(filename);
                d.graphml_output(of_graphml);
                of_graphml.close();
                std::cerr << "Saved invalid drawing to " << filename << ". Terminating execution." << std::endl;
                std::exit(1);
            }
            seen_edges.insert({u, v});
        }
    }

    // Returns a halfedge pointing to the desired source vertex, and the number 
    // of crossings the edge should have in the reduced drawing.
    //
    // Used when a halfedge has an irrelevant endpoint:
    //   vact[v] == vact.size()   iff   v is irrelevant
    template <std::size_t klim, typename F, typename A, typename V>
        inline std::pair<const HdsHalfedge*, std::size_t> prev_active(
                const HdsHalfedge* e,
                const F& face,
                const A& active,
                const V& vact) 
        {
            if (active[face[e->label]] >= 0) {
                std::size_t c = (active[face[e->twin->label]] < 0 ? klim : e->edge->ncr);
                do e = e->prev; while (vact[e->vertex->label] == vact.size());
                return std::make_pair(e, c);
            }
            e = e->twin;
            if (active[face[e->label]] < 0) throw std::runtime_error("no active face");
            while (vact[e->vertex->label] == vact.size()) e = e->next;
            return std::make_pair(e, (active[face[e->twin->label]] < 0 ? klim : e->edge->ncr));
        }

    template <std::size_t klim = 3>
        class NestedCycleSearcher {
            public:
                // Queue node wrapping drawing and explicit depth level
                struct SearchNode {
                    Drawing<klim> drawing;
                    std::size_t depth;
                };
                // Configuration options for running a search instance
                struct Config {
                    std::size_t cycle_size = 12;
                    std::size_t max_cycles = 2;
                    std::string output_dir = "./output_drawings";
                    bool export_files = true;
                    bool verbose = true;
                    bool enable_early_pruning = true;

                    // custom edge index checkpoints where early flow checks run.
                    std::vector<std::size_t> early_prune_checkpoints;

                    // Generator function for local edges added during extension step.
                    // Takes parent vertex offset (nm) and cycle size, returns vector of edges {u, v, [capacity]}.
                    std::function<std::vector<Edge>(std::size_t nm)> local_edges_builder;
                };

                // Search summary stats
                struct SearchResult {
                    std::vector<SearchNode> solutions;
                    std::size_t discarded_count = 0;
                    std::size_t full_solution_count = 0;
                    std::size_t total_processed = 0;
                    std::size_t pruned_early_count = 0;
                };

                // Constructors
                NestedCycleSearcher() = default;
                explicit NestedCycleSearcher(Config config) : config_(std::move(config)) {}


                Drawing<klim> create_base_drawing() const {
                    std::size_t cycle_size = config_.cycle_size;
                    // 12 vertices for C12_1 + 1 dummy star vertex = 15 vertices
                    Drawing<klim> d(cycle_size + 1);
                    std::vector<HdsHalfedge*> cycle(cycle_size, nullptr);

                    // build uncrossable C12 cycle by setting ncr=klim=3
                    cycle[1] = d.add_first_edge(0, 1, klim);
                    for (std::size_t i = 2; i < cycle_size; ++i) {
                        cycle[i] = d.add_edge(HdsPath({cycle[i - 1], nullptr}), i, klim);
                    }
                    cycle[0] = d.add_edge(HdsPath({cycle[cycle_size - 1], cycle[1]->twin}), 0, klim);

                    // block one side of the cycle using an uncrossable star at dummy vertex 12
                    std::size_t star_center = cycle_size;
                    auto e = d.add_edge(HdsPath({cycle[0], nullptr}), star_center, klim);
                    for (std::size_t i = cycle_size - 1; i > 0; --i) {
                        e = d.add_edge(HdsPath({cycle[i], e}), star_center, klim);
                    }

                    return d;
                }

                bool check_remaining_edges_reachability(
                        const Drawing<klim>& d,
                        const std::vector<Edge>& local_edges,
                        std::size_t next_edge_index
                        ) const {
                    // map halfedges to face indices
                    std::vector<int> face(d.halfedges.size(), -1);
                    std::size_t num_faces = 0;
                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                        if (face[i->label] != -1) continue;
                        const HdsHalfedge* j = &*i;
                        do { // walk around face
                            face[j->label] = static_cast<int>(num_faces);
                            j = j->next;
                        } while (j->label != i->label);
                        ++num_faces;
                    }
                    if (num_faces == 0) return false;

                    // dual adjacency list filtered by available crossing capacity (capacity > 0)
                    std::vector<std::vector<int>> dual_adj(num_faces);
                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                        int remaining_capacity = static_cast<int>(klim) - static_cast<int>(i->edge->ncr);
                        if (remaining_capacity > 0) {
                            int f1 = face[i->label];
                            int f2 = face[i->twin->label];
                            if (f1 >= 0 && f2 >= 0 && f1 != f2) {
                                dual_adj[f1].push_back(f2);
                            }
                        }
                    }

                    // lambda to retrieve all incident face IDs for a given vertex label
                    auto get_incident_faces = [&](std::size_t v_label) -> std::vector<int> {
                        std::vector<int> faces;
                        if (v_label >= d.vertices.size()) return faces;
                        // halfedge pointing to v_label
                        auto start_h = d.vertices[v_label].halfedge;
                        if (!start_h) return faces;

                        auto curr_h = start_h;
                        do {
                            int f = face[curr_h->label];
                            // make sure face is valid and not already in faces
                            if (f >= 0) if (std::find(faces.begin(), faces.end(), f) == faces.end())
                                    faces.push_back(f);
                            curr_h = curr_h->next->twin; // walk around halfedges incident to vertex
                        } while (curr_h != start_h);
                        return faces;
                    };

                    // test each remaining unplaced braid edge individually using BFS
                    // we use the same "visited" array for each edge, by using a different 
                    // token (rem_idx) each time
                    std::vector<int> visited_token(num_faces, -1);
                    std::queue<int> Q;
                    for (std::size_t rem_idx = next_edge_index; rem_idx < local_edges.size(); ++rem_idx) {
                        std::size_t u = local_edges[rem_idx][0];
                        std::size_t v = local_edges[rem_idx][1];
                        std::vector<int> u_faces = get_incident_faces(u);
                        std::vector<int> v_faces = get_incident_faces(v);

                        // If an endpoint has no incident face, routing is impossible
                        if (u_faces.empty() || v_faces.empty()) return false;

                        // quick target lookup mask for sink faces (v_faces)
                        std::vector<bool> is_target(num_faces, false);
                        for (int f_v : v_faces) 
                            is_target[f_v] = true;

                        // initialize queue with source
                        while (!Q.empty()) Q.pop();
                        int current_marker = static_cast<int>(rem_idx);
                        bool shares_face = false;

                        // set starting faces of BFS
                        for (int f_u : u_faces) {
                            if (is_target[f_u]) {
                                shares_face = true;
                                break;
                            }
                            visited_token[f_u] = current_marker;
                            Q.push(f_u);
                        }

                        // Vertices u and v already share a common face -> reachable with 0 crossings
                        if (shares_face) continue;

                        bool reachable = false;
                        while (!Q.empty()) {
                            int curr = Q.front();
                            Q.pop();

                            for (int neighbor : dual_adj[curr]) {
                                if (is_target[neighbor]) {
                                    reachable = true;
                                    break;
                                }
                                if (visited_token[neighbor] != current_marker) {
                                    visited_token[neighbor] = current_marker;
                                    Q.push(neighbor);
                                }
                            }
                            if (reachable) break;
                        }
                        if (!reachable) return false;
                    }

                    return true;
                }


                SearchResult run() {
                    if (!config_.local_edges_builder)
                        throw std::runtime_error("NestedCycleSearcher Error: local_edges_builder callback is not set!");

                    if (config_.export_files)
                        std::filesystem::create_directories(config_.output_dir);


                    SearchResult result;
                    result.solutions.push_back({create_base_drawing(), 0});

                    // we process the solutions in order, this is the current one we process
                    std::size_t solcount = 0;

                    while (solcount < result.solutions.size()) {
                        Drawing<klim> d_parent = result.solutions[solcount].drawing;
                        const std::size_t current_depth = result.solutions[solcount].depth;

                        if (config_.verbose) {
                            std::cout << "\n--- Processing Partial Drawing #" << solcount 
                                << " [Cycle Level " << current_depth << "]"
                                << " (Queue Size: " << result.solutions.size() << ") ---" << std::endl;
                        }

                        if (current_depth >= config_.max_cycles) {
                            if (config_.verbose) {
                                std::cout << "--> Reached MAX_CYCLES limit (" << config_.max_cycles 
                                    << "). Skipping further cycle extension for Drawing #" << solcount << std::endl;
                            }
                            ++solcount;
                            continue;
                        }

                        bool is_extensible = false;
                        // Two-Pass Loop
                        // Pass 0 (constrained == 0): Unconstrained extension & subdrawing extraction
                        // Pass 1 (constrained == 1): Constrained (uncrossable) decagon completion
                        for (int constrained = 0; constrained < 2; ++constrained) {

                            // Pass 2 (constrained == 1) only runs if Pass 1 proved extensibility
                            if (constrained == 1 && !is_extensible) {
                                if (config_.verbose) {
                                    std::cout << "--> Skipping Constrained Pass 2 (Drawing #" << solcount 
                                        << " is not extensible)" << std::endl;
                                }
                                continue;
                            }

                            if (config_.verbose) {
                                std::cout << "\n>>> Starting Pass " << (constrained == 0 ? "1 (Unconstrained)" : "2 (Constrained)")
                                    << " for Drawing #" << solcount << " <<<" << std::endl;
                            }

                            Drawing<klim> d = d_parent;
                            std::size_t nm = d.vertices.size();
                            d.add_vertices(config_.cycle_size);
                            std::vector<Edge> local_edges = config_.local_edges_builder(nm);
                            auto e = local_edges.begin();
                            uint64_t total_local_iterations = 0;

                            for (;;) {
                                ++total_local_iterations;
                                std::size_t u = (*e)[0];
                                std::size_t v = (*e)[1];

                                // In Pass 2 (constrained == 1), newly added C12 cycle edges are uncrossable (pcr = klim). 
                                // In Pass 1 (constrained == 0), edges are unconstrained (pcr = 0).
                                int pcr = (constrained == 1 && e->size() == 3) ? (*e)[2] : 0;
                                HdsPath p = d.first_path(u, v, pcr);

                                if (p.empty()) {
BACKUP:
                                    do {
                                        if (e == local_edges.begin()) {
                                            goto FINISH_PASS;
                                        }
                                        --e;
                                        u = (*e)[0];
                                        v = (*e)[1];
                                        pcr = (constrained == 1 && e->size() == 3) ? (*e)[2] : 0;
                                        p = d.edges.back().built;
                                        d.remove_edge();

                                        if (d.next_path(p, v, pcr)) {
                                            break;
                                        }
                                    } while (true);
                                }
                                d.add_edge(p, v, pcr);

                                // reachability check for remaining unplaced braid edges
                                std::size_t current_edge_idx = static_cast<std::size_t>(e - local_edges.begin());
                                if (config_.enable_early_pruning && (current_edge_idx + 1 < local_edges.size())) {
                                    bool is_checkpoint = false;
                                    is_checkpoint = (
                                        std::find(config_.early_prune_checkpoints.begin(),config_.early_prune_checkpoints.end(),current_edge_idx) 
                                        != config_.early_prune_checkpoints.end());

                                    if (is_checkpoint) {
                                        while (!check_remaining_edges_reachability(d, local_edges, current_edge_idx + 1)) {
                                            // std::string filename2 = "../quasiDrawings/failExample.graphml";
                                            // std::ofstream of_graphml(filename2);
                                            // d.graphml_output(of_graphml);
                                            // of_graphml.close();
                                            // if (config_.verbose) {
                                            //     std::cout << "  [Early Prune] Edge index " << current_edge_idx 
                                            //         << ": Remaining edges cannot be routed. Early prune count: " << result.pruned_early_count << std::endl;
                                            // }
                                            ++result.pruned_early_count;
                                            p = d.edges.back().built; // edge e fails, remove and try alt paths
                                            d.remove_edge();

                                            // keep trying next path 
                                            if (d.next_path(p, v, pcr)) d.add_edge(p, v, pcr);
                                            // until no alternatie paths
                                            else goto BACKUP;
                                        }
                                    }
                                }


                                if (++e == local_edges.end()) {
                                    if (constrained == 1)  {
                                        // =========================================================
                                        // PASS 2: CONSTRAINED COMPLETE DRAWING FOUND
                                        // =========================================================
                                        if (config_.verbose) {
                                            std::cout << ">>> [Pass 2] SUCCESS: Found Complete Drawing #" << result.full_solution_count 
                                                << " (from parent #" << solcount << ") <<<" << std::endl;
                                        }


                                        if (config_.export_files) {
                                            std::string base_path = config_.output_dir + "/solution_" + std::to_string(result.full_solution_count);

                                            std::ofstream of_gml(base_path + ".graphml");
                                            d.graphml_output(of_gml);

                                            std::ofstream of_json(base_path + ".json");
                                            nlohmann::ordered_json output_json = d.serialize_to_json();
                                            of_json << output_json.dump(4);
                                        }

                                        ++result.full_solution_count;
                                        goto BACKUP;
                                    }

                                    // map halfedges to faces
                                    std::vector<int> face(d.halfedges.size(), -1);
                                    std::size_t num_faces = 0;
                                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                                        if (face[i->label] != -1) continue;
                                        const HdsHalfedge* j = &*i;
                                        do {
                                            face[j->label] = static_cast<int>(num_faces);
                                            j = j->next;
                                        } while (j->label != i->label);
                                        ++num_faces;
                                    }
                                    if (num_faces == 0) goto BACKUP;

                                    // Flow network vertex indexing:
                                    //  - Faces: 0 ... num_faces - 1
                                    //  - Active cycle vertices: num_faces ... num_faces + config_.cycle_size - 1
                                    //  - Source node: num_faces + config_.cycle_size
                                    std::size_t source = num_faces + config_.cycle_size;
                                    DualNetwork net(source + 1);
                                    std::vector<std::vector<std::size_t>> dualg(num_faces);

                                    // source node to active cycle vertex
                                    for (std::size_t i = 0; i < config_.cycle_size; ++i)
                                        net.add_edge(source, num_faces + i, 1);

                                    // edges between faces and from (active cycle) vertices to incident faces
                                    std::size_t active_cycle_end = nm + config_.cycle_size;
                                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                                        // active cycle vertex -> incident face
                                        if (i->vertex->label >= nm && i->vertex->label < active_cycle_end)
                                            net.add_edge(num_faces + (i->vertex->label - nm), face[i->label], 1);

                                        // face -> adjacent face dual edge
                                        int remaining_capacity = static_cast<int>(klim) - static_cast<int>(i->edge->ncr);
                                        if (remaining_capacity > 0) {
                                            net.add_edge(face[i->label], face[i->twin->label], remaining_capacity);
                                            dualg[face[i->label]].push_back(face[i->twin->label]);
                                        }
                                    }

                                    // Potential final faces
                                    std::vector<int> pff;
                                    for (std::size_t f = 0; f < num_faces; ++f) {
                                        if (static_cast<std::size_t>(net.flow(static_cast<int>(source), static_cast<int>(f))) >= config_.cycle_size) {
                                            pff.push_back(static_cast<int>(f));
                                        }
                                    }
                                    if (pff.empty()) goto BACKUP;

                                    // Face classification:
                                    //   active = 3 -> possible final face (PFF)
                                    //   active = 2 -> active face
                                    //   active = 1 -> passive face
                                    //   active = 0 -> transit face
                                    //   active = -1 -> irrelevant face
                                    std::vector<int> active(num_faces, -1);
                                    for (int target_f : pff) active[target_f] = 3;

                                    // Active faces: flow >= 3 to a PFF
                                    for (std::size_t i = 0; i < num_faces; ++i) {
                                        if (active[i] != -1) continue;
                                        for (int target_f : pff) {
                                            if (net.flow(i, target_f) >= 3) {
                                                active[i] = 2; break;
                                            }
                                        }
                                    }

                                    // Compute map: face -> representative boundary halfedge
                                    std::vector<const HdsHalfedge*> fhedge(num_faces, nullptr);
                                    for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                                        if (fhedge[face[i->label]] == nullptr) fhedge[face[i->label]] = &*i;
                                    }

                                    // Passive faces: incident to active cycle vertex v and (0,1,2)-step dual path to an active face/PFF
                                    for (std::size_t i = 0; i < num_faces; ++i) {
                                        if (fhedge[i] == nullptr) throw std::runtime_error("no edge for face");
                                        if (active[i] != -1) continue;
                                        const HdsHalfedge* e_curr = fhedge[i];
                                        do {
                                            std::size_t v = e_curr->vertex->label;
                                            if (v >= nm && v < active_cycle_end) // active vertex
                                                for (const HdsHalfedge* f = e_curr->next->next; f != e_curr; f = f->next) {
                                                    if (f->edge->ncr < klim &&
                                                            f->vertex->label != v &&
                                                            f->twin->vertex->label != v) 
                                                    {
                                                        int fn = face[f->twin->label];
                                                        if (active[fn] >= 2) {active[i] = 1; break;}
                                                        for (std::size_t neighbor_f : dualg[fn])
                                                            if (active[neighbor_f] >= 2) {active[i] = 1; break;}
                                                    }
                                                }

                                            e_curr = e_curr->next;
                                        } while (active[i] == -1 && e_curr != fhedge[i]);
                                    }


                                    // Transit faces: adjacent to one active and one (active||passive) face
                                    for (std::size_t i = 0; i < num_faces; ++i) {
                                        if (active[i] > 0) continue;
                                        int an = 0, pn = 0; // active neighbor, passive neighbor count
                                        for (std::size_t neighbor_f : dualg[i])
                                            if (active[neighbor_f] >= 2) ++an;
                                            else if (active[neighbor_f] == 1) ++pn;
                                        if (an >= 2 || (an == 1 && pn >= 1))
                                            active[i] = 0;
                                    }

                                    // determine relevant vertices
                                    // ensure that the active cycle vertices offset, ..., offset+12-1, are mapped to 0...13 in relv
                                    std::vector<const HdsHalfedge*> relv(config_.cycle_size, nullptr); // relevant halfedges
                                    std::vector<const HdsHalfedge*> prelv; // possibly relevant halfedges

                                    for (auto i = d.vertices.begin(); i != d.vertices.end(); ++i) {
                                        auto e = i->halfedge;
                                        std::size_t nirf = 0; // #incident relevant faces
                                        auto a = e;
                                        do {
                                            if (active[face[e->label]] >= 0) { ++nirf; a = e; }
                                            e = e->next->twin;
                                        } while (e != i->halfedge);

                                        if (i->label >= nm && i->label < active_cycle_end) {
                                            // if no incident face is relevant, then 
                                            // this isnt a valid solution
                                            if (nirf == 0) goto BACKUP;
                                            relv[(i->label) - nm] = a;
                                        } else if (nirf >= 2)
                                            relv.push_back(a);
                                        else if (nirf == 1)
                                            prelv.push_back(a);
                                    }

                                    // ... and crossings
                                    for (auto i = d.crossings.begin(); i != d.crossings.end(); ++i) {
                                        auto e = i->halfedge;
                                        std::size_t nirf = 0;
                                        auto a = e;
                                        do {
                                            if (active[face[e->label]] >= 0) { ++nirf; a = e; }
                                            e = e->next->twin;
                                        } while (e != i->halfedge);

                                        if (nirf >= 2)
                                            relv.push_back(a);
                                        else if (nirf == 1)
                                            prelv.push_back(a);
                                    }

                                    // make sure triangles are not contracted to (non-simple) lenses
                                    std::vector<std::size_t> fsize(num_faces, 0); // size of relevant faces
                                    std::size_t total_nv = d.vertices.size() + d.crossings.size();
                                    std::vector<std::size_t> vact(total_nv, total_nv); // vertex mapping: old label -> new label
                                    std::size_t vic = 0;

                                    for (auto i = relv.begin(); i != relv.end(); ++i) {
                                        auto j = *i;
                                        do {
                                            ++fsize[face[j->label]];
                                            j = j->next->twin;
                                        } while (j != *i);
                                        vact[(*i)->vertex->label] = vic++;
                                    }

                                    for (auto i = prelv.begin(); i != prelv.end(); ++i)
                                        if (fsize[face[(*i)->label]] < 3) { // if lense, add prelv halfedge
                                            ++fsize[face[(*i)->label]];
                                            relv.push_back(*i);
                                            vact[(*i)->vertex->label] = vic++;
                                        }

                                    Drawing<klim> nd(relv.size() + 1); // extract relevant parts into pruned drawing
                                    std::deque<const HdsHalfedge*> todo(1, relv[0]);
                                    // vertex status: -1 == not considered yet, 0 == in queue and in nd, 
                                    // 1 == handled and incident edges in nd
                                    std::vector<int> done(d.vertices.size() + d.crossings.size(), -1);
                                    // We merge connected regions of irrelevant faces into "black" regions. Record them here. 
                                    std::vector<HdsHalfedge*> black;

                                    while (!todo.empty()) {
                                        auto i = todo.front();
                                        todo.pop_front();
                                        std::size_t u = i->vertex->label;

                                        if (active[face[i->label]] < 0)
                                            throw std::runtime_error("irrelevant face");

                                        HdsHalfedge* last = nd.vertices[vact[u]].halfedge;
                                        if (nd.edges.empty()) {
                                            // first edge
                                            auto p = prev_active<klim>(i, face, active, vact);
                                            std::size_t v = p.first->vertex->label;
                                            last = nd.add_first_edge(vact[u], vact[v], p.second)->twin;
                                            if (last->vertex->label != 0)
                                                throw std::runtime_error("wrong first edge in new drawing");
                                            todo.push_back(p.first);
                                            done[v] = 0;
                                        } else if (last == nullptr) {
                                            todo.push_back(i);
                                            continue;
                                        }

                                        // act[u] is now connected in nd and last points to act[u]; 
                                        // we build neighborhood of act[u] in nd
                                        done[u] = 1;
                                        std::size_t wnd = last->twin->vertex->label; // neighbor of u in nd
                                        std::size_t w = 0; // label of wnd in d
                                        auto wi = i;
                                        do {
                                            auto e_prev = prev_active<klim>(wi, face, active, vact).first;
                                            w = e_prev->vertex->label;
                                            if (vact[w] == wnd) break;
                                            wi = wi->next->twin;
                                        } while (wi != i);
                                        if (vact[w] != wnd) throw std::runtime_error("no w");

                                        // edges incident to u in d
                                        for (auto j = wi->next->twin; j != wi; j = j->next->twin) {
                                            if (active[face[j->label]] < 0 && active[face[j->twin->label]] < 0)
                                                continue;

                                            auto p = prev_active<klim>(j, face, active, vact);
                                            std::size_t v = p.first->vertex->label;
                                            if (last->next->vertex->label == vact[v]) {
                                                // edge uv already present
                                                last = last->next->twin;
                                                continue;
                                            }

                                            HdsPath path(2, last);
                                            if (!nd.find_target(path, vact[v]))
                                                throw std::runtime_error("cannot find v");

                                            auto ne = nd.add_edge(path, vact[v], p.second);
                                            if (active[face[j->label]] < 0)
                                                black.push_back(ne->twin);
                                            else if (active[face[j->twin->label]] < 0)
                                                black.push_back(ne);

                                            last = last->next->twin;
                                            if (done[v] < 0) {
                                                done[v] = 0;
                                                todo.push_back(p.first);
                                            }
                                        }
                                    }

                                    // process black faces
                                    std::vector<bool> blackdone(nd.halfedges.size(), false);
                                    std::vector<HdsHalfedge*> large_black_faces;

                                    for (auto i = black.begin(); i != black.end(); ++i) {
                                        if (blackdone[(*i)->label]) continue;
                                        auto j = *i;
                                        std::size_t bc = 0;
                                        do {
                                            blackdone[j->label] = true;
                                            ++bc;
                                            j = j->next;
                                        } while (j != *i);

                                        if (bc >= 4) large_black_faces.push_back(*i);
                                    }
                                    // allocate extra star vertices if more than one black face
                                    if (large_black_faces.size() > 1)
                                        nd.add_vertices(large_black_faces.size() - 1);

                                    // add an uncrossable star in each (large) black face
                                    std::size_t star_v = relv.size();
                                    for (auto start_edge : large_black_faces) {
                                        auto j = start_edge;
                                        // add spoke to first boundary vertex
                                        auto x = nd.add_edge(HdsPath({j, nullptr}), star_v, klim);
                                        for (;;) {
                                            j = j->next->twin->next; // jump over spoke
                                            if (j == start_edge) break;
                                            nd.add_edge(HdsPath({j, x}), star_v, klim);
                                        }
                                        ++star_v;
                                    }


                                    for (auto x = result.solutions.begin(); x != result.solutions.end(); ++x) {
                                        if (are_isomorphic(x->drawing, nd)) {
                                            if (config_.verbose) {
                                                std::cout << "--- Discard drawing #" << result.discarded_count
                                                    << ", isomorphic to solution #"
                                                    << (x - result.solutions.begin()) << std::endl;
                                            }

                                            // std::ofstream of;
                                            // std::ostringstream filename;
                                            // filename << "discard-" << discarded << ".graphml";
                                            // of.open(filename.str());
                                            // d.graphml_output(of);
                                            // of.close();

                                            ++result.discarded_count;
                                            goto BACKUP;
                                        }
                                    }

                                    // we have a new, valid solution -> record it
                                    if (!nd.is_valid()) throw std::runtime_error("nd is invalid");
                                    result.solutions.push_back({nd, current_depth+1});
                                    is_extensible = true; // Mark as extensible to enable Pass 2
                                    std::cout << "Drawing #" << result.solutions.size()-1 << ":\n"
                                        << d << std::endl;

                                    if (pff.size() > 1)
                                        std::cout << "!!! Final face is not unique for this drawing ---"
                                            << std::endl;

                                    if (config_.verbose) {
                                        std::cout << "[Pass 1] Added Intermediate Drawing #" << result.solutions.size() - 1 
                                            << " to Queue." << std::endl;
                                    }

                                    goto BACKUP; 
                                } // if (++e == local_edges.end())
                            } // for (auto e = edges.begin();;)

FINISH_PASS:
                            std::cout << "Finished Pass " << (constrained == 0 ? "1" : "2") 
                                << " for Drawing #" << solcount 
                                << " (Total DFS Iterations: " << total_local_iterations << ")" << std::endl;
                        } // for (int constrained = 0; constrained < 2; constrained++)
                        ++solcount;
                    } // while (solcount < solutions.size())

                    result.total_processed = solcount;
                    return result;
                }

            private:
                Config config_;
        };
} // namespace hds

#endif // NESTED_CYCLE_SEARCH_HPP


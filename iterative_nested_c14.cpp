#include "hds_kplanar.h"
#include "iso.h"
#include <iostream>
#include <vector>
#include <deque>
#include <sstream>
#include <fstream>
#include <cassert>
#include <nlohmann/json.hpp>

// BGL Flow
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

const std::size_t klim = 3; // 3-planar
const std::size_t CYCLE_SIZE = 14;
const std::size_t MAX_CYCLES = 2;


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

Drawing<klim> create_base_drawing() {
    // 14 vertices for C14_1 + 1 dummy star vertex = 15 vertices
    Drawing<klim> d(CYCLE_SIZE + 1);
    std::vector<HdsHalfedge*> cycle(CYCLE_SIZE, nullptr);

    // build uncrossable C14 cycle by setting ncr=klim=3
    cycle[1] = d.add_first_edge(0, 1, klim);
    for (std::size_t i = 2; i < CYCLE_SIZE; ++i) {
        cycle[i] = d.add_edge(HdsPath({cycle[i - 1], nullptr}), i, klim);
    }
    cycle[0] = d.add_edge(HdsPath({cycle[CYCLE_SIZE - 1], cycle[1]->twin}), 0, klim);

    // block one side of the cycle using an uncrossable star at dummy vertex 14
    std::size_t star_center = CYCLE_SIZE;
    auto e = d.add_edge(HdsPath({cycle[0], nullptr}), star_center, klim);
    for (std::size_t i = CYCLE_SIZE - 1; i > 0; --i) {
        e = d.add_edge(HdsPath({cycle[i], e}), star_center, klim);
    }

    return d;
}

// Returns a halfedge pointing to the desired source vertex, and the number 
// of crossings the edge should have in the reduced drawing.
//
// Used when a halfedge has an irrelevant endpoint:
//   vact[v] == vact.size()   iff   v is irrelevant
    template <typename F, typename A, typename V>
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

int main() {
    std::vector< Drawing<klim> > solutions;
    Drawing<klim> base = create_base_drawing();
    solutions.push_back(base);

    // #discarded solutions (isomorphic to a known solution)
    std::size_t discarded = 0;
    // we process the solutions in order, this is the current one we process
    std::size_t solcount = 0;
    // #full solutions
    std::size_t fullsolcount = 0;

    while (solcount < solutions.size()) {
        std::cout << "\n--- Processing Partial Drawing #" << solcount 
            << " (Queue Size: " << solutions.size() << ") ---" << std::endl;

        Drawing<klim> d_parent = solutions[solcount];
        std::size_t offset = d_parent.vertices.size();

        std::size_t current_cycle = (offset - 1) / CYCLE_SIZE;
        if (MAX_CYCLES > 0 && current_cycle >= MAX_CYCLES) {
            std::cout << "--> Reached MAX_CYCLES limit (" << MAX_CYCLES 
                      << "). Skipping further cycle extension for Drawing #" << solcount << std::endl;
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
                std::cout << "--> Skipping Constrained Pass 2 (Drawing #" << solcount 
                    << " is not extensible)" << std::endl;
                continue;
            }
            std::cout << "\n>>> Starting Pass " << (constrained == 0 ? "1 (Unconstrained)" : "2 (Constrained)")
                << " for Drawing #" << solcount << " <<<" << std::endl;

            Drawing<klim> d = d_parent;
            std::size_t nm14 = d.vertices.size();
            d.add_vertices(CYCLE_SIZE);
            std::vector<Edge> local_edges = {
                {0,nm14+11}, // one matching
                // outer cycle
                {nm14+11,nm14+12,3},
                {nm14+12,nm14+13,3},
                {nm14+13,nm14   ,3},
                {nm14   ,nm14+1 ,3},
                {nm14+1 ,nm14+2 ,3},
                {nm14+2 ,nm14+3 ,3},
                {nm14+3 ,nm14+4 ,3},
                {nm14+4 ,nm14+5 ,3},
                {nm14+5 ,nm14+6 ,3},
                {nm14+6 ,nm14+7 ,3},
                {nm14+7 ,nm14+8 ,3},
                {nm14+8 ,nm14+9 ,3},
                {nm14+9 ,nm14+10,3},
                {nm14+10,nm14+11,3},
                // rest of matching
                // EVEN: (i, nm14 + ( (i+11) % 14)) == (i, nm14 + ( (i-3) % 14))
                // ODD: (i, nm14 + ( (i+3) % 14)) == (i, nm14 + ( (i-11) % 14))
                {2,nm14+13},{4,nm14+1},{6,nm14+3},{8,nm14+5},{10,nm14+7},{12,nm14+9},
                {1,nm14+4},{3,nm14+6},{5,nm14+8},{7,nm14+10},{9,nm14+11},{11,nm14+0},{13,nm14+2}
            };

            auto e = local_edges.begin();
            uint64_t total_local_iterations = 0;

            for (;;) {
                ++total_local_iterations;
                std::size_t edge_idx = std::distance(local_edges.begin(), e);
                std::size_t u = (*e)[0];
                std::size_t v = (*e)[1];

                // In Pass 2 (constrained == 1), newly added C14 cycle edges (indices >= CYCLE_SIZE)
                // are uncrossable (pcr = klim). In Pass 1, edges are unconstrained (pcr = 0).
                std::size_t pcr = (constrained == 1 && edge_idx >= CYCLE_SIZE) ? klim : 0;
                HdsPath p = d.first_path(u, v, pcr);

                if (p.empty()) {
BACKUP:
                    do {
                        if (e == local_edges.begin()) {
                            goto FINISH_PASS;
                        }
                        --e;

                        edge_idx = std::distance(local_edges.begin(), e);
                        pcr = (constrained == 1 && edge_idx >= CYCLE_SIZE) ? klim : 0;

                        u = (*e)[0];
                        v = (*e)[1];
                        p = d.edges.back().built;
                        d.remove_edge();

                        if (d.next_path(p, v, pcr)) {
                            break;
                        }
                    } while (true);
                }
                d.add_edge(p, v, pcr);

                if (++e == local_edges.end()) {
                    if (constrained == 0) {
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
                        //  - Active cycle vertices: num_faces ... num_faces + CYCLE_SIZE - 1
                        //  - Source node: num_faces + CYCLE_SIZE
                        std::size_t source = num_faces + CYCLE_SIZE;
                        DualNetwork net(source + 1);
                        std::vector<std::vector<std::size_t>> dualg(num_faces);

                        // source node to active cycle vertex
                        for (std::size_t i = 0; i < CYCLE_SIZE; ++i)
                            net.add_edge(source, num_faces + i, 1);

                        // edges between faces and from (active cycle) vertices to incident faces
                        std::size_t active_cycle_end = offset + CYCLE_SIZE;
                        for (auto i = d.halfedges.begin(); i != d.halfedges.end(); ++i) {
                            // active cycle vertex -> incident face
                            if (i->vertex->label >= offset && i->vertex->label < active_cycle_end)
                                net.add_edge(num_faces + (i->vertex->label - offset), face[i->label], 1);

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
                            if (net.flow(static_cast<int>(source), static_cast<int>(f)) >= 14) {
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
                                if (v >= offset && v < active_cycle_end) // active vertex
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
                        // ensure that the active cycle vertices offset, ..., offset+14-1, are mapped to 0...13 in relv
                        std::vector<const HdsHalfedge*> relv(CYCLE_SIZE, nullptr); // relevant halfedges
                        std::vector<const HdsHalfedge*> prelv; // possibly relevant halfedges

                        for (auto i = d.vertices.begin(); i != d.vertices.end(); ++i) {
                            auto e = i->halfedge;
                            std::size_t nirf = 0; // #incident relevant faces
                            auto a = e;
                            do {
                                if (active[face[e->label]] >= 0) { ++nirf; a = e; }
                                e = e->next->twin;
                            } while (e != i->halfedge);

                            if (i->label >= offset && i->label < offset + CYCLE_SIZE) {
                                // if no incident face is relevant, then 
                                // this isnt a valid solution
                                if (nirf == 0) goto BACKUP;
                                relv[(i->label) - offset] = a;
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
                                auto p = prev_active(i, face, active, vact);
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
                                auto e_prev = prev_active(wi, face, active, vact).first;
                                w = e_prev->vertex->label;
                                if (vact[w] == wnd) break;
                                wi = wi->next->twin;
                            } while (wi != i);
                            if (vact[w] != wnd) throw std::runtime_error("no w");

                            // edges incident to u in d
                            for (auto j = wi->next->twin; j != wi; j = j->next->twin) {
                                if (active[face[j->label]] < 0 && active[face[j->twin->label]] < 0)
                                    continue;

                                auto p = prev_active(j, face, active, vact);
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

                        // Add an uncrossable star in each (large) black face
                        std::vector<bool> blackdone(nd.halfedges.size(), false);
                        std::size_t nbf = 0; // #black faces on >=4 vertices
                        for (auto i = black.begin(); i != black.end(); ++i) {
                            if (blackdone[(*i)->label]) continue;
                            auto j = *i;
                            std::size_t bc = 0;
                            // mark all halfedges of this blackface as done
                            do {
                                blackdone[j->label] = true;
                                ++bc;
                                j = j->next;
                            } while (j != *i);

                            if (bc >= 4) {
                                // for large black regions, add an uncrossable star center vertex
                                if (++nbf >= 2)
                                    throw std::runtime_error("too many large black regions");

                                j = *i;
                                auto x = nd.add_edge(HdsPath({j, nullptr}), nd.vertices.size() - 1, klim);
                                for (;;) {
                                    j = j->next->twin->next;
                                    if (j == *i) break;
                                    nd.add_edge(HdsPath({j, x}), nd.vertices.size() - 1, klim);
                                }
                            }
                        }

                        for (auto x = solutions.begin(); x != solutions.end(); ++x) {
                            if (are_isomorphic(*x, nd)) {
                                std::cout << "--- Discard drawing #" << discarded
                                    << ", isomorphic to solution #"
                                    << (x - solutions.begin()) << std::endl;

                                std::ofstream of;
                                std::ostringstream filename;
                                filename << "discard-" << discarded << ".graphml";
                                of.open(filename.str());
                                d.graphml_output(of);
                                of.close();

                                ++discarded;
                                goto BACKUP;
                            }
                        }

                        // we have a new, valid solution -> record it
                        if (!nd.is_valid()) throw std::runtime_error("nd is invalid");
                        solutions.push_back(nd);
                        is_extensible = true; // Mark as extensible to enable Pass 2
                        std::cout << "Drawing #" << solutions.size()-1 << ":\n"
                            << d << std::endl;

                        if (pff.size() > 1)
                            std::cout << "!!! Final face is not unique for this drawing ---"
                                << std::endl;

                        std::cout << "[Pass 1] Added Intermediate Drawing #" << solutions.size() - 1 
                                  << " to Queue." << std::endl;

                        goto BACKUP; 
                    } 
                    else {
                        // =========================================================
                        // PASS 2: CONSTRAINED COMPLETE DRAWING FOUND
                        // =========================================================
                        std::cout << ">>> [Pass 2] SUCCESS: Found Complete Drawing #" << fullsolcount 
                            << " (from parent #" << solcount << ") <<<" << std::endl;

                        // Export complete terminal solution
                        std::string gml_path = "../quasiDrawings/nested_c14/solution_" + std::to_string(fullsolcount) + ".graphml";
                        std::ofstream of_gml(gml_path);
                        d.graphml_output(of_gml);
                        of_gml.close();

                        std::string json_path = "../quasiDrawings/nested_c14/solution_" + std::to_string(fullsolcount) + ".json";
                        std::ofstream of_json(json_path);
                        nlohmann::ordered_json output_json = d.serialize_to_json();
                        of_json << output_json.dump(4);
                        of_json.close();

                        ++fullsolcount;
                        goto BACKUP;

                    }
                }
            }

FINISH_PASS:
            std::cout << "Finished Pass " << (constrained == 0 ? "1" : "2") 
                      << " for Drawing #" << solcount 
                      << " (Total DFS Iterations: " << total_local_iterations << ")" << std::endl;
        }
        ++solcount;
    }

    std::cout << "\n==================================================" << std::endl;
    std::cout << "Search Complete!" << std::endl;
    std::cout << "Total Intermediate Subdrawings in Queue: " << solutions.size() << std::endl;
    std::cout << "Total Full Terminal Solutions Found: " << fullsolcount << std::endl;
    std::cout << "==================================================" << std::endl;
    return 0;
}

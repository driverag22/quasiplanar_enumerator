#ifndef HDS_HPP
#define HDS_HPP
#include <vector>
#include <map>
#include <list>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <string>
#include <cassert> 
#include <stdint.h>
#include <nlohmann/json.hpp>

// classic Halfedge Data Structure

struct HdsVertex;
struct HdsEdge;
struct HdsHalfedge 
{
	HdsVertex* vertex;	// target vertex
	HdsHalfedge* prev;	// previous halfedge along the face
	HdsHalfedge* next;	// next halfedge along the face
	HdsHalfedge* twin;	// twin halfedge
	HdsEdge* edge;		// the underlying graph edge
	std::size_t label;	// unique label for each halfedge
}
;
struct HdsVertex 
{
	HdsVertex(HdsHalfedge* h, std::size_t l) : halfedge(h), label(l) {}
	HdsHalfedge* halfedge;	// some halfedge pointing to *this
	std::size_t label;		// unique label for each vertex
}
;
// Creation paths:
//
// To be able to remove edges, we store how they were inserted. This
// is represented by a so-called creation path, which consists of
// between 2 and kplane + 2 halfedges: The first halfedge points to the
// source vertex and the insertion position in the rotation at that
// vertex; symmetrically, the last halfedge determines the target
// vertex and the position of the inserted edge in the rotation there;
// the zero, one, two or three halfedges in between are crossed along the
// way from source to target.
//
// To insert an edge to an isolated vertex we use the convention that
// the last halfedge of the creation path is 0 (null pointer).
//
// Specifically, for a creation path p of length l, we have:
//
// (1) If l == 2, either p[1] == 0 or:
//    (a) p[0] and p[1] are on the same face of the drawing.
//    (b) p[1] is not incident to the vertex that p[0] points to.
//
// (2) If l > 2, for every intermediate crossing step i (where 1 <= i <= l - 2):
//    (a) The underlying edge of p[i] currently has strictly fewer than kplane crossings.
//    (b) The "entry door" (which is p[0] if i==1, or p[i-1]->twin otherwise) and the 
//        "exit door" p[i] are on the same face of the drawing.
//    (c) The underlying edges of all crossed halfedges (p[1] through p[i]) are distinct  
//        (i.e., a path cannot cross the same edge twice).
//
// (3) For the final step (reaching the target), either p[l-1] == 0, or:
//    (a) p[l-2]->twin and p[l-1] are on the same face of the drawing.
//    (b) p[l-1] is not incident to the vertex that p[0] points to.
//
// (1) l >= 2;
// (2) either l == 2 and p[1] == 0 or
//   (a) p[0] and p[1] are on the same face of the drawing and p[1] is
//     not incident to the vertex that p[0] points to;
// [only for kplane >= 1:]
// (3) If l >= 3, then either l == 3 and p[2] == 0 or
//   (a) the underlying edge of p[1] has strictly fewer than kplane
//       crossings in the drawing;
//   (b) p[1]->twin and p[2] are on the same face of the drawing and
//       p[2] is not incident to the vertex that p[0] points to;
// [only for kplane >= 2:]
// (4) If l >= 4, then either l == 4 and p[3] == 0 or
//   (a) the underlying edge of p[1] and p[2] has strictly fewer than kplane
//       crossings in the drawing;
//   (b) p[2]->twin and p[3] are on the same face of the drawing and
//       p[3] is not incident to the vertex that p[0] points to;
//   (c) the underlying edges of p[1] and p[2] are distinct.
// [only for kplane >= 3:]
// (5) If l >= 5, then either l == 5 and p[4] == 0 or
//   (a) the underlying edge of p[1], p[2] and p[3] has strictly fewer than kplane
//       crossings in the drawing;
//   (b) p[3]->twin and p[4] are on the same face of the drawing and
//       p[4] is not incident to the vertex that p[0] points to;
//   (c) the underlying edges of p[1], p[2] and p[3] are distinct.
//
// When building a possible creation path, we deal with objects that
// are not creation paths yet, but we want to express that, at least,
// they are valid up to this point. We call such an object a creation
// path stub; the conditions for it are the same as for a creation
// path, with one exception: the last halfedge may be incident to the
// vertex that p[0] points to.

typedef std::vector<HdsHalfedge*> HdsPath;

struct HdsEdge 
{
	HdsEdge(std::size_t u_, std::size_t v_,
		  const HdsPath& p, std::size_t c, std::size_t l_)
	    : u(u_), v(v_), ncr(c), built(p), label(l_) 
	{}
	std::size_t u, v;	// indices of the two endpoints
	std::size_t ncr;	// #crossings
	HdsPath built;		// how the edge is built into the drawing
	std::size_t label;	// unique label for each edge
}
;

// Most of the actual functionality is built into the Drawing
// class. It stores the containers of all objects and defines the
// interface to interact with them.
//
// Vertices are static and created upon initialization, although not
// yet connected by edges. Crossings may appear and disappear when
// edges are inserted and removed, respectively; hence, they are
// stored separately, in a list. Halfedges and edges are also stored
// in lists. The advantage of lists over, say, vectors/arrays is that
// pointers remain valid under insertions and removals of other
// elements. Insertion into a vector may lead to reallocation, which
// invalidates *all* pointers into it.
//
// The parameter kplane should be \geq 0.

template < int kplane >
struct Drawing {
	// n == #vertices
	Drawing(std::size_t n) 
	{
		vertices.reserve(n);
		for (std::size_t i = 0; i < n; ++i)
		    vertices.emplace_back((HdsHalfedge*)(0), i);
    cross_mat.assign(64, 0); // init to 64 edges
	}
	
	// copy constructor; many pointers, we have to recompute all of them
	Drawing(const Drawing& d)
	    : vertices(d.vertices),
	      crossings(d.crossings),
	      halfedges(d.halfedges),
	      edges(d.edges),
        cross_mat(d.cross_mat)
	{
		// build indices
		std::vector<HdsHalfedge*> ind(halfedges.size());
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
		    ind[i->label] = &*i;
		std::vector<HdsEdge*> ind_e(edges.size());
		for (auto i = edges.begin(); i != edges.end(); ++i) 
		    ind_e[i->label] = &*i;
		std::vector<HdsVertex*> ind_c(crossings.size());
		for (auto i = crossings.begin(); i != crossings.end(); ++i) 
		    ind_c[i->label - vertices.size()] = &*i;
		    
		// fix pointers
		for (std::size_t i = 0; i < vertices.size(); ++i) 
		{
			if (vertices[i].label != i)
				throw std::runtime_error("wrong vertex label");
			if (vertices[i].halfedge != 0)
				vertices[i].halfedge = ind[vertices[i].halfedge->label];
		}
		std::size_t count = 0;
		for (auto i = crossings.begin(); i != crossings.end(); ++i,++count) 
		{
			if (i->label < vertices.size())
				throw std::runtime_error("crossing label too small");
			if (i->label - vertices.size() != count)
				throw std::runtime_error("wrong crossing label");
			if (ind_c[i->label - vertices.size()] != &*i)
				throw std::runtime_error("pointer error");
			if (i->halfedge == 0)
				throw std::runtime_error("crossing without halfedge");
			i->halfedge = ind[i->halfedge->label];
		}
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
		{
			if (i->vertex->label < vertices.size())
				i->vertex = &(vertices[i->vertex->label]); 
			else 
				i->vertex = ind_c[i->vertex->label - vertices.size()];
			i->prev = ind[i->prev->label];
			i->next = ind[i->next->label];
			i->twin = ind[i->twin->label];
			i->edge = ind_e[i->edge->label];
		}
		for (auto i = edges.begin(); i != edges.end(); ++i) 
		    for (auto x = i->built.begin(); x != i->built.end(); ++x)
				if (*x != 0) 
					*x = ind[(*x)->label];
	}

    // Construct drawing from json recipe
    Drawing(const nlohmann::json& root, const std::size_t target_n = 0) {
        if (!root.contains("kplane") || !root.contains("drawing_recipe") || !root.contains("num_vertices")) {
            throw std::runtime_error("Deserialization Constructor Error: Missing essential JSON keys.");
        }

        int json_kplane = root["kplane"];
        if (json_kplane > kplane) {
            throw std::runtime_error("Deserialization Constructor Error: Compile-time kplane (" + 
                    std::to_string(kplane) + ") is strictly larger than JSON kplane (" + 
                    std::to_string(json_kplane) + ").");
        }

        std::size_t num_v = root["num_vertices"];
        if (target_n > num_v) num_v = target_n;
        for (std::size_t i = 0; i < num_v; ++i) {
            vertices.push_back(HdsVertex(nullptr, i));
        }

        cross_mat.assign(64, 0);
        std::map<std::size_t, HdsEdge*> label_to_edge_map; // map JSON labels to edges

        const auto& recipe = root["drawing_recipe"];
        for (const auto& step : recipe) {
            std::size_t u = step["u"];
            std::size_t v = step["v"];
            std::size_t edge_label = step["edge_label"];

            // first edge
            if (edge_label == 0) {
                add_first_edge(u, v);
                label_to_edge_map[edge_label] = &(edges.back());
                continue;
            }

            std::size_t start_after = step["start_after_edge"];

            // find the anchor edge
            auto it = label_to_edge_map.find(start_after);
            if (it == label_to_edge_map.end()) {
                throw std::runtime_error("Deserialization error: Anchor edge label " + std::to_string(start_after) + " not mapped.");
            }
            HdsEdge* anchor_edge = it->second;

            // find the incoming halfedge segment belonging to anchor_edge incident to vertex u
            HdsHalfedge* p0 = nullptr;
            for (auto& he : halfedges) {
                if (he.edge == anchor_edge && he.vertex->label == u) {
                    p0 = &he;
                    break;
                }
            }

            if (!p0) {throw std::runtime_error("Deserialization error: could not find start_after_edge ID " + 
                    std::to_string(start_after) + " at vertex " + std::to_string(u));}

            HdsPath p;
            p.push_back(p0);

            // reconstruct intermediate crossings using correct face-walking logic
            const auto& crossed = step["crossed_edges"];
            HdsHalfedge* face_runner = p0; // start from p0

            for (const auto& crossed_label_json : crossed) {
                std::size_t crossed_label = crossed_label_json;

                auto cross_it = label_to_edge_map.find(crossed_label);
                if (cross_it == label_to_edge_map.end()) {
                    throw std::runtime_error("Deserialization error: Crossed edge label " + std::to_string(crossed_label) + " not mapped.");
                }
                HdsEdge* target_cross_edge = cross_it->second;

                HdsHalfedge* found_crossing = nullptr;
                HdsHalfedge* start_face = face_runner;

                // walk around face until we find edge to cross
                do {
                    face_runner = face_runner->next;
                    if (face_runner->edge == target_cross_edge) {
                        found_crossing = face_runner;
                        break;
                    }
                } while (face_runner != start_face);

                if (!found_crossing) {throw std::runtime_error("Deserialization error: could not find crossed edge ID " + 
                        std::to_string(crossed_label) + " along the current face boundary.");}

                p.push_back(found_crossing);
                // after crossing, enter the adjacent face via the twin halfedge
                face_runner = found_crossing->twin;
            }

            // reconstruct target pointer p[l-1]
            if (vertices[v].halfedge == nullptr) {
                p.push_back(nullptr); // target vertex is isolated
            } else {
                HdsHalfedge* found_target = nullptr;
                HdsHalfedge* start_face = face_runner;

                // walk around face until we find target vertex
                do {
                    face_runner = face_runner->next;
                    if (face_runner->vertex->label == v) {
                        found_target = face_runner;
                        break;
                    }
                } while (face_runner != start_face);

                if (!found_target) {throw std::runtime_error("Deserialization error: could not find target vertex " + 
                        std::to_string(v) + " in the current face, from: " + std::to_string(u));}
                p.push_back(found_target);
            }

            add_edge(p, v, 0);
            label_to_edge_map[edge_label] = &(edges.back());
        }
    }

	
	// no copy-assignment allowed
	Drawing& operator=(const Drawing&) = delete;
	
	// add x new vertices to the drawing
	// PRE: The drawing must not contain any crossings.
	void add_vertices(std::size_t x) 
	{
		if (!crossings.empty())
		    throw std::runtime_error("no crossings allowed when adding vertices");
		
		// As vertices are stored in a vector, inserting into it
		// potentially invalidates all vertex pointers. Thus, we have to
		// recompute them.
		
		// store vertex labels for halfedges
		std::vector<std::size_t> hind(halfedges.size());
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i)
			hind[i->label] = i->vertex->label;
		    
		// expand the drawing
		std::size_t n = vertices.size();
		vertices.reserve(n+x);
		for (std::size_t i = 0; i < x; ++i)
		    vertices.emplace_back((HdsHalfedge*)(0), n+i);
		    
		// build index
		std::vector<HdsVertex*> vind(vertices.size(), 0);
		for (auto i = vertices.begin(); i != vertices.end(); ++i)
		    vind[i->label] = &*i;
		    
		// consolidate pointers
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
		    i->vertex = vind[hind[i->label]];
	}
	
	// create and return a new pair of twin halfedges that correspond to
	// the edge e
	HdsHalfedge* new_twin(HdsEdge* e) 
	{
		halfedges.emplace_back();
		halfedges.back().label = halfedges.size()-1;
		HdsHalfedge* nhe = &(halfedges.back());
		halfedges.emplace_back();
		halfedges.back().label = halfedges.size()-1;
		HdsHalfedge* nhet = &(halfedges.back());
		nhe->edge = nhet->edge = e;
		nhe->twin = nhet;
		return nhet->twin = nhe;
	}
	
	// split h into a path of two halfedges by inserting a crossing on
	// it, so that h points to this new crossing. The underlying edge of
	// the graph remains the same (for both halfedges along the path).
	void split_edge(HdsHalfedge* h) 
	{
		crossings.emplace_back(h, vertices.size() + crossings.size());
		HdsVertex* nv = &(crossings.back());
		HdsHalfedge* nhe = new_twin(h->edge);
		HdsHalfedge* nhet = nhe->twin;
		nhe->vertex = h->vertex;
		nhe->prev = h;
		nhet->vertex = nv;
		nhet->next = h->twin;
		if (h->vertex->halfedge == h) 
			h->vertex->halfedge = nhe;
		h->vertex = nv;
		nv->halfedge = h;
		if (h->next == h->twin) 
		{
			nhe->next = nhet;
			nhet->prev = nhe;
		} 
		else 
		{
			nhe->next = h->next;
			nhe->next->prev = nhe;
			nhet->prev = h->twin->prev;
			nhet->prev->next = nhet;
		}
		nhe->prev->next = nhe;
		nhet->next->prev = nhet;
		if (++(h->edge->ncr) > kplane)
			throw std::runtime_error("too many crossings per edge");
	}
	
	// opposite of split_edge above: contract the path of two halfedges
	// formed by h and h->next.
	// PRE: h points to the last crossing created.
	void join_vertex(HdsHalfedge* h) 
	{
		if (h->edge->ncr == 0) 
			throw std::runtime_error("internal error");
		if (h->vertex != &(crossings.back()))
		    throw std::runtime_error("only the last crossing can be removed");
		HdsHalfedge* nhe = h->next;
		HdsHalfedge* nhet = nhe->twin;
		if (nhe != &(halfedges.back()) && nhet != &(halfedges.back()))
		    throw std::runtime_error("only the last halfedge can be removed");
		h->vertex = nhe->vertex;
		if (nhe->vertex->halfedge == nhe) 
			nhe->vertex->halfedge = h;
		if (nhe->next == nhet) 
		{
			h->next = h->twin;
			h->twin->prev = h;
		} 
		else 
		{
			h->next = nhe->next;
			h->next->prev = h;
			h->twin->prev = nhet->prev;
			h->twin->prev->next = h->twin;
		}
		--(h->edge->ncr);
		crossings.pop_back();
		halfedges.pop_back();
		halfedges.pop_back();
	}
	
private:
	
	// helper function to insert edges: add and return a new halfedge
	// stub that corresponds to the edge e, and add it in the rotation
	// of s->vertex right behind s.
	//
	// It is a "stub" because it is only properly connected to the rest
	// of the HDS at its source s, the other endpoint must be addressed
	// by the calling function.
	HdsHalfedge* create_starting_stub_at(HdsHalfedge* s, HdsEdge* e) 
	{
		// create new edge twin
		HdsHalfedge* nhe = new_twin(e);
		HdsHalfedge* nhet = nhe->twin;
		// insert at source
		nhet->vertex = s->vertex;
		nhet->next = s->next;
		nhet->next->prev = nhet;
		nhe->prev = s;
		return nhe->prev->next = nhe;
	}
	
	// helper function to remove edges and the inverse of
	// create_starting_stub_at above: detach the edge (stub) right
	// behind s in the rotation of s->vertex
	void detach_stub(HdsHalfedge* s) 
	{
		s->next = s->next->twin->next;
		s->next->prev = s;
	}
	
public:
	// add a first edge to the drawing, with pcr prescribed (artifical)
	// crossings; return the halfedge that points to v
	HdsHalfedge* add_first_edge(std::size_t u,
				      std::size_t v,
				      std::size_t pcr = 0) 
	{
		if (!edges.empty())
		      throw std::runtime_error("inserting another first edge");
		// create new edge twin (the creation path for this edge will be
		// built manually once the two halfedges have been created)
		edges.emplace_back(u, v, HdsPath(), pcr, 0);
		HdsHalfedge* nhe = new_twin(&(edges.back()));
		HdsHalfedge* nhet = nhe->twin;
		nhe->prev = nhe->next = nhet;
		nhet->prev = nhet->next = nhe;
		nhe->vertex = &(vertices[u]);
		nhet->vertex = &(vertices[v]);
		nhe->vertex->halfedge = nhe;
		nhet->vertex->halfedge = nhet;
		edges.back().built.emplace_back(nhe);
		edges.back().built.emplace_back(nhet);
		return nhet;
	}
	
	// add an edge described by p (which determines source and potential
	// crossings) with target v to the drawing, with pcr prescribed
	// (artifical) crossings; return the halfedge that points to v
	HdsHalfedge* add_edge(const HdsPath& p,
				std::size_t v,
				std::size_t pcr = 0) 
	{
		if (p.size() > 2 + pcr + kplane)
			throw std::runtime_error("inserting an edge with too many crossings");
		if (p.size() < 2) 
			throw std::runtime_error("edge not properly described");
		//std::cerr << "add: " << p[0]->vertex->label << ">" << v << " (cr = "
		// << pcr << ")" << std::endl;
		edges.emplace_back(p[0]->vertex->label, v, p, pcr, edges.size());

    // update cross_mat
    std::size_t new_label = edges.back().label;
    for (std::size_t i = 1; i < p.size(); ++i) {
      if (i < p.size() - 1 && p[i] != 0) {
        std::size_t crossed_label = p[i]->edge->label;
        cross_mat[new_label] |= (1ULL << crossed_label);
        cross_mat[crossed_label] |= (1ULL << new_label);
      }
    }

		HdsHalfedge* nhe = create_starting_stub_at(p[0], &(edges.back()));
		for (std::size_t i = 1; i < p.size(); ++i) 
		{
			HdsHalfedge* nhet = nhe->twin;
			if (i < p.size()-1) 
			{
				split_edge(p[i]);
				++(edges.back().ncr);
			} 
			else if (p[i] == 0) 
			{
				// this is a new vertex
				if (vertices[v].halfedge != 0)
					throw std::runtime_error("vertex is already present");
				nhe->next = nhet;
				nhet->prev = nhe;
				nhe->vertex = &(vertices[v]);
				nhe->vertex->halfedge = nhe;
				edges.back().built.back() = nhe;
				break;
			}
			nhe->next = p[i]->next;
			nhe->next->prev = nhe;
			nhet->prev = p[i];
			nhet->prev->next = nhet;
			nhe->vertex = p[i]->vertex;
			if (i < p.size()-1) 
				nhe = create_starting_stub_at(nhe->next->twin, &(edges.back()));
		}
		return nhe;
	}
	
	// remove the last edge (i.e., the edge that was inserted last)
	void remove_edge() 
	{
    // cross_mat cleanup
    std::size_t removed_label = edges.back().label;
    uint64_t crossed = cross_mat[removed_label];
    for (int b = 0; b < 64; ++b) { // hard-coded 64 edge limit
        if ((crossed >> b) & 1) {
            cross_mat[b] &= ~(1ULL << removed_label);
        }
    }
    cross_mat[removed_label] = 0; // wipe the removed edge's row

		HdsPath p = edges.back().built;
		if (p.size() < 2) 
			throw std::runtime_error("no creation path");
		if (p[0]->next->edge != &(edges.back()))
		    throw std::runtime_error("creation path mismatch");
		HdsVertex* v = p.back()->vertex;
		if (v->halfedge->next == v->halfedge->twin) 
		    // this is the only edge incident to v
			v->halfedge = 0; 
		else 
		    detach_stub(p.back());
		for (std::size_t i = p.size()-2;; --i) 
		{
			if (i == 0) 
			{
				detach_stub(p[i]);
				halfedges.pop_back();
				halfedges.pop_back();
				break;
			}
			detach_stub(p[i]->next->twin->next->twin);
			halfedges.pop_back();
			halfedges.pop_back();
			detach_stub(p[i]);
			join_vertex(p[i]);
		}
		edges.pop_back();
	}
	
	// Try to extend the creation path stub p by walking p.back() around
	// its face using ->next, so as to reach v; return true iff
	// successful (p.back() is updated along the way, regardless).
	// PRE: p.size() >= 2 and p is a creation path stub.
	bool find_target(HdsPath& p, std::size_t v) 
	{
		assert(p.size() >= 2);
		if (vertices[v].halfedge == 0) 
		{
			// new vertex, we can just insert it into this face
			p.back() = 0;
			return true;
		}
		// walk p.back() around its face
		std::size_t i = p.size()-1;
		HdsHalfedge* end = (i == 1 ? p[0] : p[i-1]->twin);
    std::vector<uint8_t> jumped(i, 0);

		for (;;) 
		{
			p[i] = p[i]->next;
      std::fill(jumped.begin(), jumped.end(), 0);

      // avoid self-cross
			bool changed = true;

			while (changed && i > 1) {
				changed = false;

				for (std::size_t x = i - 1; x > 0; --x) {
          if (jumped[x] == 1) continue;

					HdsHalfedge* v_start = (x == 1) ? p[0] : p[x-1]->twin;
					HdsHalfedge* v_end = p[x];

          if (p[i] == v_end) { 
            p[i] = v_start; 
            jumped[x] = 1;
            changed = true; 
            break;
          } else if (p[i] == v_start) { 
            p[i] = v_end; 
            jumped[x] = 1;
            changed = true; 
            break;
          }
				}
			}

			if (p[i] == end) 
				return false;
			if (p[i]->vertex->label == v) 
				return true;
		}
		return false;	// avoid compiler warning
	}
	
private:
	// General recursive pathfinding function
	bool find_crossing(HdsPath& p, std::size_t v, std::size_t ci, std::size_t max_crossings) 
	{
		assert(p.size() >= ci + 2);
		std::size_t u = p[0]->vertex->label;
		HdsHalfedge* stop_condition = (ci == 1) ? p[0] : p[ci-1]->twin;
		
    std::vector<uint8_t> jumped(ci, 0);
		for (;;) 
		{
			p[ci] = p[ci]->next;
      std::fill(jumped.begin(), jumped.end(), 0);
			
			// "virtual" face boundaries (self-cross)
			bool changed = true;

			while (changed && ci > 1) {
				changed = false;

				for (std::size_t x = ci - 1; x > 0; --x) {
          if (jumped[x] == 1) continue;

					HdsHalfedge* v_start = (x == 1) ? p[0] : p[x-1]->twin;
					HdsHalfedge* v_end = p[x];

          if (p[ci] == v_end) { 
            p[ci] = v_start; 
            jumped[x] = 1; // Mark as jumped
            changed = true; 
            break; // Restart the while loop safely
          }
          else if (p[ci] == v_start) { 
            p[ci] = v_end; 
            jumped[x] = 1; // Mark as jumped
            changed = true; 
            break; // Restart the while loop safely
          }
				}
			}
			
			if (p[ci] == stop_condition) break;
			

      // avoid crossing edges with overlapping nodes
			if (p[ci]->edge->u == u || p[ci]->edge->u == v ||
				p[ci]->edge->v == u || p[ci]->edge->v == v ||
        p[ci]->edge->ncr >= max_crossings) {
        continue;
			}

			bool invalid = false;
      // dont cross edge twice
			for (std::size_t x = 1; x < ci; ++x) {
				if (p[ci]->edge->label == p[x]->edge->label) {
					invalid = true; break;
				}
			}
			if (invalid) continue;

			// quasiplanarity check
			bool qp_violation = false;
			std::size_t e_new = p[ci]->edge->label; // edge we cross
			for (std::size_t x = 1; x < ci; ++x) { // prev edgges
				std::size_t e_old = p[x]->edge->label;
				if ((cross_mat[e_new] & (1ULL << e_old)) != 0) {
          qp_violation = true; 
          break;
        }
			}
			if (qp_violation) continue;

			// move to next step
			p[ci+1] = p[ci]->twin;
			
			if (p.size() == ci + 2) { // last crossing of this edge
				if (find_target(p, v)) return true;
			} else if (p.size() > ci + 2) { // recurse
				if (find_crossing(p, v, ci + 1, max_crossings)) return true;
			}
		}
		return false;
	}

    // Helper backtracking function for the matching verification
    bool find_matching_dfs(const std::vector<const HdsEdge*>& uncrossed, 
            std::vector<bool>& vertex_used, 
            std::size_t edge_idx, 
            std::size_t current_size, 
            std::size_t target_size) const {
        if (current_size == target_size) return true; // base case
        if (edge_idx >= uncrossed.size()) return false; // out of edges

        // if the remaining edges aren't enough to reach the target, return
        if (current_size + (uncrossed.size() - edge_idx) < target_size) return false;

        const HdsEdge* e = uncrossed[edge_idx];

        // option 1: include edge in the matching (if both endpoints are free)
        if (!vertex_used[e->u] && !vertex_used[e->v]) {
            vertex_used[e->u] = true;
            vertex_used[e->v] = true;

            if (find_matching_dfs(uncrossed, vertex_used, edge_idx + 1, current_size + 1, target_size)) return true;

            // Backtrack
            vertex_used[e->u] = false;
            vertex_used[e->v] = false;
        }

        // option 2: skip edge
        return find_matching_dfs(uncrossed, vertex_used, edge_idx + 1, current_size, target_size);
    }

public:
	// try to build a first creation path for an edge from u to v with c
	// prescribed (artificial) crossings; return the path, or an empty
	// path if no such creation path exists.
	// PRE: u has degree at least one in the drawing
	HdsPath first_path(std::size_t u, std::size_t v, std::size_t pcr = 0) 
	{
		HdsHalfedge* end = vertices[u].halfedge;
		if (end == 0)
			throw std::runtime_error("don't start an edge from an isolated vertex");
		HdsPath p;
		p.emplace_back(end);
		p.emplace_back(end);
		// try to get to v with no crossings
		do 
		{
			if (find_target(p, v)) return p;
			p[1] = p[0] = p[0]->twin->prev;
		}while (p[0] != end);
    // continue exploring with crossings, up to the (heuristic) depth limit
		for(std::size_t cr = 1; cr + pcr <= kplane; cr++) {
			p.emplace_back();
			do {
				if (find_crossing(p, v, 1, kplane)) return p;
				p[1] = p[0] = p[0]->twin->prev;
			} while (p[0] != end);
		}
		return HdsPath();
	}
	
	// try to find another creation path that extends the creation path
	// stub p to reach v, with pcr prescribed (artificial) crossings on
	// this edge. Return true iff another path was found.
	// PRE: (1) p.size() >= 2 and (2) p[0],...,p[p.size()-2] is a
	//      creation path stub.
	//
	// We enumerate creation paths the following way:
	// (1) shorter paths (i.e., with fewer crossings) before longer
	//     paths;
	// (2) explore the possibilities along each face boundary using
	//     ->next, until reaching the starting point again.
	bool next_path(HdsPath& p, std::size_t v, std::size_t pcr = 0) 
	{
		assert(p.size() >= 2);
		HdsHalfedge* end = p[0]->vertex->halfedge;
		bool v_is_isolated = (vertices[v].halfedge == 0);

    // try to find another option on the same face
		if (!v_is_isolated && find_target(p, v)) return true;
		
		// iteratively backtrack through the crossings
		for (int ci = p.size() - 2; ci >= 1; --ci) {
			if (find_crossing(p, v, ci, kplane)) return true;
		}
		
		// all crossings failed for current start edge, try a new start edge
		if (p.size() > 2) {
			for (;;) {
				p[1] = p[0] = p[0]->twin->prev;
				if (p[0] == end) break;
				if (find_crossing(p, v, 1, kplane)) return true;
			}
		} else {
			for (;;) {
				p[1] = p[0] = p[0]->twin->prev;
				if (p[0] == end) break;
				if (find_target(p, v)) return true;
			}
		}
		
		// if all start edges failed at current depth, increment crossing limit
		for (std::size_t cr = p.size() - 1; cr + pcr <= kplane; cr++) {
			p.emplace_back();
			do {
				if (find_crossing(p, v, 1, kplane)) {
          return true;
        }
				p[1] = p[0] = p[0]->twin->prev;
			} while (p[0] != end);
		}
		
		return false;
	}
	
	
	// debug output of all halfedges
	void dump() const {
		std::cerr << "{ #vertices = " << vertices.size()
			      << ", #edges = " << edges.size()
			      << ", #halfedges = " << halfedges.size()
			      << ", #crossings = " << crossings.size() << "\n";
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
		      std::cerr << i->label << "\n"
				<< "\tprev = " << i->prev->label << "\n"
				<< "\tnext = " << i->next->label << "\n"
				<< "\ttwin = " << i->twin->label << "\n"
				<< "\tvert = " << i->vertex->label << "\n"
				<< "\tedge = " << i->edge->u << "-" << i->edge->v << "\n";
		std::size_t c = 0;
		for (auto i = edges.begin(); i != edges.end(); ++i, ++c)
		      std::cerr << "e" << c << " " << i->u << " - " << i->v << "\n";
		std::cerr << "}\n";
	}
	// output a human readable recipe to draw the drawing, edge by edge
	void construction() {
		for (auto i = edges.begin(); i != edges.end(); ++i) {
			HdsPath& p = i->built;
			std::cerr << "draw an edge leaving " << i->u
					<< " cw after the edge to "
					<< (p[0]->edge->u == i->u ? p[0]->edge->v : p[0]->edge->u);
			for (std::size_t j = 1; j < p.size()-1; ++j) {
				// determine crossing orientation
				HdsHalfedge* src = p[j];
				while (src->vertex->label >= vertices.size())
					  src = src->next->twin->next;
				HdsHalfedge* tgt = p[j]->twin;
				while (tgt->vertex->label >= vertices.size())
					  tgt = tgt->next->twin->next;
				std::cerr << "\n\tacross the edge " << src->vertex->label
						  << "-" << tgt->vertex->label;
			}
			std::cerr << "\n\tto " << i->v << "\n";
      std::cerr << "\n --------------------------- \n";
		}
	}

	// graphml output of the drawing, to be viewed with yEd, for instance
	std::ostream& graphml_output(std::ostream& o) const 
	{
		std::size_t total_nv = vertices.size() + crossings.size();
		std::vector<std::size_t> vcol(total_nv, total_nv);
		for (std::size_t i = 0; i < vertices.size(); ++i) vcol[i] = 0;
		return graphml_output(o, vcol);
	}
	
	std::ostream& graphml_output(std::ostream& o,
				       const std::vector<std::size_t>& vcol) const 
	{
		return graphml_output_footer(graphml_output_header(o, vcol));
	}
	
	std::ostream& graphml_output_header(std::ostream& o,
					      const std::vector<std::size_t>& vcol) const 
	{
		o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		      << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
		      << "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
		      << "    xmlns:y=\"http://www.yworks.com/xml/graphml\"\n"
		      << "    xsi:schemaLocation"
		      << "=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
		      << "  <key for=\"node\" id=\"d5\" yfiles.type=\"nodegraphics\"/>\n"
		      << "  <key for=\"edge\" id=\"d9\" yfiles.type=\"edgegraphics\"/>\n"
		      << "  <graph id=\"G\" edgedefault=\"undirected\">\n";
		for (std::size_t i = 0; i < vertices.size(); ++i)
		      o << "    <node id=\"v" << i << "\">\n"
			<< "      <data key=\"d5\">\n"
		        << "        <y:ShapeNode>\n"
			<< "          <y:NodeLabel>" << i << "</y:NodeLabel>\n"
			<< "          <y:Geometry height=\"30.0\" width=\"30.0\"/>\n"
			<< "          <y:Fill color=\"#"
			<< (vcol[i] < vcol.size() ? "99CCFF" : "C0C0C0")
			<< "\"/>\n"
			<< "          <y:Shape type=\"ellipse\"/>\n"
		        << "        </y:ShapeNode>\n"
			<< "      </data>\n"
			<< "    </node>\n";
		for (std::size_t i = 0; i < crossings.size(); ++i)
		      o << "    <node id=\"c" << i << "\">\n"
			<< "      <data key=\"d5\">\n"
		        << "        <y:ShapeNode>\n"
			<< "          <y:NodeLabel>" << i << "</y:NodeLabel>\n"
			<< "          <y:Geometry height=\"15.0\" width=\"15.0\"/>\n"
			<< "          <y:Fill color=\"#"
			<< (vcol[i+vertices.size()] < vcol.size() ? "FF6600" : "C0C0C0")
			<< "\"/>\n"
			<< "          <y:Shape type=\"rectangle\"/>\n"
		        << "        </y:ShapeNode>\n"
			<< "      </data>\n"
			<< "    </node>\n";
		std::vector<std::size_t> edgeparts(edges.size(), 0);
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) 
		      if (i->label > i->twin->label) ++edgeparts[i->edge->label];
		std::size_t ec = 0;
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
			if (i->label > i->twin->label) continue;
			std::size_t u = i->vertex->label;
			std::size_t v = i->twin->vertex->label;
			std::size_t cr = edgeparts[i->edge->label] - 1;
			o << "    <edge id=\"e" << ec++ << "\" source=\"";
			if (u < vertices.size())
				o << "v" << u; else o << "c" << u - vertices.size();
			o << "\" target=\"";
			if (v < vertices.size())
				o << "v" << v; else o << "c" << v - vertices.size();
			o << "\"";
			std::string color; // >= 3 crossings: blue
			std::string width = "1.5";
      if (cr == 0) {
        color = "00CC00"; // Green
      } else if (cr == 1) {
        color = "FFCC00"; // Yellow
      } else if (cr == 2) {
        color = "FFA500"; // Orange
      } else if (cr == 3) {
        color = "FF0000"; // Red
      } else if (cr == 4) {
        color = "9933CC"; // Purple
      } else {
        color = "0000FF"; // Blue
      }

			o << ">\n"
				<< "      <data key=\"d9\">\n"
				<< "        <y:PolyLineEdge>\n"
				<< "          <y:LineStyle color=\"#" << color << "\" "
				<<                "type=\"line\" width=\"" << width << "\"/>\n"
				<< "          <y:BendStyle smoothed=\"true\"/>\n"
				<< "        </y:PolyLineEdge>\n"
				<< "      </data>\n"
				<< "    </edge>\n";
		}
		return o;
	}
	std::ostream& graphml_output_footer(std::ostream& o) const {
		return o << "  </graph>\n</graphml>\n";
	}
	bool is_valid() const {
		// build index
		std::vector<const HdsHalfedge*> ind(halfedges.size());
		std::size_t count = 0;
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
			ind[i->label] = &*i;
			if (i->label != count++) return false;
		}
		std::vector<const HdsEdge*> ind_e(edges.size());
		count = 0;
		for (auto i = edges.begin(); i != edges.end(); ++i) {
			ind_e[i->label] = &*i;
			if (i->label != count++) return false;
		}
		std::vector<const HdsVertex*> ind_c(crossings.size());
		count = 0;
		for (auto i = crossings.begin(); i != crossings.end(); ++i) {
			ind_c[i->label - vertices.size()] = &*i;
			if (i->label - vertices.size() != count++) return false;
		}
		count = 0;
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
			if (i->vertex->label < vertices.size() &&
				  i->vertex != &(vertices[i->vertex->label]))
				return false;
			if (i->vertex->label >= vertices.size() &&
				  i->vertex != ind_c[i->vertex->label - vertices.size()])
				return false;
			if (i->prev != ind[i->prev->label]) return false;
			if (i->next != ind[i->next->label]) return false;
			if (i->twin != ind[i->twin->label]) return false;
			if (i->edge != ind_e[i->edge->label]) return false;
			if (i->next->prev != &*i) return false;
			if (i->prev->next != &*i) return false;
			if (i->twin->twin != &*i) return false;
		}
		count = 0;
		for (auto i = vertices.begin(); i != vertices.end(); ++i,++count) {
			if (i->label != count) return false;
			//if (i->halfedge == 0) return false;
			if (i->halfedge != 0) {
				if (i->halfedge->vertex->label != i->label) return false;
				if (i->halfedge != ind[i->halfedge->label]) return false;
			}
		}
		for (auto i = crossings.begin(); i != crossings.end(); ++i) {
			if (i->halfedge == 0 || i->halfedge->vertex->label != i->label)
				return false;
			if (i->halfedge != ind[i->halfedge->label]) return false;
		}
		// check that halfedges correctly implement their underlying edges
		for (auto i = halfedges.begin(); i != halfedges.end(); ++i) {
			std::size_t u = i->edge->u;
			if (i->vertex->label != u) continue;
			std::size_t v = i->edge->v;
			std::size_t iel = i->edge->label;
			auto e = i->twin;
			if (e->edge->label != iel) return false;
			if (e->vertex->label == v) continue;
			if (e->vertex->label < vertices.size()) return false;
			e = e->next->twin->next;
			if (e->edge->label != iel) return false;
			if (e->vertex->label == v) continue;
			if (e->vertex->label < vertices.size()) return false;
			e = e->next->twin->next;
			if (e->edge->label != iel) return false;
			if (e->vertex->label != v) return false;
		}
		return true;
	}

  // Verifies if the current drawing is 3-quasiplanar, without using cross_mat.
  // Returns true if NO three edges pairwise cross (triangle-free crossing graph).
  bool verify_quasiplanarity() const {
    std::size_t num_edges = edges.size();

    // crossing_graph[i][j] == true if edge i and edge j cross
    std::vector<std::vector<bool>> crossing_graph(num_edges, std::vector<bool>(num_edges, false));

    for (const auto& crossing : crossings) {
        HdsHalfedge* e0 = crossing.halfedge;
        if (!e0 || !e0->edge) continue;
        
        HdsHalfedge* e1 = e0->twin->prev; // the intersecting path choice
        if (!e1 || !e1->edge) continue;

        std::size_t label1 = e0->edge->label;
        std::size_t label2 = e1->edge->label;

        if (label1 < num_edges && label2 < num_edges) {
            crossing_graph[label1][label2] = true;
            crossing_graph[label2][label1] = true;
        }
    }

    // check for triangle
    for (std::size_t i = 0; i < num_edges; ++i) {
      for (std::size_t j = i + 1; j < num_edges; ++j) {
        if (crossing_graph[i][j]) {
          for (std::size_t k = j + 1; k < num_edges; ++k) {
            if (crossing_graph[i][k] && crossing_graph[j][k]) {
                auto it_i = std::next(edges.begin(), i);
                auto it_j = std::next(edges.begin(), j);
                auto it_k = std::next(edges.begin(), k);
                std::cout << "Quasiplanarity violation by " << i << " " << j <<  " " << k << ":\n";
                std::cout << it_i->u << " to " << it_i->v << "\n";
                std::cout << it_j->u << " to " << it_j->v << "\n";
                std::cout << it_k->u << " to " << it_k->v << std::endl;
                return false; 
            }
          }
        }
      }
    }

    return true;
  }

  // Returns true if all vertices are incident on an uncrossed edge
  bool verify_vertex_non_crossing_edge() const {
    std::size_t n = vertices.size();
    std::vector<bool> has_uncrossed(n, false);

    for (const auto& edge : edges) {
        if (edge.ncr == 0) {
            has_uncrossed[edge.u] = 1;
            has_uncrossed[edge.v] = 1;
        }
    }
    std::size_t covered_vertices = std::count(has_uncrossed.begin(), has_uncrossed.end(), true);

    return covered_vertices == n;
  }

  // Returns true if the current drawing contains a matching of non-crossing edges.
  bool verify_non_crossing_matching() const {
    std::size_t n = vertices.size();
    std::size_t target_size = n/2;

    std::vector<const HdsEdge*> uncrossed_edges;
    // get uncrossed edges
    for (const auto& edge : edges) {
        if (edge.ncr == 0) uncrossed_edges.push_back(&edge);
    }

    if (uncrossed_edges.size() < target_size) return false;

    std::vector<bool> vertex_used(n,false);

    // backtracking to find matching
    return find_matching_dfs(uncrossed_edges, vertex_used, 0, 0, target_size);
  }

  // Extract abstract graph as a list of canonical, sorted edges
  nlohmann::json extract_abstract_graph() const {
    nlohmann::json j_graph = nlohmann::json::array();
    for (const auto& edge : edges) {
      std::size_t low = std::min(edge.u, edge.v);
      std::size_t high = std::max(edge.u, edge.v);
      j_graph.push_back(nlohmann::json::array_t{low, high});
    }
    return j_graph;
  }

  // Extract the order and topological context of edge insertions
  nlohmann::json extract_drawing_recipe() const {
    nlohmann::json j_recipe = nlohmann::json::array();

    for (const auto& edge : edges) {
      nlohmann::json step;
      step["edge_label"] = edge.label;
      step["u"] = edge.u;
      step["v"] = edge.v;

      if (edge.label != 0) {
        // edge.built[0] points to the halfedge we insert after clockwise
        if (edge.built.empty() || edge.built[0] == nullptr || edge.built[0]->edge == nullptr) {
          throw std::runtime_error("Malformed creation path at edge " + std::to_string(edge.label));
        }
        step["start_after_edge"] = edge.built[0]->edge->label;

        // get crossed halfedges (from 1 to edge.built.size()-1)
        auto crossed_arr = nlohmann::json::array();
        for (std::size_t j = 1; j < edge.built.size() - 1; ++j) {
          if (edge.built[j] != nullptr && edge.built[j]->edge != nullptr) {
            crossed_arr.push_back(edge.built[j]->edge->label);
          }
        }
        step["crossed_edges"] = crossed_arr;
      }
      j_recipe.push_back(step);
    }
    return j_recipe;
  }

  nlohmann::ordered_json serialize_to_json() const {
    nlohmann::ordered_json root;

    root["kplane"] = kplane;
    root["num_vertices"] = vertices.size();
    root["abstract_graph"] = extract_abstract_graph();
    root["drawing_recipe"] = extract_drawing_recipe();
    return root;
  }

  std::vector<HdsVertex> vertices;
  std::list<HdsVertex> crossings;
  std::list<HdsHalfedge> halfedges;
  std::list<HdsEdge> edges;

  std::vector<uint64_t> cross_mat; // edge crossing matrix, assumes at most 64 edges
}
;

// text output of rotation system
template < int kplane >
std::ostream& operator<<(std::ostream& o, const Drawing<kplane>& d) 
{
	o << "{ #vertices = " << d.vertices.size()
	    << ", #edges = " << d.edges.size()
	    << ", #halfedges = " << d.halfedges.size()
	    << ", #crossings = " << d.crossings.size() << "\n";
	for (auto i = d.vertices.begin(); i != d.vertices.end(); ++i) 
	{
		o << "\t" << i->label << " : ";
		HdsHalfedge* e = i->halfedge;
		if (e != 0)
		do 
		{
			std::size_t target = (e->edge->u == i->label ? e->edge->v : e->edge->u);
			o << target;
			HdsHalfedge* f = e->twin;
			if (f->vertex->label != target) {
				o << "(";
				for (; f->vertex->label != target; f = f->next->twin->next) {
					// determine crossing orientation
					HdsHalfedge* g = f->next;
					std::size_t c1 = g->edge->u;
					std::size_t c2 = g->edge->v;
					for (;;) {
						if (g->vertex->label == c1) break;
						if (g->vertex->label == c2) {
							std::swap(c1, c2);
							break;
						}
						g = g->next->twin->next;
					}
					o << "x" << c1 << "-" << c2;
				}
				o << ")";
			}
			o << " ";
			e = e->twin->prev;
		}
		while (e != i->halfedge);
		o << "\n";
	}
	return o << "}\n";
}
#endif

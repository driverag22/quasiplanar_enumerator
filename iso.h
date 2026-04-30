#ifndef ISO_HPP
#define ISO_HPP

#include "hds_3_planar.h"
#include <vector>
#include <stdexcept>
#include <queue>
#include <tuple>

// return true iff the two given drawings d1 and d2 are strongly
// isomorphic, i.e., both the graphs.h and their planarizations (where
// each crossing is replaced by a degree four vertex) are
// isomorphic. If the parameter mirror is set, then all rotations in
// one of the drawings are considered reversed.
// PRE: The drawing is connected.
template < int kplane >
bool are_isomorphic_oriented(const Drawing<kplane>& d1, const Drawing<kplane>& d2, bool mirror) 
{
  	typedef std::vector<std::size_t> SV;
  	typedef const HdsHalfedge* HEP;
  	// to get out of current attempt as soon as a problem is detected
  	struct Contradiction {}; 

  	// first check the raw numbers, which have to agree
  	std::size_t c = d1.crossings.size();
  	if (c != d2.crossings.size()) 
		return false;
  	std::size_t np = d1.vertices.size();
  	if (np != d2.vertices.size()) 
	  	return false;
  	std::size_t m = d1.edges.size();
  	if (m != d2.edges.size())
		return false;
  	std::size_t n = np + c;

  	// choosing the image for vertex 0 and one of its incident halfedges
  	// fixes the whole map; just try all options
  	for (std::size_t j = 0; j < np; ++j) 
	{
    	HEP h1 = d1.vertices[0].halfedge;
    	// try all halfedges pointing to j as an image for h0
    	HEP h2 = d2.vertices[j].halfedge;
    	HEP s = h2;
    	do 
		{
      		// initialize map
      		SV phi(n, n); // map d1 -> d2
      		SV psi(n, n); // == phi^{-1}
      		phi[0] = j;
      		psi[j] = 0;
      
      		// now do a simultaneous bfs from vertex 0 and j; for each
      		// vertex, queue the halfedge that reaches it first and its
      		// image, to be processed
      		std::queue<std::tuple<const HdsHalfedge*,const HdsHalfedge*> > queue;
      		queue.emplace(h1, s);
      		try 
			{
				while (!queue.empty()) 
				{
					const HdsHalfedge* s1 = std::get<0>(queue.front());
					const HdsHalfedge* s2 = std::get<1>(queue.front());
					queue.pop();
					const HdsHalfedge* t1 = s1;
					const HdsHalfedge* t2 = s2;
	  				do 
					{
	    				// crossings along matched edges must be the same
	    				if (t1->edge->ncr != t2->edge->ncr) 
						{ 
							throw Contradiction();
						}
	    				// look at the other endpoints
	    				std::size_t i1 = t1->twin->vertex->label;
	    				std::size_t i2 = t2->twin->vertex->label;
	    				// vertices to vertices and crossings to crossings
	    				if ((i1 < np && i2 >= np) || (i1 >= np && i2 < np)) 
						{ 
							throw Contradiction(); 
						}
	    				// if they don't have an image so far, now they do
	    				if (phi[i1] == n) 
						{
	      					if (psi[i2] < n) 
							{
								throw Contradiction(); 
							}
	      					phi[i1] = i2;
	      					psi[i2] = i1;
	      					queue.emplace(t1->twin, t2->twin);
	    				} 
						else if (phi[i1] != i2 || psi[i2] != i1)
						{
							throw Contradiction();
						}
	    				t1 = (mirror ? t1->twin->prev : t1->next->twin);
	    				t2 = t2->next->twin;
	  				} while (t1 != s1 && t2 != s2);
	  				if (t1 != s1 || t2 != s2) 
					{ 
						throw Contradiction(); 
					}
				}
				// found an isomorphism, make a sanity check
				for (std::size_t z = 0; z < n; ++z)
	  				if (phi[z] == n || psi[phi[z]] != z)
	    				throw std::runtime_error("internal error in isomorphy test");
				return true;
      		} 
			catch (Contradiction) 
			{}
      		// try next edge incident to vertex j
      		s = s->next->twin;
    	} while (s != h2);
  	}
  	return false;
}

template < int kplane >
inline bool are_isomorphic(const Drawing<kplane>& d1, const Drawing<kplane>& d2) 
{
	return are_isomorphic_oriented(d1, d2, true) || are_isomorphic_oriented(d1, d2, false);
}

#endif // ISO_HPP

#include "hds_3_planar.h"
#include "iso.h"
#include <sstream>
#include <fstream>

typedef std::vector<std::size_t> Edge;
typedef std::vector<Edge> Edges;

// Example graph K_9
// K_9 is not 3-planar
const std::size_t n = 5;
const Edges edges =
{
    {0,1},{0,2},{0,3},{0,4},//{0,5},{0,6},{0,7},{0,8},
    {1,2},{1,3},{1,4},//{1,5},{1,6},{1,7},{1,8},
    {2,3},{2,4},//{2,5},{2,6},{2,7},{2,8},
    {3,4},//{3,5},{3,6},{3,7},{3,8},
    //{4,5},{4,6},{4,7},{4,8},
    //{5,6},{5,7},{5,8},
    //{6,7},{6,8},
    //{7,8},
};
const std::size_t klim = 2;

int main()
{
    std::size_t minimal_cr = 0x3f3f3f3f;
    std::vector< Drawing<klim> > solutions;
    Drawing<klim> d(n);
    d.add_first_edge(edges[0][0], edges[0][1]);
    for (auto e = edges.begin() + 1;;)
    {
        std::size_t u = (*e)[0];
        std::size_t v = (*e)[1];
        HdsPath p = d.first_path(u, v);
        if (p.empty())
        {
            BACKUP:
            // no way to add uv -> do previous edges differently
            do
            {
                if (--e == edges.begin())
                    goto END;
                u = (*e)[0];
                assert(u == d.edges.back().u);
                v = (*e)[1];
                assert(v == d.edges.back().v);
                p = d.edges.back().built;
                d.remove_edge();
            }while (!d.next_path(p, v));
        }
        d.add_edge(p, v);
        if (++e == edges.end())
        {
            bool is_new = true;
            if (is_new)
            {

                solutions.push_back(d);
                minimal_cr = std::min(minimal_cr,d.crossings.size());
            }
            // goto BACKUP;
            goto END;
        }
    }
    END:
    std::cout << "Found " << solutions.size() << " drawings in total." << std::endl;
    if(solutions.size() == 0)
    {
    	return 0;
	}
    
    std::size_t cnt = 0;
    std::vector< Drawing<klim> > solutions_mincr_uni;
    std::vector< Drawing<klim>> solutions_mincr;
    // Assume the graph has no more than 100 different crossing minimal drawings
    std::vector<std::size_t> d_cnt(100,1);
    for (auto it = solutions.begin();it!=solutions.end();it++)
    {
		// Output all drawings with minimal crossings
        if(it->crossings.size() == minimal_cr)
        {
            cnt++;
            solutions_mincr.push_back((*it));
            bool is_unique = true;
            std::size_t d_ind = 0;
            for(auto iit = solutions_mincr_uni.begin();iit!=solutions_mincr_uni.end();iit++)
            {
                d_ind++;
                if(are_isomorphic((*it),(*iit)))
                {
                    is_unique = false;
                    d_cnt[d_ind]++;
                    break;
                }
            }
            if(is_unique)
            {
                solutions_mincr_uni.push_back((*it));
                std::ofstream of;
                std::ostringstream filename;
                filename << "drawings/Drawing_" << klim << "_" << n << "_" << solutions_mincr_uni.size() << ".graphml";
                of.open(filename.str());
                (*it).graphml_output(of);
                of.close();
            }
        }
    }
    std::cout << "Found " << cnt << " drawings with Minimal Crossing in total." << std::endl;
    std::cout << "Minimal Crossing Number is "<<minimal_cr<<std::endl;
    std::cout << "Found " << solutions_mincr_uni.size() << " unique drawings with Minimal Crossing in total." << std::endl;
    for (std::size_t i = 1; i <= solutions_mincr_uni.size(); i++)
    {
        std::cout << "Drawing-" << i << " has " << d_cnt[i] << " isomorphic drawings" << std::endl;
    }
    return 0;
}

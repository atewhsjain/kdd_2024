#ifndef ESCAPE_CYCLE_HELPER_H_
#define ESCAPE_CYCLE_HELPER_H_

#include <algorithm>
#include <chrono>
#include <queue>
#include <time.h>

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Utils.h"
#include "JointSort.h"

using namespace Escape;
using namespace std;
using namespace std::chrono;

#define NNBORS 1000
/* find neighbors nbr of ver such that nbr > v and nbr is not blocked. 
Note: if r is a neighbor that will also be included */
vector<VertexIdx>* find_neighbors(CDAG &DG, VertexIdx ver, VertexIdx v, vector<VertexIdx> &blocked)
{
    // cout << "In find_neighbors" << endl;
    vector<VertexIdx> *neighbors = new vector<VertexIdx>();
    neighbors->reserve(NNBORS);
    for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
    {
        VertexIdx nbr = DG.outlist.nbors[i];
        if ((blocked[nbr] == 0) && (nbr>v))
        {
            neighbors->push_back(nbr);
        }
    }
    for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
    {
        VertexIdx nbr = DG.inlist.nbors[i];
        if ((blocked[nbr] == 0) && (nbr>v))
        {
            neighbors->push_back(nbr);
        }
    }
    // cout << "returning from find_neighbors" << endl;
    return neighbors;
}
//returns non blocked out nbrs. This can include r.
vector<VertexIdx>* blockNeighbors(CDAG &DG, VertexIdx ver, VertexIdx v, vector<int> &blocked)
{
    vector<VertexIdx> *neighbors = new vector<VertexIdx>();
    neighbors->reserve(NNBORS);
    for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
    {
        VertexIdx nbr = DG.outlist.nbors[i];
        // if (nbr > v)
        // {
            if ((blocked[nbr] == 0) && nbr > v) neighbors->push_back(nbr);
            blocked[nbr]++;
        // }
    }
    for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
    {
        VertexIdx nbr = DG.inlist.nbors[i];
        // if (nbr > v)
        // {
            if ((blocked[nbr] == 0) && nbr > v) neighbors->push_back(nbr);
            blocked[nbr]++;
        // }
    }
    return neighbors;
}

void unblockNeighbors(CDAG &DG, VertexIdx ver, VertexIdx v, vector<int> &blocked)
{
    for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
    {
        VertexIdx nbr = DG.outlist.nbors[i];
        // if (nbr > v) 
            blocked[nbr]--;
    }
    for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
    {
        VertexIdx nbr = DG.inlist.nbors[i];
        // if (nbr > v) 
            blocked[nbr]--;
    }
}

bool isConnectedBFS(CDAG &DG, VertexIdx src, VertexIdx dest, VertexIdx v, vector<int> &blocked)
{
    queue<VertexIdx> q; 
    q.push(src);

    vector<VertexIdx> toUnblock;
    toUnblock.reserve(NNBORS);
    bool rFound = false;
    while(!q.empty())
    {
        VertexIdx ver = q.front();
        q.pop();

        for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
        {
            VertexIdx nbr = DG.outlist.nbors[i];
            if (nbr == dest)
            {
                rFound = true;
                break;
            }
            if (blocked[nbr] == 0)
            {
                if (nbr < v) cout << "nbr < v. This shouldn't happen" << endl;
                q.push(nbr);
                blocked[nbr] = 1;
                toUnblock.push_back(nbr);
            }
        }
        for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
        {
            VertexIdx nbr = DG.inlist.nbors[DG.inlist.offsets[ver]];
            if (nbr == dest)
            {
                rFound = true;
                break;
            }
            if (blocked[nbr] == 0)
            {
                if (nbr < v) cout << "nbr < v. This shouldn't happen" << endl;
                q.push(nbr);
                blocked[nbr] = 1;
                toUnblock.push_back(nbr);
            }
        }
    }
    for (int i=0; i<toUnblock.size(); i++)
    {
        blocked[toUnblock[i]] = 0;
    }
    if (rFound == true)
        return true;
    else return false;
    // return false;
}

int get_num_shortest_paths(CDAG &DG, VertexIdx v, VertexIdx l, vector<int> &B, int shortest_path_length)
{
    static VertexIdx n = DG.outlist.nVertices;
    static vector<int> count_paths(n);
    static vector<int> A(n);
    A = B;

    count_paths.assign(n,0);

    static vector<int> levels(n);
    levels.assign(n,numeric_limits<int>::infinity()-100);

    queue<pair<VertexIdx, int>> q; 
    // q.clear();
    q.push(make_pair(v,1));
    count_paths[v] = 1;
    levels[v] = 1;

    for (int i=DG.outlist.offsets[v]; i<(int)DG.outlist.offsets[v+1]; i++)
    {
        if (A[DG.outlist.nbors[i]] >= 1)
        {
            q.push(make_pair(DG.outlist.nbors[i],2));
            levels[DG.outlist.nbors[i]] = 2;
            count_paths[DG.outlist.nbors[i]] = 1;
            A[DG.outlist.nbors[i]] = 0;
            // B[DG.outlist.nbors[i]] = 1;
        }
    }

    while(!q.empty())
    {
        pair<VertexIdx,int> p = q.front();
        int level = p.second;
        if (level > shortest_path_length) break; //TODO:changed
        // if (level > max_cycle_size) return;

        VertexIdx ver = p.first;
        q.pop();
        // B[ver] = 1;
        levels[ver] = level;

        for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
        {
            if (A[DG.outlist.nbors[i]] >= 1)
            {
                q.push(make_pair(DG.outlist.nbors[i],level+1));
                A[DG.outlist.nbors[i]] = 0;
                count_paths[DG.outlist.nbors[i]] += count_paths[ver];
                // B[DG.outlist.nbors[i]] = 1;
            }
            // else if ((levels[DG.outlist.nbors[i]] != numeric_limits<int>::infinity()-100) && ())
            else if (levels[DG.outlist.nbors[i]] == level+1)
            {
                count_paths[DG.outlist.nbors[i]] += count_paths[ver];
            }

        }

        for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
        {
            if (A[DG.inlist.nbors[i]] >= 1)
            {
                q.push(make_pair(DG.inlist.nbors[i],level+1));
                A[DG.inlist.nbors[i]] = 0;
                count_paths[DG.inlist.nbors[i]] += count_paths[ver];

                // B[DG.outlist.nbors[i]] = 1;
            }
            else if (levels[DG.inlist.nbors[i]] == level+1)
            {
                count_paths[DG.inlist.nbors[i]] += count_paths[ver];
            }
        }
    }
    return count_paths[l];

}

// returns the sphere of v with old blocking
// void find_l_sphere(CDAG &DG, VertexIdx v, vector<int> &B, int max_cycle_size, vector<int> &levels)
// {
//     VertexIdx n = DG.outlist.nVertices;
//     static vector<int> A(n);
//     A = B;

//     queue<pair<VertexIdx, int>> q; 
//     q.push(make_pair(v,1));
//     A[v] = 0;
//     levels[v] = 1;
//     // VertexIdx n = DG.outlist.nVertices;

//     // vector<int> B = new vector<int>(n,0);

//     for (int i=DG.outlist.offsets[v]; i<(int)DG.outlist.offsets[v+1]; i++)
//     {
//         if (A[DG.outlist.nbors[i]] >= 1)
//         {
//             q.push(make_pair(DG.outlist.nbors[i],2));
//             levels[DG.outlist.nbors[i]] = 2;
//             A[DG.outlist.nbors[i]] = 0;
//             // B[DG.outlist.nbors[i]] = 1;
//         }
//     }

//     while(!q.empty())
//     {
//         pair<VertexIdx,int> p = q.front();
//         int level = p.second;
//         if (level > (max_cycle_size /2)+1) return; //TODO:changed
//         // if (level > max_cycle_size) return;

//         VertexIdx ver = p.first;
//         q.pop();
//         // B[ver] = 1;
//         levels[ver] = level;

//         for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
//         {
//             if (A[DG.outlist.nbors[i]] >= 1)
//             {
//                 q.push(make_pair(DG.outlist.nbors[i],level+1));
//                 A[DG.outlist.nbors[i]] = 0;
//                 // B[DG.outlist.nbors[i]] = 1;
//             }
//         }

//         for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
//         {
//             if (A[DG.inlist.nbors[i]] >= 1)
//             {
//                 q.push(make_pair(DG.inlist.nbors[i],level+1));
//                 A[DG.inlist.nbors[i]] = 0;
//                 // B[DG.outlist.nbors[i]] = 1;
//             }
//         }
//     }
//     return;
// }

// returns the sphere of v in B - with parents
bool find_l_sphere(CDAG &DG, VertexIdx v, vector<int> &B, int max_cycle_size, vector<int> &levels)
{
    VertexIdx n = DG.outlist.nVertices;
    bool found_cycle = false;
    static vector<int> A(n);
    static vector<VertexIdx> ancestors(n, -1);
    static vector<VertexIdx> parents(n, -1); 
    A = B; // B is what is available coming in, A is a copy of B to start with, 
    // as we do the BFS we mark a done vertex with 1 in A so equivalent to blocking

    queue<pair<VertexIdx, int>> q; 
   
    A[v] = 1;
    levels[v] = 1;
    ancestors[v] = v;
    parents[v] = v;
   
    for (int i=DG.outlist.offsets[v]; i<(int)DG.outlist.offsets[v+1]; i++)
    {
        if (A[DG.outlist.nbors[i]] == 0)
        {
            q.push(make_pair(DG.outlist.nbors[i],2));
            levels[DG.outlist.nbors[i]] = 2;
            A[DG.outlist.nbors[i]] = 1;
            ancestors[DG.outlist.nbors[i]] = DG.outlist.nbors[i];
            parents[DG.outlist.nbors[i]] = DG.outlist.nbors[i];
        }
    }

    while(!q.empty())
    {
        pair<VertexIdx,int> p = q.front();
        int level = p.second;
        if (level > ((max_cycle_size /2)+1))
        {
            return found_cycle; 
        }

        VertexIdx ver = p.first;
        q.pop();
        // B[ver] = 1;
        levels[ver] = level;

        for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
        {
            VertexIdx nbr = DG.outlist.nbors[i];
            if (nbr <= v) continue; 

            if ((A[nbr] == 1) // nbr has been seen by the BFS before
                && (B[nbr] == 0) // nbr was originally available. this means it cannot be v
                && (parents[ver]!= nbr) // nbr is not the immediate parent of ver
                && (!found_cycle) // cycle has not yet been found so it is worth checking
                && ((levels[ver] != 2) || (levels[nbr] != 2)) // it is not the case that both are of level 2. at least one is of level >2
                && (levels[nbr] < max_cycle_size)
                && (levels[nbr] > 1) 
                && (nbr > v) // we only need to do BFS with vertices > v
                // thus, here we have seen nbr before and it is not the immediate parent of ver, nor v, and at least 1 of them is not level 2
                // && (parents[ver] <)
                && (ancestors[ver] != ancestors[nbr]) // ver and nbr don't have the same parent
                && (levels[ver] + levels[nbr] -1 <= max_cycle_size) // the edge (ver,nbr) closes a cycle that is within the size range
                )
            {
                found_cycle = true;
            }
            else if (A[nbr] == 0)
            {
                q.push(make_pair(nbr,level+1));
                A[nbr] = 1;
                ancestors[nbr] = ancestors[ver];
                parents[nbr] = ver;
                levels[nbr] = level+1;
            }
        }

        for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
        {
            VertexIdx nbr = DG.inlist.nbors[i];
            if (nbr <= v) continue; 
            if ((A[nbr] == 1) // nbr has been seen by the BFS before
                && (B[nbr] == 0) // nbr was originally available. this means it cannot be v
                && (parents[ver]!= nbr) // nbr is not the immediate parent of ver
                && (!found_cycle) // cycle has not yet been found so it is worth checking
                && ((levels[ver] != 2) || (levels[nbr] != 2)) // it is not the case that both are of level 2. at least one is of level >2
                && (levels[nbr] < max_cycle_size)
                && (levels[nbr] > 1) 
                && (nbr > v) // we only need to do BFS with vertices > v
                // thus, here we have seen nbr before and it is not the immediate parent of ver, nor v, and at least 1 of them is not level 2
                && (ancestors[ver] != ancestors[nbr]) // ver and nbr don't have the same parent
                && (levels[ver] + levels[nbr] -1 <= max_cycle_size) // the edge (ver,nbr) closes a cycle that is within the size range
                )
            {
                found_cycle = true;                
            }
            else if (A[nbr] == 0)
            {
                q.push(make_pair(nbr,level+1));
                A[nbr] = 1;
                ancestors[nbr] = ancestors[ver];
                parents[nbr] = ver;
                levels[nbr] = level+1;
            }
        }
    }
    return found_cycle;
}



// returns the sphere of v in B
bool isConnectedShortestPath(CDAG &DG, VertexIdx src, VertexIdx dest, vector<int> A, int max_path_size)
{
    queue<pair<VertexIdx, int>> q; 
    q.push(make_pair(src,1));
    A[src] = 0;
    VertexIdx n = DG.outlist.nVertices;
    // vector<int> B = new vector<int>(n,0);

    for (int i=DG.outlist.offsets[src]; i<(int)DG.outlist.offsets[src+1]; i++)
    {
        if (A[DG.outlist.nbors[i]] == 1)
        {
            q.push(make_pair(DG.outlist.nbors[i],2));
            A[DG.outlist.nbors[i]] = 0;
            // B[DG.outlist.nbors[i]] = 1;
        }
    }

    while(!q.empty())
    {
        pair<VertexIdx,int> p = q.front();
        int level = p.second;
        if (level > max_path_size) return false;
        VertexIdx ver = p.first;
        q.pop();
        A[ver] = 1;

        if (ver == dest) return true;

        for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
        {
            if (A[DG.outlist.nbors[i]] >= 1)
            {
                q.push(make_pair(DG.outlist.nbors[i],level+1));
                A[DG.outlist.nbors[i]] = 0;
                // B[DG.outlist.nbors[i]] = 1;
            }
        }

        for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
        {
            if (A[DG.inlist.nbors[i]] >= 1)
            {
                q.push(make_pair(DG.inlist.nbors[i],level+1));
                A[DG.inlist.nbors[i]] = 0;
                // B[DG.outlist.nbors[i]] = 1;
            }
        }
    }
    return false;
}

bool isConnectedBFSOld(CDAG &DG, VertexIdx src, VertexIdx dest, VertexIdx v, vector<int> A)
{
    queue<VertexIdx> q; 
    q.push(src);
    A[src] = 0;
    VertexIdx n = DG.outlist.nVertices;

    // cout << "In isConnectedBFSOld. src = " << src << " dest = " << dest << " v = " << v << endl;
    // vector<VertexIdx> toUnblock;
    // toUnblock.reserve(NNBORS);
    // bool rFound = false;
    while(!q.empty())
    {
        VertexIdx ver = q.front();
        q.pop();


        if ((DG.isEdge(dest,ver) != -1) || (ver == dest)) return true;

        for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
        {
            if (A[DG.outlist.nbors[i]] >= 1)
            {
                q.push(DG.outlist.nbors[i]);
                A[DG.outlist.nbors[i]] = 0;
            }
        }

        for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
        {
            if (A[DG.inlist.nbors[i]] >= 1)
            {
                q.push(DG.inlist.nbors[i]);
                A[DG.inlist.nbors[i]] = 0;
            }
        }

        // for (int i=v+1; i<(int)n; i++)
        // {
        //     if ((DG.isEdge(i,ver) != -1) && (A[i] == 1))
        //     {
        //         q.push(i);
        //         A[i] = 0;
        //     }
        // }
    }


        // for (int i=DG.outlist.offsets[ver]; i<(int)DG.outlist.offsets[ver+1]; i++)
        // {
        //     VertexIdx nbr = DG.outlist.nbors[i];
        //     if (nbr == dest)
        //     {
        //         return true;
        //     }
        //     if (A[nbr] == 1)
        //     {
        //         if (nbr < v) cout << "nbr < v. This shouldn't happen" << endl;
        //         q.push(nbr);
        //         A[nbr] = 0;
        //         // toUnblock.push_back(nbr);
        //     }
        // }
        // for (int i=DG.inlist.offsets[ver]; i<(int)DG.inlist.offsets[ver+1]; i++)
        // {
        //     VertexIdx nbr = DG.inlist.nbors[i];
        //     if (nbr == dest)
        //     {
        //         return true;
        //     }
        //     if (A[nbr] == 1)
        //     {
        //         if (nbr < v) cout << "nbr < v. This shouldn't happen" << endl;
        //         q.push(nbr);
        //         A[nbr] = 0;
        //         // toUnblock.push_back(nbr);
        //     }
        // }
    // }
    // for (int i=0; i<toUnblock.size(); i++)
    // {
    //     blocked[toUnblock[i]] = 0;
    // }
    // if (rFound == true)
    //     return true;
    // else return false;
    return false;
}
#endif

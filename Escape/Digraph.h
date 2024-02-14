#ifndef ESCAPE_DIGRAPH_H_
#define ESCAPE_DIGRAPH_H_

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
// #include "Escape/ErrorCode.h"
#include <algorithm>
#include <map>
#include <vector>
#include <unordered_set>

using namespace std;
using namespace Escape;

static bool isBlankLine( const char* line )
{
  while( *line )
  {
    if( !isspace( *line ) )
      return false;
    ++line;
  }
  return true;
}

// DAG structure has two pointers, one to the 
// adjacency list of outedge, one to the adjacency list of inedges
struct CDAG
{
    CGraph outlist;
    CGraph inlist;

    // Checks if edge (v1, v2) is present
    int isEdge(VertexIdx v1, VertexIdx v2);
    int isEdgeBinary(VertexIdx v1, VertexIdx v2);
};

// Checks if edge (v1, v2) is present in CGraph
// If edge is present: Return index (in nbors) of v2 in neighbors of v1
// If edge is not present: Return -1
int CDAG::isEdge(VertexIdx v1, VertexIdx v2)
{
    if (v1 >= outlist.nVertices)
        return -1;
    for (EdgeIdx i=outlist.offsets[v1]; i < outlist.offsets[v1+1]; ++i)
        if (outlist.nbors[i] == v2)
            return i;

    if (v2 >= outlist.nVertices)
        return -1;
    for (EdgeIdx i=outlist.offsets[v2]; i < outlist.offsets[v2+1]; ++i)
        if (outlist.nbors[i] == v1)
            return i;

    return -1;
}

//Note: functionality is slightly different from the 
// other isEdgeBinary functions.
int CDAG::isEdgeBinary(VertexIdx v1, VertexIdx v2)
{
    if (v2 < v1)
    {
        VertexIdx swp = v1;
        v1 = v2;
        v2 = swp;
    }

    // can assume that v1 < v2 
    /* if (v1,v2) edge exists then v2 exists in v1's outnbrs
    and v1 exists in v2's innbrs
    */
    VertexIdx outdegv1 = outlist.offsets[v1+1] - outlist.offsets[v1];
    VertexIdx indegv2 = inlist.offsets[v2+1] - inlist.offsets[v2];

    if (outdegv1 < indegv2) // search in outnbrs of v1
    {
        EdgeIdx low = outlist.offsets[v1];
        EdgeIdx high = outlist.offsets[v1+1]-1;
        EdgeIdx mid;

        while(low <= high)
        {
            mid = (low+high)/2;

            if (outlist.nbors[mid] == v2)
                return 1;
            if (outlist.nbors[mid] > v2)
                high = mid-1;
            else
                low = mid+1;
        }
        return -1;
    }
    else
    {
        EdgeIdx low = inlist.offsets[v2];
        EdgeIdx high = inlist.offsets[v2+1]-1;
        EdgeIdx mid;

        while(low <= high)
        {
            mid = (low+high)/2;

            if (inlist.nbors[mid] == v1)
                return 1;
            if (inlist.nbors[mid] > v1)
                high = mid-1;
            else
                low = mid+1;
        }
        return -1;
    }

    // if(deg2 < deg1)
    // {
    //     VertexIdx swp = v1;
    //     v1 = v2;
    //     v2 = swp;
    // }

    
}

void delCDAG(CDAG DAG)
{
  delCGraph(DAG.outlist);
  delCGraph(DAG.inlist);
}

// Structure for comparing nodes according to their degree.
// So u < v if degree of u less than that of v in graph g.

struct DegreeComp
{
    CGraph *g;
    DegreeComp(CGraph *g) { this->g = g;}

    bool operator () (VertexIdx u, VertexIdx v)
    {
        VertexIdx degu = g->offsets[u+1] - g->offsets[u];  // Degree of u
        VertexIdx degv = g->offsets[v+1] - g->offsets[v];  // Degree of v
    
        if (degu < degv || (degu == degv && u < v))    // Comparing degrees and breaking ties by id
            return true;
        else
            return false;
    }
};

//Construct DAG based on degree ordering
//
// Input: Pointer for CGraph g
// Output: CDAG for degree ordering in g
//         This is the CDAG for the DAG where each edge points from lower degree endpoint to higher degree endpoint.
//
//
//         The outlist in CDAG is guaranteed to be sorted by degrees. This means that the neighbors
//         of every vertex in the outlist are sorted by their degrees in g. This is quite useful in
//         further processing.
CDAG degreeOrdered(CGraph *g)
{
    CDAG ret;     // CDAG to be returned
    printf("In degreeOrdered. number of edges= %ld", g->nEdges);
    CGraph outdag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};  // Initialize DAG of out-edges
    CGraph indag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};   // Initialize DAG of in-edges
    EdgeIdx outcur = 0;
    EdgeIdx incur = 0;
    VertexIdx dest;
    VertexIdx degi;
    VertexIdx degdest;

    outdag.offsets[0] = 0;
    indag.offsets[0] = 0;
    for (VertexIdx i=0; i < g->nVertices; ++i)   // Looping over all vertices in g
    {
        for (EdgeIdx j = g->offsets[i]; j < g->offsets[i+1]; ++j)   // Looping over neighbors of i in g
        {
            dest = g->nbors[j];     // We are now looking at edge (i,dest)
            degi = g->offsets[i+1] - g->offsets[i];   // Degree of i
            degdest = g->offsets[dest+1]- g->offsets[dest];   // Degree of dest
            //printf("i=%ld dest=%ld degi=%ld degdest=%ld\n",i,dest,degi,degdest);

            //We now orient the edge depending of degi vs degdest.
            // We break ties according to vertex id.
            // In the output, the g-edge (i,dest) is either pointing to dest (in if condition) or pointing to i (in else condition).

            if (degi < degdest || (degi == degdest && i < dest))   
            {
                outdag.nbors[outcur] = dest;   // We want point edge from i to dest. So this directed edge is added to outdag.
                ++outcur;                      // Increment pointer in outdag.nbors and the number of edges in outdag.
                ++outdag.nEdges;
            }
            else
            {
                indag.nbors[incur] = dest;     // We point edge from dest to i. So this edge goes into indag.
                ++incur;                       // Pointer and number of edges incremented
                ++indag.nEdges;
            }
        }
        outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
        indag.offsets[i+1] = incur;
    }

    for (VertexIdx i=0; i < g->nVertices;++i)  // Loops over vertices
        std::sort(outdag.nbors+outdag.offsets[i], outdag.nbors+outdag.offsets[i+1], DegreeComp(g)); // In outdag, sort all neighbors of i according to their degree. Note that DegreeComp gives the desired comparator.

    ret.outlist = outdag;
    ret.inlist = indag;

    return ret;
}


//Construct DAG based on degeneracy ordering. This only makes sense for an undirected graph
//
// Input: Pointer for CGraph g
// Output: CDAG for degeneracy ordering in g
//         This is the CDAG for the DAG where each edge points from lower order vertex endpoint to higher order endpoint according to the degree ordering.



// Note: This function is not doing what it is supposed to!
CDAG degeneracyOrdered(CGraph *g, VertexIdx *ordering)
{
    CDAG ret;     // CDAG to be returned
    CGraph outdag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};  // Initialize DAG of out-edges
    CGraph indag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};   // Initialize DAG of in-edges
    EdgeIdx outcur = 0;
    EdgeIdx incur = 0;
    VertexIdx dest;
    VertexIdx degi;
    VertexIdx degdest;

    outdag.offsets[0] = 0;
    indag.offsets[0] = 0;

    VertexIdx n = g->nVertices;
    VertexIdx min_deg = n;

    map <VertexIdx, unordered_set<VertexIdx> > deg_list;  
    map <VertexIdx, unordered_set<VertexIdx> > in_nbrs;  
    map <VertexIdx, unordered_set<VertexIdx> > out_nbrs;  
    map <VertexIdx, bool> touched;  
    vector<VertexIdx> cur_degs;
    cur_degs.resize(n);

    for (VertexIdx i=0; i<n; i++) // populate the map of deg -> list of vertices with that degree
    {
        VertexIdx deg = g->degree(i);
        deg_list[deg].insert(i);
        cur_degs[i] = deg;
        if (min_deg > deg) min_deg = deg;
    }

    for(VertexIdx i=0; i<n; i++)
    {
        while(deg_list[min_deg].size() == 0)
            min_deg++;
        unordered_set<VertexIdx>::iterator it = deg_list[min_deg].begin();
        VertexIdx source = *it; // current vertex
        touched[source] = true;
        ordering[i] = source;


        deg_list[min_deg].erase(it);
        
        for (int j=g->offsets[source]; j<g->offsets[source+1]; j++)
        {
            VertexIdx nbr = g->nbors[j];
            // if (touched[nbr] == true) // already processed. Add to in_nbrs of source
            if (touched.find(nbr) != touched.end()) // already processed. Add to in_nbrs of source
            {
                in_nbrs[source].insert(nbr);
                continue;
            }
            else
            {
                out_nbrs[source].insert(nbr);
                VertexIdx deg = cur_degs[nbr];
                deg_list[deg].erase(nbr);
                deg_list[deg-1].insert(nbr);
                if (deg-1 < min_deg)
                    min_deg = deg-1;
                cur_degs[nbr] = deg - 1;    
            }
        }
    }

    for(VertexIdx i=0; i<n; i++)
    {
        
        for (unordered_set<VertexIdx>::iterator in_it = in_nbrs[i].begin(); in_it!=in_nbrs[i].end(); ++in_it)
        {
            indag.nbors[incur] = *in_it;     // We point edge from dest to i. So this edge goes into indag.
            ++incur;                       // Pointer and number of edges incremented
            ++indag.nEdges;
        }
        for (unordered_set<VertexIdx>::iterator out_it = out_nbrs[i].begin(); out_it!=out_nbrs[i].end(); ++out_it)
        {
            outdag.nbors[outcur] = *out_it;     // We point edge from dest to i. So this edge goes into indag.
            ++outcur;                       // Pointer and number of edges incremented
            ++outdag.nEdges;
        }
        outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
        indag.offsets[i+1] = incur;
    }
    
    ret.outlist = outdag;
    ret.inlist = indag;

    return ret;
}

VertexIdx * getDegenWedgeOrdering(CGraph &cg)
{
    map <VertexIdx, unordered_set<VertexIdx> > induced_wedge_list;
    VertexIdx n = cg.nVertices;
    VertexIdx *num_induced_wedges = (VertexIdx *) calloc(n, sizeof(VertexIdx));
    VertexIdx *rank = (VertexIdx *) calloc(n, sizeof(VertexIdx));
    VertexIdx min_num_wedges = n*n;
    // map <VertexIdx, bool> touched; 
    vector<VertexIdx> touched(n,0); 
    vector<VertexIdx> cur_num_wedges(n,0);
    // cur_num_wedges.resize(n);

    for (VertexIdx i=0; i<n; i++)
    {
        for (EdgeIdx jpos=cg.offsets[i]; jpos < cg.offsets[i+1]; jpos++)
        {
            VertexIdx j = cg.nbors[jpos];
            for (VertexIdx kpos = cg.offsets[j]; kpos < cg.offsets[j+1]; kpos++)
            {
                VertexIdx k = cg.nbors[kpos];
                if ((i<k) && (cg.isEdgeBinary(i,k) == false))
                {
                    num_induced_wedges[j]++;
                }
            }
        }
    }

    for (VertexIdx i=0; i<n; i++)
    {
        induced_wedge_list[num_induced_wedges[i]].insert(i);

        if (min_num_wedges > num_induced_wedges[i])
            min_num_wedges = num_induced_wedges[i];
    }

    for(VertexIdx i=0; i<n; i++)
    {
        while(induced_wedge_list[min_num_wedges].size() == 0)
            min_num_wedges++;
        unordered_set<VertexIdx>::iterator it = induced_wedge_list[min_num_wedges].begin();
        VertexIdx source = *it; // current vertex
        touched[source] = 1;
        rank[source] = i;

        induced_wedge_list[min_num_wedges].erase(it);
        
        for (EdgeIdx jpos=cg.offsets[source]; jpos<cg.offsets[source+1]; jpos++)
        {
            VertexIdx j = cg.nbors[jpos];
            VertexIdx cur_wedges = num_induced_wedges[j];
            if (touched[j] == 1) // already processed. Add to in_nbrs of source
            {
                continue;
            }
            else
            {
                
                for (VertexIdx kpos = cg.offsets[j]; kpos < cg.offsets[j+1]; kpos++)
                {
                    VertexIdx k = cg.nbors[kpos];
                    if (touched[k] == 1)
                        continue;
                    else // k is not yet touched so update its count of wedges
                    {
                        if ((source != k) && (cg.isEdgeBinary(source,k) == false))
                        {
                            num_induced_wedges[j]--;
                        }
                    }
                }
            }

            induced_wedge_list[cur_wedges].erase(j);
            induced_wedge_list[num_induced_wedges[j]].insert(j);
            if (num_induced_wedges[j] < min_num_wedges)
                min_num_wedges = num_induced_wedges[j];  
        }
    }

    return rank;

}

VertexIdx * getReverseDegenWedgeOrdering(CGraph &cg)
{
    map <VertexIdx, unordered_set<VertexIdx> > induced_wedge_list;
    VertexIdx n = cg.nVertices;
    VertexIdx *num_induced_wedges = (VertexIdx *) calloc(n, sizeof(VertexIdx));
    VertexIdx *rank = (VertexIdx *) calloc(n, sizeof(VertexIdx));
    VertexIdx max_num_wedges = 0;
    vector<VertexIdx> touched(n,0); 
    vector<VertexIdx> cur_num_wedges(n,0);

    for (VertexIdx i=0; i<n; i++)
    {
        for (EdgeIdx jpos=cg.offsets[i]; jpos < cg.offsets[i+1]; jpos++)
        {
            VertexIdx j = cg.nbors[jpos];
            for (VertexIdx kpos = cg.offsets[j]; kpos < cg.offsets[j+1]; kpos++)
            {
                VertexIdx k = cg.nbors[kpos];
                if ((i<k) && (cg.isEdgeBinary(i,k) == false))
                {
                    num_induced_wedges[j]++;
                }
            }
        }
    }

    for (VertexIdx i=0; i<n; i++)
    {
        induced_wedge_list[num_induced_wedges[i]].insert(i);

        if (max_num_wedges < num_induced_wedges[i])
            max_num_wedges = num_induced_wedges[i];
    }

    for(VertexIdx i=0; i<n; i++)
    {
        while(induced_wedge_list[max_num_wedges].size() == 0)
            max_num_wedges--;
        unordered_set<VertexIdx>::iterator it = induced_wedge_list[max_num_wedges].begin();
        VertexIdx source = *it; // current vertex
        touched[source] = 1;
        rank[source] = i;

        induced_wedge_list[max_num_wedges].erase(it);
        
        for (EdgeIdx jpos=cg.offsets[source]; jpos<cg.offsets[source+1]; jpos++)
        {
            VertexIdx j = cg.nbors[jpos];
            VertexIdx cur_wedges = num_induced_wedges[j];
            if (touched[j] == 1) // already processed. Add to in_nbrs of source
            {
                continue;
            }
            else
            {
                
                for (VertexIdx kpos = cg.offsets[j]; kpos < cg.offsets[j+1]; kpos++)
                {
                    VertexIdx k = cg.nbors[kpos];
                    if (touched[k] == 1)
                        continue;
                    else // k is not yet touched so update its count of wedges
                    {
                        if ((source != k) && (cg.isEdgeBinary(source,k) == false))
                        {
                            num_induced_wedges[j]--;
                        }
                    }
                }
            }

            induced_wedge_list[cur_wedges].erase(j);
            induced_wedge_list[num_induced_wedges[j]].insert(j);
        }
    }

    return rank;

}



///////////////


#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

typedef struct {
    VertexIdx s;
    VertexIdx t;
} edge;

typedef struct {
    VertexIdx node;
    VertexIdx deg;
} nodedeg ;

typedef struct {
    VertexIdx n;//number of nodes
    VertexIdx e;//number of edges
    edge *edges;//list of edges
    VertexIdx *rank;//ranking of the nodes according to degeneracy ordering
    //VertexIdx *map;//oldID newID correspondance NOT USED IN THIS VERSION
} edgelist;

typedef struct {
    VertexIdx n;
    VertexIdx *cd;//cumulative degree: (starts with 0) length=n+1
    VertexIdx *adj;//truncated list of neighbors
    VertexIdx core;//core value of the graph
} graph;

typedef struct {
    VertexIdx *n;//n[l]: number of nodes in G_l
    VertexIdx **d;//d[l]: degrees of G_l
    VertexIdx *adj;//truncated list of neighbors
    unsigned char *lab;//lab[i] label of node i
    VertexIdx **nodes;//sub[l]: nodes in G_l
    VertexIdx core;
} subgraph;

void free_edgelist(edgelist *el){
    free(el->edges);
    free(el->rank);
    free(el);
}

void free_graph(graph *g){
    free(g->cd);
    free(g->adj);
    free(g);
}

void free_subgraph(subgraph *sg, unsigned char k){
    unsigned char i;
    free(sg->n);
    for (i=2;i<k;i++){
        free(sg->d[i]);
        free(sg->nodes[i]);
    }
    free(sg->d);
    free(sg->nodes);
    free(sg->lab);
    free(sg->adj);
    free(sg);
}


//Compute the maximum of three unsigned integers.
// inline unsigned int max3(unsigned int a,unsigned int b,unsigned int c){
inline VertexIdx max3(VertexIdx a,VertexIdx b,VertexIdx c){
    a=(a>b) ? a : b;
    return (a>c) ? a : c;
}

edgelist* readedgelist(const char* input){
    VertexIdx e1=NLINKS;
    edgelist *el= (edgelist *) malloc(sizeof(edgelist));
    FILE *file;
    char line[1024];

    el->n=0;
    el->e=0;
    file=fopen(input,"r");

    if (!file)
    {
        fprintf(stderr, "could not open file %s\n", input);
        // return ecInvalidInput;
    }
    el->edges= (edge *) malloc(e1*sizeof(edge));
    // while (fscanf(file,"%u %u", &(el->edges[el->e].s), &(el->edges[el->e].t))==2) {//Add one edge
    

    // while (fscanf(file,"%ld\t%ld", &(el->edges[el->e].s), &(el->edges[el->e].t))==2) {//Add one edge
    // fgets(line, sizeof(line), file);
    while (fgets(line, sizeof(line), file))
    {
    //Ignore comment lines.
        if ((line[0] == '#') || (line[0] == '%'))
            continue;
        else if (isBlankLine(line))
            continue;
        else
            break;
    }
    int64_t i1, i2;
    sscanf(line, "%ld%ld", &i1, &i2);
    while (fgets(line, sizeof(line), file))
    {
    //Ignore comment lines.
        if ((line[0] == '#') || (line[0] == '%'))
        continue;
        
      
        // int64_t i1, i2;
        sscanf(line, "%ld%ld", &i1, &i2);

        el->edges[el->e].s = i1;
        el->edges[el->e].t = i2;

        // cout << el->edges[el->e].s << "," << el->edges[el->e].t << endl;
        el->n=max3(el->n,el->edges[el->e].s,el->edges[el->e].t);
        el->e++;
        if (el->e==e1) {
            e1+=NLINKS;
            el->edges= (edge *) realloc(el->edges,e1*sizeof(edge));
        }
    }
    fclose(file);
    el->n++;

    el->edges= (edge *) realloc(el->edges,el->e*sizeof(edge));

    return el;
}

void relabel(edgelist *el){
    VertexIdx i, source, target, tmp;

    for (i=0;i<el->e;i++) {
        // cout << "Edge " << el->edges[i].s << "," << el->edges[i].t << endl;
        source=el->n - 1 - el->rank[el->edges[i].s];
        target=el->n - 1 - el->rank[el->edges[i].t];
        if (source>target){
            tmp=source;
            source=target;
            target=tmp;
        }
        el->edges[i].s=source;
        el->edges[i].t=target;
    }

}

///// CORE ordering /////////////////////

typedef struct {
    VertexIdx key;
    VertexIdx value;
} keyvalue;

typedef struct {
    VertexIdx n_max; // max number of nodes.
    VertexIdx n; // number of nodes.
    VertexIdx *pt;   // pointers to nodes.
    keyvalue *kv; // nodes.
} bheap;


bheap *construct(VertexIdx n_max){
    VertexIdx i;
    bheap *heap= (bheap *) malloc(sizeof(bheap));

    heap->n_max=n_max;
    heap->n=0;
    heap->pt= (VertexIdx *) malloc(n_max*sizeof(VertexIdx));
    for (i=0;i<n_max;i++) heap->pt[i]=-1;
    heap->kv=(keyvalue *) malloc(n_max*sizeof(keyvalue));
    return heap;
}

void swap(bheap *heap,VertexIdx i, VertexIdx j) {
    keyvalue kv_tmp=heap->kv[i];
    VertexIdx pt_tmp=heap->pt[kv_tmp.key];
    heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
    heap->kv[i]=heap->kv[j];
    heap->pt[heap->kv[j].key]=pt_tmp;
    heap->kv[j]=kv_tmp;
}

void bubble_up(bheap *heap,VertexIdx i) {
    VertexIdx j=(i-1)/2;
    while (i>0) {
        if (heap->kv[j].value>heap->kv[i].value) {
            swap(heap,i,j);
            i=j;
            j=(i-1)/2;
        }
        else break;
    }
}

void bubble_down(bheap *heap) {
    VertexIdx i=0,j1=1,j2=2,j;
    while (j1<heap->n) {
        j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
        if (heap->kv[j].value < heap->kv[i].value) {
            swap(heap,i,j);
            i=j;
            j1=2*i+1;
            j2=j1+1;
            continue;
        }
        break;
    }
}

void insert(bheap *heap,keyvalue kv){
    heap->pt[kv.key]=(heap->n)++;
    heap->kv[heap->n-1]=kv;
    bubble_up(heap,heap->n-1);
}

void update(bheap *heap,VertexIdx key){
    VertexIdx i=heap->pt[key];
    if (i!=-1){
        ((heap->kv[i]).value)--;
        bubble_up(heap,i);
    }
}

void update_max_heap(bheap *heap,VertexIdx key){
    VertexIdx i=heap->pt[key];
    if (i!=-1){
        ((heap->kv[i]).value)++;
        bubble_up(heap,i);
    }
}

keyvalue popmin(bheap *heap){
    keyvalue min=heap->kv[0];
    heap->pt[min.key]=-1;
    heap->kv[0]=heap->kv[--(heap->n)];
    heap->pt[heap->kv[0].key]=0;
    bubble_down(heap);
    return min;
}

// keyvalue popmax(bheap *heap){
//     keyvalue max=heap->kv[0];
//     heap->pt[max.key]=-1;
//     heap->kv[0]=heap->kv[--(heap->n)];
//     heap->pt[heap->kv[0].key]=0;
//     bubble_down(heap);
//     return max;
// }

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(VertexIdx n,VertexIdx *v){
    VertexIdx i;
    keyvalue kv;
    bheap* heap=construct(n);
    for (i=0;i<n;i++){
        kv.key=i;
        kv.value=v[i];
        insert(heap,kv);
    }
    return heap;
}

bheap* mkmaxheap(VertexIdx n,VertexIdx *v){
    VertexIdx i;
    keyvalue kv;
    bheap* heap=construct(n);
    for (i=0;i<n;i++){
        kv.key=i;
        kv.value=-1*v[i];
        insert(heap,kv);
    }
    return heap;
}

void freeheap(bheap *heap){
    free(heap->pt);
    free(heap->kv);
    free(heap);
}

//computing degeneracy ordering and core value
void ord_core(edgelist* el){
    VertexIdx i,j,r=0,n=el->n,e=el->e;
    keyvalue kv;
    bheap *heap;

    VertexIdx *d0= (VertexIdx *) calloc(el->n,sizeof(VertexIdx));
    VertexIdx *cd0= (VertexIdx *) malloc((el->n+1)*sizeof(VertexIdx));
    VertexIdx *adj0= (VertexIdx *) malloc(2*el->e*sizeof(VertexIdx));
    for (i=0;i<e;i++) {
        d0[el->edges[i].s]++;
        d0[el->edges[i].t]++;
    }
    cd0[0]=0;
    for (i=1;i<n+1;i++) {
        cd0[i]=cd0[i-1]+d0[i-1];
        d0[i-1]=0;
    }
    for (i=0;i<e;i++) {
        adj0[ cd0[el->edges[i].s] + d0[ el->edges[i].s ]++ ]=el->edges[i].t;
        adj0[ cd0[el->edges[i].t] + d0[ el->edges[i].t ]++ ]=el->edges[i].s;
    }

    heap=mkheap(n,d0);

    el->rank= (VertexIdx *) malloc(n*sizeof(VertexIdx));
    for (i=0;i<n;i++){
        kv=popmin(heap);
        el->rank[kv.key]=n-(++r);
        for (j=cd0[kv.key];j<cd0[kv.key+1];j++){
            update(heap,adj0[j]);
        }
    }
    freeheap(heap);
    free(d0);
    free(cd0);
    free(adj0);
}

// void ord_induced_wedges(edgelist* el){
//     VertexIdx i,j,r=0,n=el->n,e=el->e;
//     keyvalue kv;
//     bheap *heap;

//     VertexIdx *d0= (VertexIdx *) calloc(el->n,sizeof(VertexIdx));
//     VertexIdx *cd0= (VertexIdx *) malloc((el->n+1)*sizeof(VertexIdx));
//     VertexIdx *adj0= (VertexIdx *) malloc(2*el->e*sizeof(VertexIdx));
//     for (i=0;i<e;i++) {
//         d0[el->edges[i].s]++;
//         d0[el->edges[i].t]++;
//     }
//     //d0 has the degrees of each vertex
//     //cd0 stores teh cumulative degrees
//     cd0[0]=0;
//     for (i=1;i<n+1;i++) {
//         cd0[i]=cd0[i-1]+d0[i-1];
//         d0[i-1]=0;
//     }
//     for (i=0;i<e;i++) {
//         adj0[ cd0[el->edges[i].s] + d0[ el->edges[i].s ]++ ]=el->edges[i].t;
//         adj0[ cd0[el->edges[i].t] + d0[ el->edges[i].t ]++ ]=el->edges[i].s;
//     }

//     heap=mkheap(n,d0);

//     el->rank= (VertexIdx *) malloc(n*sizeof(VertexIdx));
//     for (i=0;i<n;i++){
//         kv=popmin(heap);
//         el->rank[kv.key]=n-(++r);
//         for (j=cd0[kv.key];j<cd0[kv.key+1];j++){
//             update(heap,adj0[j]);
//         }
//     }
//     freeheap(heap);
//     free(d0);
//     free(cd0);
//     free(adj0);
// }

void ord_degree(edgelist* el){
    VertexIdx i,j,r=0,n=el->n,e=el->e;
    keyvalue kv;
    bheap *heap;

    VertexIdx *d0= (VertexIdx *) calloc(el->n,sizeof(VertexIdx));
    VertexIdx *cd0= (VertexIdx *) malloc((el->n+1)*sizeof(VertexIdx));
    VertexIdx *adj0= (VertexIdx *) malloc(2*el->e*sizeof(VertexIdx));
    for (i=0;i<e;i++) {
        d0[el->edges[i].s]++;
        d0[el->edges[i].t]++;
    }
    cd0[0]=0;
    for (i=1;i<n+1;i++) {
        cd0[i]=cd0[i-1]+d0[i-1];
        d0[i-1]=0;
    }
    for (i=0;i<e;i++) {
        adj0[ cd0[el->edges[i].s] + d0[ el->edges[i].s ]++ ]=el->edges[i].t;
        adj0[ cd0[el->edges[i].t] + d0[ el->edges[i].t ]++ ]=el->edges[i].s;
    }

    heap=mkmaxheap(n,d0);

    el->rank= (VertexIdx *) malloc(n*sizeof(VertexIdx));
    for (i=0;i<n;i++){
        kv=popmin(heap);
        el->rank[kv.key]=n-(++r);
        for (j=cd0[kv.key];j<cd0[kv.key+1];j++){
            update_max_heap(heap,adj0[j]);
        }
    }
    freeheap(heap);
    free(d0);
    free(cd0);
    free(adj0);
}


//////////////////////////

graph* mkgraph(edgelist *el){
    VertexIdx i,max;
    VertexIdx *d;
    graph* g = (graph *) malloc(sizeof(graph));

    d= (VertexIdx *) calloc(el->n,sizeof(VertexIdx));

    for (i=0;i<el->e;i++) {
        d[el->edges[i].s]++;
    }

    g->cd = (VertexIdx *) malloc((el->n+1)*sizeof(VertexIdx));
    g->cd[0]=0;
    max=0;
    for (i=1;i<el->n+1;i++) {
        g->cd[i]=g->cd[i-1]+d[i-1];
        max=(max>d[i-1])?max:d[i-1];
        d[i-1]=0;
    }
    // printf("core value (max truncated degree) = %u\n",max);

    g->adj= (VertexIdx *) malloc(el->e*sizeof(VertexIdx));

    for (i=0;i<el->e;i++) {
        g->adj[ g->cd[el->edges[i].s] + d[ el->edges[i].s ]++ ]=el->edges[i].t;
    }

    free(d);
    g->core=max;
    g->n=el->n;
    return g;
}



CDAG convertEdgeListToCDAG(edgelist *el, VertexIdx *ordering)
{
    graph *g = mkgraph(el);


    for (int i=0;i<el->n;i++) 
    {
        // ordering[el->rank[i]] = i;
        ordering[i] = el->n - 1 - i;
        // ordering[i] = i;
    }

    CDAG ret;     // CDAG to be returned
    CGraph outdag = {el->n, el->e, new EdgeIdx[el->n+1], new VertexIdx[el->e+1]};  // Initialize DAG of out-edges
    outdag.offsets = g->cd;
    outdag.nbors = g->adj;     // We point edge from dest to i. So this edge goes into indag.
//             ++outcur;                       // Pointer and number of edges incremented
//             ++outdag.nEdges;
//         }
//         outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
//         indag.offsets[i+1] = incur;
//     }
    
    CGraph indag = {el->n, el->e, new EdgeIdx[el->n+1], new VertexIdx[el->e+1]};  // Initialize DAG of out-edges

  
    edgelist *el_in= (edgelist *) malloc(sizeof(edgelist));
    el_in->n = el->n;
    el_in->e = el->e;
    el_in->edges= (edge *) malloc(el->e*sizeof(edge));

    for (int i=0; i<el->e; i++)
    {
        // cout << i << endl;
        el_in->edges[i].s = el->edges[i].t;
        el_in->edges[i].t = el->edges[i].s;
    }

    graph *g_in = mkgraph(el_in);
    indag.offsets = g_in->cd;
    indag.nbors = g_in->adj;

    ret.inlist = indag;
    ret.inlist.nVertices = el->n;
    ret.inlist.nEdges = el->e;

    ret.outlist = outdag;
    ret.outlist.nVertices = el->n;
    ret.outlist.nEdges = el->e;

    ret.outlist.sortById();
    ret.inlist.sortById();

//     ret.inlist = indag;

    return ret;
}






    // EdgeIdx outcur = 0;
  
    // outdag.offsets[0] = 0;

    // for (i=0;i<el->e;i++) {
    //     source=el->rank[el->edges[i].s];
    //     target=el->rank[el->edges[i].t];
    //     if (source<target){
    //         tmp=source;
    //         source=target;
    //         target=tmp;
    //     }
    //     el->edges[i].s=source;
    //     el->edges[i].t=target;
    // }

    /* 

    ret.outlist = outdag;

    return ret;
}*/
#endif




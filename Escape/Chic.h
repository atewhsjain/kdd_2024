#ifndef ESCAPE_CYCLES_H_
#define ESCAPE_CYCLES_H_

#include "Escape/CycleHelper.h"
#include <cassert>
#include <vector>
#include <random>
#include <ctgmath>
#include <iomanip>
#include <algorithm>

using namespace Escape;
using namespace std;

bool DEBUG = false;


void extend_path(CDAG &DG, VertexIdx l, VertexIdx v, vector<VertexIdx> &I, vector<int> &A, vector<int> &active_nbrs_v, vector<int> &levels, int num_active_neighbors_of_v, vector<long long> &counts_cycles, vector<long long> &counts_paths, long long &num_rec_calls, int cycle_size, bool &flag_vertex_cycle_found)
{
	num_rec_calls++;
	// cout << "length of path = " << I.size() + 2 << endl;

	counts_paths[I.size()+2]++;

	if (I.size()+2 > cycle_size) 
		return;

	if ((levels[l] < cycle_size) && (I.size()+ levels[l] -1 > cycle_size)) 
		return;


	assert(A[l] == 1);
	
	VertexIdx n = DG.outlist.nVertices;


	if (num_active_neighbors_of_v <= 0) return;

	for (VertexIdx i=DG.outlist.offsets[l]; i<(int)DG.outlist.offsets[l+1]; i++)
	{
		VertexIdx ver = DG.outlist.nbors[i];
		A[ver]++;
	}

	for (VertexIdx i=DG.inlist.offsets[l]; i<(int)DG.inlist.offsets[l+1]; i++)
	{
		VertexIdx ver = DG.inlist.nbors[i];
		A[ver]++;
	}


	I.push_back(l);

	for (VertexIdx i=DG.outlist.offsets[l]; i<(int)DG.outlist.offsets[l+1]; i++)
	{
		VertexIdx ver = DG.outlist.nbors[i];
		if (A[ver]==1)
		{
			if (active_nbrs_v[ver] == 2)
			{
				counts_cycles[(int)I.size()+1]++; //x and l
				num_active_neighbors_of_v--;
				flag_vertex_cycle_found = true;
				continue;
			}
			else if ((levels[ver] + I.size() -1 <= cycle_size) && (num_active_neighbors_of_v > 0))
				extend_path(DG, ver, v, I, A, active_nbrs_v, levels, num_active_neighbors_of_v, counts_cycles, counts_paths, num_rec_calls, cycle_size, flag_vertex_cycle_found);

		}
	}

	for (VertexIdx i=DG.inlist.offsets[l]; i<(int)DG.inlist.offsets[l+1]; i++)
	{
		VertexIdx ver = DG.inlist.nbors[i];
		if (A[ver]==1)
		{
			// B[ver] = 1;
			if (active_nbrs_v[ver] == 2)
			{
				counts_cycles[(int)I.size()+1]++; //x and l
				num_active_neighbors_of_v--;
				flag_vertex_cycle_found = true;
				continue;
			}
			else if ((levels[ver] + I.size() -1 <= cycle_size) && (num_active_neighbors_of_v > 0))
				extend_path(DG, ver, v, I, A, active_nbrs_v, levels, num_active_neighbors_of_v, counts_cycles, counts_paths, num_rec_calls, cycle_size, flag_vertex_cycle_found);

		}
	}

	for (VertexIdx i=DG.outlist.offsets[l]; i<(int)DG.outlist.offsets[l+1]; i++)
	{
		VertexIdx ver = DG.outlist.nbors[i];
		A[ver]--;
	}

	for (VertexIdx i=DG.inlist.offsets[l]; i<(int)DG.inlist.offsets[l+1]; i++)
	{
		VertexIdx ver = DG.inlist.nbors[i];
		A[ver]--;
	}
	I.pop_back();
}

void chic(CDAG &DG, string gname, int cycle_size)
{

	high_resolution_clock::time_point cyclecount_s = high_resolution_clock::now();
	VertexIdx n = DG.outlist.nVertices;
	// vector<int> A(n,1);
	vector<int> A(n,0);
	vector<int> levels;
	levels.reserve(n);
	vector<int> active_nbrs_v(n,0);
	vector<long long> counts_cycles(n+1,0);
	vector<long long> counts_paths(n+1,0);
	long long num_rec_calls = 0;

	counts_paths[1] = n;
	vector<VertexIdx> I;
	I.reserve(cycle_size);
	int num_vertices_cycle_found = 0;
	int num_vertices_active_outnbrs = 0;
	bool found_cycle = true;
	int num_vertices_found_cycle_bfs = 0;

	for (VertexIdx v=0; v<n; v++)
	{
		A[v] = 1;
		levels.assign(n,10000);
		bool flag_vertex_cycle_found = false;
		bool flag_vertex_active_outnbrs = false;

		found_cycle = find_l_sphere(DG, v, A, cycle_size, levels);
		if (found_cycle) 
		{
			num_vertices_found_cycle_bfs++;
		}

		I.push_back(v);
		
		for (VertexIdx i=DG.outlist.offsets[v] ; i<DG.outlist.offsets[v+1]; i++)
		{
	        active_nbrs_v[DG.outlist.nbors[i]] = 2;
		}

		for (VertexIdx i=DG.outlist.offsets[v] ; i<DG.outlist.offsets[v+1]; i++)
	    {
	    	VertexIdx u = DG.outlist.nbors[i];
	    	A[u] = 1;

	    	int num_active_neighbors_of_v = DG.outlist.offsets[v+1] - i - 1;
	    	for (VertexIdx j=i+1 ; j<DG.outlist.offsets[v+1]; j++)
		    {
		    	VertexIdx w = DG.outlist.nbors[j];
		    	if (DG.isEdgeBinary(u,w) != -1)
		    	{
		    		// cout << "Found cycle of length 3\n";
		    		counts_cycles[3]++;
		    		A[w] = 1;
		    		num_active_neighbors_of_v--;
		    		continue;
		    	}
		    }
		    if (num_active_neighbors_of_v > 0) flag_vertex_active_outnbrs = true;
		    if ((num_active_neighbors_of_v > 0) && (found_cycle))
		    {
			    extend_path(DG, u, v, I, A, active_nbrs_v, levels, num_active_neighbors_of_v, counts_cycles, counts_paths, num_rec_calls, cycle_size, flag_vertex_cycle_found);
		    }
		    for (VertexIdx j=i+1 ; j<DG.outlist.offsets[v+1]; j++)
		    {
		    	VertexIdx w = DG.outlist.nbors[j];
		    	A[w] = 0;
		    }
		}
		if (flag_vertex_active_outnbrs) num_vertices_active_outnbrs++;
		if (flag_vertex_cycle_found) num_vertices_cycle_found++;


		for (VertexIdx i=DG.outlist.offsets[v] ; i<DG.outlist.offsets[v+1]; i++)
		{
	        A[DG.outlist.nbors[i]] = 0;
	        active_nbrs_v[DG.outlist.nbors[i]] = 0;
		}
		I.pop_back();
	}

	high_resolution_clock::time_point cyclecount_e = high_resolution_clock::now();

	auto duration_cyclecount = std::chrono::duration_cast<std::chrono::microseconds>( cyclecount_e - cyclecount_s ).count();

	cout << "Output files:" << endl;
	string pathsdatafile(gname + "_" + to_string(cycle_size) + "_pathcounts_chic");
    std::string pdfname = "../results/" + pathsdatafile;
    std::cout << pdfname << std::endl;
    std::ofstream pdf;
    pdf.open(pdfname);
    if (!pdf.is_open())
    {
        std::cout << "Could not open output file." << std::endl;
        exit(1);
    }
    pdf << "size;count" << endl;

    string cyclesdatafile(gname + "_" + to_string(cycle_size) + "_cyclecounts_chic");
    std::string cdfname = "../results/" + cyclesdatafile;
    std::cout << cdfname << std::endl;
    std::ofstream cdf;
    cdf.open(cdfname);
    if (!cdf.is_open())
    {
        std::cout << "Could not open output file." << std::endl;
        exit(1);
    }
    cdf << "size;count" << endl;

    string paramsfile(gname + "_" + to_string(cycle_size) + "_meta_chic");
    std::string pfname = "../results/" + paramsfile;
    std::cout << pfname << std::endl;
    std::ofstream paramsf;
    paramsf.open(pfname);
    if (!paramsf.is_open())
    {
        std::cout << "Could not open output file." << std::endl;
        exit(1);
    }

    cout << "\nCycle counts:" << endl;
    for (VertexIdx i=0; i<n+1; i++)
	{
		if (counts_cycles[i] != 0)
		{
			cout << "size = " << i << " count = " << counts_cycles[i] << endl;
			cdf << i << ";" << counts_cycles[i] << endl;
		}
		if (counts_paths[i] != 0)
		{
			pdf << i << ";" << counts_paths[i] << endl;
		}
	}

	cout << endl;
    paramsf << "time;num_rec_calls;num_vertices_active_outnbrs;num_vertices_cycle_found;num_vertices_active_outnbrs_but_no_cycle_found" << endl;
    paramsf << duration_cyclecount << ";" << num_rec_calls << ";" << num_vertices_active_outnbrs << ";" << num_vertices_cycle_found << ";" << num_vertices_active_outnbrs - num_vertices_cycle_found << endl;
    cout << "Time = " << duration_cyclecount << " milliseconds" << endl;
	cout << "Number of recursive calls = " << num_rec_calls << endl;
	cout << endl;

	// cout << "\n\nnum_vertices_active_outnbrs = " << num_vertices_active_outnbrs << endl;
	// cout << "\n\nnum_vertices_active_outnbrs_but_no_cycle_found = " << num_vertices_active_outnbrs - num_vertices_cycle_found << endl;
	// // cout << "\n\nnum_vertices that have active outnbrs but lead to no cycles = " << num_vertices_no_cycle - num_vertices_no_active_outnbrs << endl;
	// cout << "\n\nnum_vertices that lead to cycles = " <<  num_vertices_cycle_found << endl;


	paramsf.close();
	pdf.close();
	cdf.close();

}

#endif
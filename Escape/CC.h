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


void extend_path(CDAG &DG, VertexIdx l, VertexIdx v, vector<VertexIdx> &I, VertexIdx r, vector<int> &A, vector<long long> &counts_cycles, vector<long long> &counts_paths, long long &num_rec_calls, int cycle_size, bool &flag_vertex_cycle_found)
{
	num_rec_calls++;

	counts_paths[I.size()+2]++;

	if (I.size()+2 >= cycle_size) 
		return;


	assert(A[l] == 1);
	assert(A[r] == 1);
	
	VertexIdx n = DG.outlist.nVertices;

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
			if (DG.isEdgeBinary(ver,r) != -1)
			{
				counts_cycles[(int)I.size()+2]++; //x and l
				flag_vertex_cycle_found = true;
				
				continue;
			}
			extend_path(DG, ver, v, I, r, A, counts_cycles, counts_paths, num_rec_calls, cycle_size, flag_vertex_cycle_found);
		}
	}

	for (VertexIdx i=DG.inlist.offsets[l]; i<(int)DG.inlist.offsets[l+1]; i++)
	{
		VertexIdx ver = DG.inlist.nbors[i];
		if (A[ver]==1)
		{
			if (DG.isEdgeBinary(ver,r) != -1)
			{
				counts_cycles[(int)I.size()+2]++; //x and l
				flag_vertex_cycle_found = true;
				
				continue;
			}
			extend_path(DG, ver, v, I, r, A, counts_cycles, counts_paths, num_rec_calls, cycle_size, flag_vertex_cycle_found);

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

void cc(CDAG &DG, string gname, int cycle_size)
{

	high_resolution_clock::time_point cyclecount_s = high_resolution_clock::now();
	VertexIdx n = DG.outlist.nVertices;
	vector<int> A(n,0);
	vector<long long> counts_cycles(n+1,0);
	vector<long long> counts_paths(n+1,0);
	long long num_rec_calls = 0;
	int num_vertices_cycle_found = 0;
	int num_vertices_active_outnbrs = 0;
	

	counts_paths[1] = n;
	vector<VertexIdx> I;
	I.reserve(cycle_size);

	for (VertexIdx v=0; v<n; v++)
	{
		A[v] = 1;
		I.push_back(v);
		bool flag_vertex_cycle_found = false;
		bool flag_vertex_active_outnbrs = false;
		

		for (VertexIdx i=DG.outlist.offsets[v] ; i<DG.outlist.offsets[v+1]; i++)
	        A[DG.outlist.nbors[i]] = 1;

		for (VertexIdx i=DG.outlist.offsets[v] ; i<DG.outlist.offsets[v+1]; i++)
	    {
	    	VertexIdx u = DG.outlist.nbors[i];
	    	for (VertexIdx j=i+1 ; j<DG.outlist.offsets[v+1]; j++)
		    {
		    	VertexIdx w = DG.outlist.nbors[j];
		    	if (DG.isEdgeBinary(u,w) != -1)
		    	{
		    		// cout << "Found cycle of length 3\n";
		    		counts_cycles[3]++;
		    		continue;
		    	}
		    	flag_vertex_active_outnbrs = true;
		    	extend_path(DG, u, v, I, w, A, counts_cycles, counts_paths, num_rec_calls, cycle_size, flag_vertex_cycle_found);
		    }
		}

		if (flag_vertex_active_outnbrs) num_vertices_active_outnbrs++;
		if (flag_vertex_cycle_found) num_vertices_cycle_found++;


		for (VertexIdx i=DG.outlist.offsets[v] ; i<DG.outlist.offsets[v+1]; i++)
	        A[DG.outlist.nbors[i]] = 0;
		I.pop_back();
	}

	high_resolution_clock::time_point cyclecount_e = high_resolution_clock::now();

	auto duration_cyclecount = std::chrono::duration_cast<std::chrono::microseconds>( cyclecount_e - cyclecount_s ).count();

	
	cout << "Output files:" << endl;
	string pathsdatafile(gname + "_" + to_string(cycle_size) + "_pathcounts_cc");
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

    string cyclesdatafile(gname + "_" + to_string(cycle_size) + "_cyclecounts_cc");
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

    string paramsfile(gname + "_" + to_string(cycle_size) + "_meta_cc");
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
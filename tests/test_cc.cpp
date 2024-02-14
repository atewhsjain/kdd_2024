#include <iostream>
#include <fstream>
#include <array>
#include <unistd.h>
#include <cstring>
#include "Escape/GraphIO.h"
#include "Escape/Digraph.h"
#include "Escape/CC.h"


using namespace Escape;

void usage() {
    printf("\n./test_cc \t-i input_file_1[,input_file_2,...] \t-k cycle_size  \n");
    printf("\nEg: ./test_cc \n\t-i oz-oz.edges \n\t-k 10  \n");
}

int main(int argc, char *argv[])
{

    if(argc <5) {
        usage();
        exit(1);
    }

    // std::ofstream of; // output file
    std::vector<std::string> graphs; // input files
	srand (time(NULL));
    int cycle_size = 0;
    char c;
    printf("*****************************\n");
    while ((c = getopt(argc, argv, "i:k:")) != -1) {
        switch (c) {
            case 'i': {
                // multiple input files
                char *token = std::strtok(optarg, ",");
                printf("Input Files: ");
                while (token != NULL) {
                    graphs.push_back(std::string(token));
                    printf("%s,",token);
                    token = std::strtok(NULL, ",");
                }
                printf("\n");
            }
            break;

            case 'k':
                cycle_size = atoi(optarg);
                cout << "Cycle size: " << cycle_size << endl;
            break;
		}
    }
    printf("*****************************\n\n");
    
    
    for (std::vector<std::string>::iterator it=graphs.begin(); it!=graphs.end(); ++it)
    {
        std::string fname = "../graphs/" + *it;
        edgelist* el=readedgelist(fname.c_str());
        printf("Number of nodes = %u\n",el->n);
        printf("Number of edges = %u\n",el->e);


        printf("Building the graph structure\n");
        ord_core(el);
        // ord_degree(el);
        relabel(el);
        
        VertexIdx *ordering = (VertexIdx *) calloc(el->n,sizeof(VertexIdx)); 
        printf("Creating DAG\n");
        CDAG dag = convertEdgeListToCDAG(el, ordering);

        free_edgelist(el);

        string f = *it;
        char *token = std::strtok(strdup(f.c_str()), ".");
        string gname(token);
        cc(dag, gname, cycle_size);
    }

}

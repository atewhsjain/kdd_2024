This code accompanies the paper titled "A fast and efficient algorithm for counting cycles". 

The "graphs" folder stores all the graphs used in the experiments in the paper. The graphs are stored in a custom format with the extension ".edges". Each ".edges" file consists of a list of pairs of numbers on each line separated by a space. The first number on the first line represents the number of vertices in the graph and the 2nd number on the 1st line represents the number of edges. Every subsequent line represents an edge in the graph. All these graphs have been obtained from Konect (http://konect.cc/) and from SNAP (https://snap.stanford.edu/data/index.html). The graphs have been made simple and undirected.

The result files are stored in the "results" folder. This code uses clang.

To compile the code:
cd to the directory that has this README.
From this directory, run 

	make clean
	make

cd into the "tests" directory.
From the "tests" directory, run

	make clean
	make

*********Chic**********
To run the Chic algorithm (proposed in the paper), run

	./test_chic -i input_file_1 -k cycle_size 

For eg:
	./test_chic -i oz-oz.edges -k 6

The code expects that the "oz-oz.edges" file must be in the "graphs" folder (which is at the same level as "tests").
The above command outputs on the command line the number of cycles found by the Chic algorithm in the oz-oz.edges graph for each size from 3 to 6 (cycle_size), the time required to find these cycles in milliseconds and the number of recursive calls made by the algorithm.


*********CC**********
To run the CC algorithm (proposed in the paper), run

	./test_cc -i input_file_1 -k cycle_size 

For eg:
	./test_cc -i oz-oz.edges -k 6

The code expects that the "oz-oz.edges" file must be in the "graphs" folder (which is at the same level as "tests").
The above command outputs on the command line the number of cycles found by the CC algorithm in the oz-oz.edges graph for each size from 3 to 6 (cycle_size), the time required to find these cycles in milliseconds and the number of recursive calls made by the algorithm.




Our Paper Title:
  This code is the official implementation of my paper titled 
    “A Novel Approach for Constrained Optimization in Graphical Models”

Requirements for Compilation:
	— Access to SCIP, Gurobi or CBC Mixed Integer Linear Programming (MILP) solvers.
		Scip: https://www.scipopt.org/
		Gurobi: https://www.gurobi.com/
		CBC: https://projects.coin-or.org/Cbc
	- Google OR tools with C++ interface
		https://developers.google.com/optimization
	- To compile the code, you can use the provided CMakeLists.txt file as a reference
	- The code has four executables and main*.cpp files associated with the executables
		(1) CMPE: optimization algorithm (Algorithm-CMPE) described in the paper
		(2) generate_mps: Convert CMPE to MILP format for use by MILP solvers
		(3) generate_q: Generate "q" values used in the paper
		(4) MCKP_Greedy: Test code for Greedy solver for MCKP
	- If you prefer compiling it via commandline, use the following commands:
		To compile CMPE use the following command:
		 - g++ -O3 -std=c++11 main.cpp Knapsack.cpp MN.cpp -o CMPE
		To compile generate_mps use the following two commands:
		 - g++ -c -I <ortools-include-dir> -O3 -std=c++11 main_generate_mps.cpp\
			MN.cpp Knapsack.cpp  -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP
		 - g++ -o generate_mps -L <ortools-lib-dir> main_generate_mps.o MN.o\
			Knapsack.o -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc\
			 -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools
		To compile generate_q use the following command:
		 - g++ -O3 -std=c++11 main_generate_q.cpp Knapsack.cpp MN.cpp -o generate_q
		To compile MCKP_Greedy use the following two commands:
		 - g++ -c -I <ortools-include-dir> -O3 -std=c++11 main_mckp_test.cpp\
			Knapsack.cpp  -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP
		 - g++ -o generate_mps -L <ortools-lib-dir> main_mckp_test.o\
			Knapsack.o -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc\
			 -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools

Requirements for Running the Algorithm:
	- Once compiled, you can download the Markov network files from the UAI 2014 competition
		http://www.hlt.utdallas.edu/~vgogate/uai14-competition/
	- After compiling you can get help using: "<executable> -h". For example,
		./CMPE -h
	   will generate the following output
		Usage: ./CMPE -m1 <uaifilename1> -m2 <uaifilename2> -o <outfilename> -q <q-value>
	 	Other Options: [-t -k -s -si ]
		-------------------------------------------------------------------------
		 Details on Required Option
			 uaifilename1 and uaifilename2: are evidence instantiated Markov networks in UAI format
			 outfilename: Results of experiments will be stored here
			 q-value: (Real number): constraint on weight of the assignment in CMPE
		-------------------------------------------------------------------------
		 Details on Other Options and Default values
			 -t     [int]: max time for which each k is run; default 2
			 -k     [int]: max k-seperator. The code iterates from 1 to k in increments of 2
			 -s     [int]: Seed for Repeatability; default 1000000L
			 -si    [int]: print status every integer seconds; default 1
		
Results and Evaluation:
To get results described in the paper, run the following two commands:
	./generate_q -m <uaifilename> -o <q-filename>
	./CMPE -m1 <uaifilename1> -m2 <uaifilename2> -q <one-q-from-q-file> -o <stats-filename> -t 1200 -k 11
	The stats-file contains the following information in comma separated format which can be
	plotted using any plotting software.
	
		q-used,k-used,current-value-of-objective-function,time-in-seconds
	
	For example, the two types of plots in the paper are:
		(1) For each k: value of objective function as a function of time for a given q
		(2) For each k: value of objective function as a function of q for a given time
	
	Note: If uaifilename1=uaifilename2 then we get subset-sum type of problems.
	Otherwise, we get Knapsack-type problems.

	To generate mps-files that can be used for MILP solvers, use:
	./generate_mps -m1 <uaifilename1> -m2 <uaifilename2> -q <one-q-from-q-file> -o <mps-filename>

Code is released under:
    MIT License.







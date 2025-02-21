SCHEDULER = 8

SpMV:
	g++ -std=c++14 -Ofast -I include main.cpp -D SPMV -D SCHEDULER='$(SCHEDULER)' -lm -fopenmp -g -o dist/Sim

SpMSpV:
	g++ -std=c++14 -Ofast -I include main.cpp -D SPMSpV -D SCHEDULER='$(SCHEDULER)' -lm -fopenmp -g -o dist/Sim

SpMM:
	g++ -std=c++14 -Ofast -I include main.cpp -D SPMM -D SCHEDULER='$(SCHEDULER)' -lm -fopenmp -g -o dist/Sim

SpGEMM:
	g++ -std=c++14 -Ofast -I include main.cpp -D SPGEMM -D SCHEDULER='$(SCHEDULER)' -lm -fopenmp -g -o dist/Sim

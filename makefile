hellomake: psort.cpp
	mpic++ psort.cpp -o psort.o
	ar cr libpsort.a psort.o

my:psort.cpp a3.cpp
	mpic++ -c psort.cpp -o psort.o
	mpic++ -c a3.cpp -o a3.o
	ar cr libpsort.a psort.o
	mpic++ a3.o libpsort.a -o a.out
	# time mpirun -np 4 ./a.out
	mpirun -np 4 ./a.out
	cat output_dir/pivoted_0.txt output_dir/ipvt_0.txt
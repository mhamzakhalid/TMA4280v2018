# Run 1D mesh generation
mesh: $(1dmesh.cpp)
	mpiCC 1dmesh.cpp -lm -o mesh
	mpirun -np 4 ./mesh

#Run 1D for assembly
assembly: $(assembly.cpp)
	mpiCC assembly.cpp -lm -o matrix
	mpirun -np 4 ./matrix


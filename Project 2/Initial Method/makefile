# Run 1D case
1D: $(1d.cpp)
	mpiCC 1d.cpp -lm -o 1d
	mpirun -np 4 ./1d

# Run 2D case
2D: $(2d.cpp)
	mpiCC 2d.cpp -lm -o 2d
	mpirun -np 4 ./2d

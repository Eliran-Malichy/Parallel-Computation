build:
	mpicxx -fopenmp -c main.c -o main.o
	mpicxx -fopenmp -c Result.c -o Result.o
	mpicxx -fopenmp -c cFunctions.c -o cFunctions.o
	nvcc -Xcompiler -fopenmp -I./inc -c cudaFunctions.cu -o cudaFunctions.o
	mpicxx -fopenmp -lgomp -o finalProject  main.o Result.o cFunctions.o cudaFunctions.o  /usr/local/cuda/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./mpiCudaOpemMP

run:
	mpiexec -np 3 ./finalProject < input.txt > output.txt

runH:
	mpiexec -hostfile hosts -np 3 ./finalProject < input.txt > output.txt

run1:
	mpiexec -np 1 ./finalProject < input.txt > output.txt

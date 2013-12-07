#include <curand_kernel.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define BLOCK_SIZE 16 
#define E_CONST 2.71828182845904523536f

//! Perform the Ising simulation
//! @param lattice A pointer to the lattice of atoms.
//! @param height The number of rows of the input lattice.
//! @param width The number of columns of the input lattice.
//! @param T The temperature, in units of epsilon/k. Epsilon is the exchange energy and k is the boltzmann constant.
//! #param iterations The number of Metropolis iterations to perform.
__global__ void ising(int * lattice, int height, int width, float T, long iterations, unsigned long long seed) {
	__shared__ int slattice[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ int sneighbors[4 * BLOCK_SIZE];

	int tx = blockIdx.x * BLOCK_SIZE + threadIdx.x;
	int ty = blockIdx.y * BLOCK_SIZE + threadIdx.y;
	int deltaU, top, bottom, left, right;
	curandStatePhilox4_32_10_t state;

	//Initialize random number generator
	/* From the curand library guide:
	   Each experiment should be assigned a unique seed.
	   Each thread should have a unique sequenc number. 
	*/
	//TODO: It may be beneficial to separate calls to curand_init() and curand() into separate kernels for maximum performance.
	//See 3.5: Performance notes of curand library guide
	curand_init(seed, ty * width + tx, 0, &state); //seed, sequence, offset

		//Perform simulation
		//Each turn of the loop performs BLOCK_SIZE^2 iterations of the Metropolis algorithm
		for(int k = 0; k < iterations; k++) {

	//Load sublattice into shared memory
	if(tx < width && ty < height) {
		slattice[threadIdx.y][threadIdx.x] = lattice[ty * width + tx];
		if(threadIdx.y == 0) {
			if(ty == 0)
				sneighbors[threadIdx.x] = lattice[(height - 1) * width + tx];
			else
				sneighbors[threadIdx.x] = lattice[(ty - 1) * width + tx];
		}

		if(threadIdx.y == (BLOCK_SIZE - 1)) {
			if(ty == (height - 1))
				sneighbors[2 * BLOCK_SIZE + threadIdx.x] = lattice[tx];
			else
				sneighbors[2 * BLOCK_SIZE + threadIdx.x] = lattice[(ty + 1) * width + tx];
		}

		if(threadIdx.x == 0) {
			if(tx == 0)
				sneighbors[3 * BLOCK_SIZE + threadIdx.y] = lattice[ty * width + (width - 1)];
			else
				sneighbors[3 * BLOCK_SIZE + threadIdx.y] = lattice[ty * width + (tx - 1)];
		}

		if(threadIdx.x == (BLOCK_SIZE - 1)) {
			if(tx == (width - 1))
				sneighbors[BLOCK_SIZE + threadIdx.y] = lattice[ty * width];
			else
				sneighbors[BLOCK_SIZE + threadIdx.y] = lattice[ty * width + (tx + 1)];
		}
		__syncthreads();

			for(int i = 0; i < 2; i ++) {
				//Checkerboard
				if((threadIdx.x + threadIdx.y) % 2 == i) {
					if(threadIdx.y == 0)
						top = sneighbors[threadIdx.x]; 
					else
						top = slattice[threadIdx.y - 1][threadIdx.x];
					
					if(threadIdx.x == 0)
						left = sneighbors[3 * BLOCK_SIZE + threadIdx.y];
					else
						left = slattice[threadIdx.y][threadIdx.x - 1];

					if(threadIdx.y == (BLOCK_SIZE - 1))
						bottom = sneighbors[2 * BLOCK_SIZE + threadIdx.x];
					else
						bottom = slattice[threadIdx.y + 1][threadIdx.x];

					if(threadIdx.x == (BLOCK_SIZE - 1))
						right = sneighbors[BLOCK_SIZE + threadIdx.y]; 
					else
						right = slattice[threadIdx.y][threadIdx.x + 1];

					//Calculate change in energy if dipole were flipped
					deltaU = 2 * slattice[threadIdx.y][threadIdx.x] * (top + bottom + left + right);

					//If the energy would decrease, flip the dipole
					if(deltaU <= 0)
						slattice[threadIdx.y][threadIdx.x] *= -1;
					else {
						float rand = curand_uniform(&state);
						//Else the probability of a flip is given by the Boltzmann factor
						if(rand < powf(E_CONST, -deltaU/T))
							slattice[threadIdx.y][threadIdx.x] *= -1;
					}
					__syncthreads(); //Is this needed?
				}
			}
		}
	lattice[ty * width + tx] = slattice[threadIdx.y][threadIdx.x]; 
	__syncthreads();
		}

}

//! Generate an image file of the lattice and write it to disk
//! @param lattice A 2D array of ints, whose values are either 1 or -1.
//! @param height The number of rows of the lattice
//! @param width The number of colums of the lattice
//! @param filename The name of the file to write to.
void print(int * lattice, int height, int width, char * filename) {

}

int main(int argc, char ** argv) {
	//The height and width of the lattice
	int height = 400;
	int width = 400;

	int * lattice = (int *)malloc(sizeof(int) * height * width);

	//Seed random number generator for kernel call
	srand(time(NULL));

	if(argc < 3) {
                printf("Usage: %s [iterations] [temperature]\n", argv[0]);
                return 0;
        }
	long n = strtol(argv[1], NULL, 0);
	float T = strtof(argv[2], NULL);

	//Initialize input lattice
	for(int i = 0; i < (height * width); i++) {
		//lattice[i] = (rand() % 2 ? 1 : -1);
		lattice[i] = 1; 
	}

	int * lattice_d = NULL;
	cudaMalloc((void **)&lattice_d, sizeof(int) * height * width);
	cudaMemcpy(lattice_d, lattice, height * width * sizeof(int), cudaMemcpyHostToDevice);

	dim3 blockDim(BLOCK_SIZE, BLOCK_SIZE, 1);
	dim3 gridDim(ceil(height/BLOCK_SIZE), ceil(width/BLOCK_SIZE), 1);
	unsigned long long seed = rand();
	ising<<<gridDim, blockDim>>>(lattice_d, height, width, T, n, seed);
	
	cudaMemcpy(lattice, lattice_d, height * width * sizeof(int), cudaMemcpyDeviceToHost);

	float magnetization = 0;
	for(int i = 0; i < (height * width); i++) {
		magnetization += lattice[i];
	}
	magnetization /= (height * width);

	/*
	for(int i = 0; i < height; i++) {
		for(int j = 0; j < width; j++) {
			printf("%d ", lattice[i * width + j]);
		}
		printf("\n");
	}
	*/
	printf("%f\n", magnetization);
	return 0;
}

#define BLOCK_SIZE 16

//! Perform the Ising simulation
//! @param lattice A pointer to the lattice of atoms.
//! @param height The number of rows of the input lattice.
//! @param width The number of columns of the input lattice.
//! @param T The temperature, in units of epsilon/k. Epsilon is the exchange energy and k is the boltzmann constant
__global__ void ising(int * lattice, int height, int width, float T) {
	__shared__ int slattice[BLOCK_SIZE][BLOCK_SIZE];
	__constant__ int neighbors[4 * BLOCK_SIZE];

	int tx = blockIdx.x * BLOCK_SIZE + threadIdx.x;
	int ty = blockIdx.y * BLOCK_SIZE + threadIdx.y;
	int deltaU, top, bottom, left, right;

	//Load sublattice into shared memory
	if(tx < width && ty < height) {
		slattice[threadIdx.y][threadIdx.x] = lattice[ty * width + threadIdx.x];
		//TODO: Figure out how to use constant memory
		//Maybe put the entire lattice into constant memory?
		#if 0
		if(threadIdx.y == 0) {
			if(ty == 0)
				neighbors[threadIdx.x] = lattice[(height - 1) * width + tx];
			else
				neighbors[threadIdx.x] = lattice[(ty - 1) * width + tx];
		}
		#endif
		__syncthreads();

	//Perform simulation
	for(int i = 0; i < 2; i ++) {
		//Checkerboard
		if((threadIdx.x + threadIdx.y) % 2 == i) {
			if(threadIdx.y == 0)
				top =  
			else
				top = slattice[threadIdx.y - 1][threadIdx.x];
			
			if(threadIdx.x == 0)
				left = 
			else
				left = slattice[threadIdx.y][threadIdx.x - 1];

			if(threadIdx.y == (BLOCK_SIZE - 1))
				bottom =
			else
				bottom = slattice[threadIdx.y + 1][threadIdx.x];

			if(threadIdx.x == (BLOCK_SIZE - 1))
				right =
			else
				right = slattice[threadIdx.y][threadIdx.x + 1];

			__syncthreads(); //Is this needed?

			//Calculate change in energy if dipole were flipped
			deltaU = 2 * slattice[threadIdx.y][threadIdx.x] * (top + bottom + left + right);

			//If the energy would decrease, flip the dipole
			if(deltaU <= 0)
				slattice[threadIdx.y][threadIdx.x] *= -1;
			else {
				rand = curand_uniform(randState);
				//TODO: Put in a more precise value of e
				//Else the probability of a flip is given by the Boltzmann factor
				if(rand < powf(2.71f, -deltaU/T))
					slattice[threadIdx.y][threadIdx.x] *= -1;
			}
		}
	}
	}

}



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h> // Must be compiled with -lrt flag
#include <math.h>

//Forward declaration
void simulate(int *, long, float, int);

int main(int argc, char ** argv) {
	//The side length of the lattice
	int length;
	double total_time;
	timespec ts, te;
	
	if(argc < 3) {
		printf("Usage: %s [iterations] [side length]\n", argv[0]);
		return 0;
	}

	long long unsigned iterations = strtol(argv[1], NULL, 0);
	length = (int)strtol(argv[2], NULL, 0);

	iterations = iterations * length * length; // Use "GPU Iterations"

	int * lattice = (int *)malloc(sizeof(int) * length * length);


	for(float T = 0.0f; T < 3.0f; T += 0.1) {
		//Initialize input lattice
		for(int i = 0; i < (length * length); i++) {
			lattice[i] = 1;
		}
		
		//Monte carlo simulation
		clock_gettime(CLOCK_REALTIME, &ts);
		simulate(lattice, iterations, T, length);
		clock_gettime(CLOCK_REALTIME, &te);

		total_time += (((double)te.tv_sec) - ((double)ts.tv_sec));

		float magnetization = 0;
		for(int i = 0; i < (length * length); i++) {
			magnetization += lattice[i];
		}

		magnetization /= (length * length);
		printf("%f\t%f\n", T, magnetization);
	}

	printf("Average Execution Time: %.2f s\n", total_time/30);
	return 0;
}

void simulate(int * lattice, long iterations, float T, int length) {
        int row, col, deltaU, top, bottom, left, right;
	long i = 0;
        float E = 2.718281828459045235f; //E = 2.7182...

	//Seed random nuber generator
        srand(time(NULL));

	while (i != iterations) {
            row = rand() % length; //Pick a random dipole
            col = rand() % length;

            //If the dipole is at the edge, pick the neighbor on the other side
            if (row == 0) {
                top = lattice[length * (length - 1) + col];
            } else {
                top = lattice[length * (row - 1) + col];
            }

            if (col == 0) {
                left = lattice[length * row + (length - 1)];
            } else {
                left = lattice[length * row + (col - 1)];
            }

            if (row == length - 1) {
                bottom = lattice[length * 0 + col];
            } else {
                bottom = lattice[length * (row + 1) + col];
            }

            if (col == length - 1) {
                right = lattice[length * row + 0];
            } else {
                right = lattice[length * row + (col + 1)];
            }

            /*
             * Calculates the change in energy if the dipole were flipped. If
             * two dipoles are pointing in the same direction, they have energy
             * -ε. If they are pointing in opposite directions (antiparallel),
             * they have energy ε.
             */
            deltaU = 2 * lattice[length * row + col] * (top + bottom + left + right);
            //If the energy would decrease, flip the dipole
            if (deltaU <= 0) {
                lattice[length * row + col] *= -1;
            } //Else the probability of a flip is exp[-U/kT]
            else if (((float)rand())/RAND_MAX < powf(E, -(float)deltaU/T)) {
                lattice[length * row + col] *= -1;
            }
            i++;
        }
}

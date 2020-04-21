#include "common.h"
#include <omp.h>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

typedef struct
{   int* indeces;
    int bin_size;
    int capacity;
    int nei_num;
    int* neighbor;
}bin_t;

static int bin_i, bin_j;
static bin_t *bin_list; 
static double bin_x, bin_y;
static int num_bins;
static omp_lock_t * bin_locks;
static omp_lock_t * part_locks;


void init_grid(int num_bins, bin_t *bin_list)
{
	int x, y, cur_x, cur_y, cur_id;
	int dx[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
    int dy[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
    for(int i = 0; i < num_bins; i ++)
    {
        bin_list[i].capacity = 20;
        bin_list[i].bin_size = 0;
        bin_list[i].indeces = (int*)malloc(bin_list[i].capacity*sizeof(int));
        bin_list[i].nei_num = 0;
        bin_list[i].neighbor = (int*)malloc(9*sizeof(int));
        x = i % bin_i;
        y = i / bin_i;
        for(int j = 0; j < 9; j++){
        	cur_x = x + dx[j];
        	cur_y = y + dy[j];
        	if(cur_x >= 0 && cur_x < bin_i && cur_y >=0 && cur_y < bin_j){
        		cur_id = cur_y * bin_i + cur_x;
        		bin_list[i].neighbor[bin_list[i].nei_num] = cur_id;
        		bin_list[i].nei_num++;
        	}
        }
    }
}

void add_particle(bin_t *bin_list, int i, int j)
{ 
    if (bin_list[j].bin_size == bin_list[j].capacity)
    {
        bin_list[j].indeces = (int*) realloc(bin_list[j].indeces, 2*bin_list[j].capacity*sizeof(int));
        bin_list[j].capacity  *= 2;
    }

    bin_list[j].indeces[bin_list[j].bin_size] = i;
    bin_list[j].bin_size ++;
}

void bin_particles(int n, particle_t *particles, int num_bins, bin_t *bin_list, double bin_x, double bin_y, int num_rows)
{
    int i;
    for(i=0;i<num_bins;i++) bin_list[i].bin_size = 0;
    for(i=0;i<n;i++)
    {   
        int x = particles[i].x / bin_x, y = particles[i].y / bin_y;
        int index = y+x*num_rows;
        add_particle(bin_list, i, index);
    }
}

void remove_particle(bin_t *bin_list, int i, int j)
{
    for (int k = 0; k < bin_list[j].bin_size; k++)
    {
        if (bin_list[j].indeces[k] == i) 
        {
            for (int l = k; l < bin_list[j].bin_size; l ++)
                bin_list[j].indeces[l] = bin_list[j].indeces[l+1];
            bin_list[j].bin_size--;
            break;
        }
    }
    
}
// Put any static global variables here that you will use throughout the simulation.
// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
    neighbor.ax += -coef * dx;
    neighbor.ay += -coef * dy;
}


// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

void init_simulation(particle_t* particles, int n, double grid_size) {
	// You can use this space to initialize data objects that you may need
	// This function will be called once before the algorithm begins
	// Do not do any particle simulation here
	bin_i = (int)ceil(grid_size / cutoff);
    bin_j = (int)ceil(grid_size / cutoff);
    num_bins = bin_i * bin_j;
    bin_list = (bin_t*) malloc(num_bins * sizeof(bin_t));

    bin_x = grid_size / bin_i, bin_y = grid_size / bin_j;

    init_grid(num_bins, bin_list);
    bin_particles(n, particles, num_bins, bin_list, bin_x, bin_y, bin_j);

    bin_locks = new omp_lock_t[num_bins];
    for(int i = 0; i < num_bins; ++i)
        omp_init_lock(&bin_locks[i]);

    part_locks = new omp_lock_t[n];
    for(int i = 0; i < n; i++){
    	omp_init_lock(&part_locks[i]);
    }
}

void simulate_one_step(particle_t* particles, int n, double size){

	#pragma omp for
	for(int i = 0; i < n; i++){
		particles[i].ax = particles[i].ay = 0;
	}

	#pragma omp for
	for(int i = 0; i < num_bins; i++){
		bin_t* cur_bin = bin_list + i;
		bin_t* nei_bin;
		int cur_p, cur_n;
		for(int j = 0; j < cur_bin->bin_size; j++){
			for(int k = 0; k < cur_bin->nei_num; k++){
				nei_bin = bin_list + cur_bin->neighbor[k];
				for(int h = 0; h < nei_bin->bin_size; h++){
					cur_p = cur_bin->indeces[j];
					cur_n = nei_bin->indeces[h];
					if(cur_p > cur_n){
						omp_set_lock(&part_locks[cur_p]);
						omp_set_lock(&part_locks[cur_n]);
						apply_force(particles[cur_p], particles[cur_n]);
						omp_unset_lock(&part_locks[cur_n]);
						omp_unset_lock(&part_locks[cur_p]);
					}
				}
			}
		}
	}

	#pragma omp for
	for(int i = 0; i < n; i++){
		int r_old = particles[i].y / bin_y, c_old = particles[i].x / bin_x;
		move( particles[i] ,size);
		int r = particles[i].y / bin_y, c = particles[i].x / bin_x;
		if (r != r_old || c != c_old)
        {
        	omp_set_lock(&bin_locks[r_old + c_old * bin_j]);
            remove_particle(bin_list, i, r_old + c_old*bin_j);
            omp_unset_lock(&bin_locks[r_old + c_old * bin_j]);

            omp_set_lock(&bin_locks[r + c * bin_j]);
            add_particle(bin_list, i, r + c*bin_j);
            omp_unset_lock(&bin_locks[r + c * bin_j]);
        }
	}
	
} 

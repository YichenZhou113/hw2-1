#include "common.h"
#include <cmath>
#include <vector>

class bin_t {
  public:
    std::vector<particle_t*> particles;
    bin_t() { particles.resize(100); }
};

bin_t* bins;
int num_bins;
double bin_size;
int num_bins_per_side;

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
}

void apply_force_for_particle_in_bin(particle_t& particle, int bin_id) {
    int bin_x = bin_id / num_bins_per_side;
    int bin_y = bin_id % num_bins_per_side;

    particle.ax = particle.ay = 0;
    for (int j = 0; j < bins[bin_id].particles.size(); j++) {
        apply_force(particle, *(bins[bin_id].particles[j]));
    }

    if (bin_x != 0) {
        int neighbor_id = bin_id - num_bins_per_side;
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
    if (bin_x != num_bins_per_side - 1) {
        int neighbor_id = bin_id + num_bins_per_side;
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
    if (bin_y != 0) {
        int neighbor_id = bin_id - 1;
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
    if (bin_y != num_bins_per_side - 1) {
        int neighbor_id = bin_id + 1;
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
    if (bin_x != 0 && bin_y != 0) {
        int neighbor_id = bin_id - (num_bins_per_side + 1);
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
    if (bin_x != 0 && bin_y != num_bins_per_side - 1) {
        int neighbor_id = bin_id - (num_bins_per_side - 1);
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
    if (bin_x != num_bins_per_side - 1 && bin_y != 0) {
        int neighbor_id = bin_id + num_bins_per_side - 1;
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
    if (bin_x != num_bins_per_side - 1 && bin_y != num_bins_per_side - 1) {
        int neighbor_id = bin_id + num_bins_per_side + 1;
        for (int j = 0; j < bins[neighbor_id].particles.size(); j++) {
            apply_force(particle, *(bins[neighbor_id].particles[j]));
        }
    }
}

void bin_particles(particle_t* particles, int n) {
    for (int i = 0; i < num_bins_per_side; i++) {
        for (int j = 0; j < num_bins_per_side; j++) {
            bins[i * num_bins_per_side + j].particles.clear();
        }
    }
    for (int i = 0; i < n; i++) {
        int bin_x = floor(particles[i].x / bin_size);
        int bin_y = floor(particles[i].y / bin_size);
        int bin_id = bin_x * num_bins_per_side + bin_y;
        bins[bin_id].particles.push_back(particles + i);
    }
}

void apply_force_in_bin(int bin_id) {
    for (int i = 0; i < bins[bin_id].particles.size(); i++) {
        apply_force_for_particle_in_bin(*(bins[bin_id].particles[i]), bin_id);
    }
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

void init_simulation(particle_t* parts, int num_parts, double size) {
    // You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
    bin_size = cutoff;
    num_bins_per_side = ceil(size / bin_size);
    num_bins = num_bins_per_side * num_bins_per_side;
    bins = new bin_t[num_bins];
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {

    bin_particles(parts, num_parts);

    for (int i = 0; i < num_bins; i++) {
        apply_force_in_bin(i);
    }

    for (int i = 0; i < num_parts; i++)
        move(parts[i], size);
}
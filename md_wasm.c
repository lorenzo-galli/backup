#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <emscripten/emscripten.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Costanti fisiche di default (Argon: Angstrom, ps, u, K) */
#define DEFAULT_EPSILON 99.365319
#define DEFAULT_SIGMA 3.4
#define DEFAULT_MASS 39.9
#define KB 0.831446
#define MAX_FORCE 8000.0

typedef struct {
    double x;
    double y;
} Vec2D;

typedef struct {
    int N;
    double density;
    double L;
    double T_target;
    double dt;
    int mode_nvt;
    
    /* Parametri fisici della specie atomica selezionata */
    double sigma;
    double epsilon;
    double mass;
    
    Vec2D *pos;
    Vec2D *vel;
    Vec2D *force;
    
    int M;
    double cell_size;
    int **head;
    int *list;
    
    double time;
    int step;
    
    double E_kin;
    double E_pot;
    double E_tot;
    double T_inst;
} System;

static System *g_sys = NULL;

#define GR_NBINS 200
static double g_r_histogram[GR_NBINS];
static double g_r_axis[GR_NBINS];
static double g_r_result[GR_NBINS];
static int g_gr_active = 0;           /* 0: idle, 1: sampling, 2: completed */
static int g_gr_target_steps = 0;
static int g_gr_collected_steps = 0;
static double g_gr_dr = 0.0;

static double rand_gaussian(double mean, double stddev) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 < 1e-10) u1 = 1e-10;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mean + z0 * stddev;
}

static inline void wrap_positions(System *sys) {
    for (int i = 0; i < sys->N; i++) {
        while (sys->pos[i].x < 0.0) sys->pos[i].x += sys->L;
        while (sys->pos[i].x >= sys->L) sys->pos[i].x -= sys->L;
        while (sys->pos[i].y < 0.0) sys->pos[i].y += sys->L;
        while (sys->pos[i].y >= sys->L) sys->pos[i].y -= sys->L;
    }
}

static void build_cell_list(System *sys) {
    for (int i = 0; i < sys->M; i++) {
        for (int j = 0; j < sys->M; j++) {
            sys->head[i][j] = -1;
        }
    }
    
    for (int i = 0; i < sys->N; i++) {
        double px = sys->pos[i].x;
        double py = sys->pos[i].y;
        
        while (px < 0) px += sys->L;
        while (px >= sys->L) px -= sys->L;
        while (py < 0) py += sys->L;
        while (py >= sys->L) py -= sys->L;
        
        int cx = (int)(px / sys->cell_size);
        int cy = (int)(py / sys->cell_size);
        
        if (cx < 0) cx = 0;
        if (cx >= sys->M) cx = sys->M - 1;
        if (cy < 0) cy = 0;
        if (cy >= sys->M) cy = sys->M - 1;
        
        sys->list[i] = sys->head[cx][cy];
        sys->head[cx][cy] = i;
    }
}

/* Calcolo delle Forze Ottimizzato O(N) con Linked Cells e Cutoff Distance */
static void compute_forces_cell(System *sys) {
    for (int i = 0; i < sys->N; i++) {
        sys->force[i].x = 0.0;
        sys->force[i].y = 0.0;
    }
    sys->E_pot = 0.0;
    
    build_cell_list(sys);
    
    double r2_cut = (2.5 * sys->sigma) * (2.5 * sys->sigma);
    if (r2_cut > (sys->L * 0.5) * (sys->L * 0.5)) {
        r2_cut = (sys->L * 0.5) * (sys->L * 0.5);
    }
    
    double r2_min = (0.65 * sys->sigma) * (0.65 * sys->sigma);
    double sig2 = sys->sigma * sys->sigma;
    
    /* 9 celle vicine (inclusa la cella corrente) */
    const int neighbor_dx[9] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
    const int neighbor_dy[9] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
    
    for (int i = 0; i < sys->N; i++) {
        double px = sys->pos[i].x;
        double py = sys->pos[i].y;
        
        int cx = (int)(px / sys->cell_size);
        int cy = (int)(py / sys->cell_size);
        if (cx < 0) cx = 0; if (cx >= sys->M) cx = sys->M - 1;
        if (cy < 0) cy = 0; if (cy >= sys->M) cy = sys->M - 1;
        
        for (int n = 0; n < 9; n++) {
            int ncx = (cx + neighbor_dx[n] + sys->M) % sys->M;
            int ncy = (cy + neighbor_dy[n] + sys->M) % sys->M;
            
            int j = sys->head[ncx][ncy];
            while (j != -1) {
                if (j > i) {
                    double rx = px - sys->pos[j].x;
                    double ry = py - sys->pos[j].y;
                    
                    rx -= sys->L * round(rx / sys->L);
                    ry -= sys->L * round(ry / sys->L);
                    
                    double r2 = rx * rx + ry * ry;
                    if (r2 < r2_cut) {
                        double r2_safe = (r2 < r2_min) ? r2_min : r2;
                        double sr2 = sig2 / r2_safe;
                        double sr6 = sr2 * sr2 * sr2;
                        double sr12 = sr6 * sr6;
                        
                        sys->E_pot += 4.0 * sys->epsilon * (sr12 - sr6);
                        
                        double f_mag = 24.0 * sys->epsilon / r2_safe * (2.0 * sr12 - sr6);
                        if (f_mag > MAX_FORCE) f_mag = MAX_FORCE;
                        if (f_mag < -MAX_FORCE) f_mag = -MAX_FORCE;
                        
                        sys->force[i].x += f_mag * rx;
                        sys->force[i].y += f_mag * ry;
                        sys->force[j].x -= f_mag * rx;
                        sys->force[j].y -= f_mag * ry;
                    }
                }
                j = sys->list[j];
            }
        }
    }
}

static void compute_energies(System *sys) {
    double v2_sum = 0.0;
    for (int i = 0; i < sys->N; i++) {
        v2_sum += sys->vel[i].x * sys->vel[i].x + sys->vel[i].y * sys->vel[i].y;
    }
    sys->E_kin = 0.5 * sys->mass * v2_sum;
    sys->E_tot = sys->E_kin + sys->E_pot;
    double dof = 2.0 * (sys->N - 1);
    sys->T_inst = 2.0 * sys->E_kin / (dof * KB);
}

/* VELOCITY RESCALING */
static void rescale_velocities(System *sys) {
    if (sys->T_inst > 1.0 && !isnan(sys->T_inst)) {
        double factor = sqrt(sys->T_target / sys->T_inst);
        if (factor < 0.5) factor = 0.5;
        if (factor > 2.0) factor = 2.0;
        for (int i = 0; i < sys->N; i++) {
            sys->vel[i].x *= factor;
            sys->vel[i].y *= factor;
        }
        compute_energies(sys);
    }
}

/* ==============================================================================
 * FUNZIONI ESPORTATE WEBASSEMBLY
 * ============================================================================== */

EMSCRIPTEN_KEEPALIVE
void free_simulation() {
    if (g_sys != NULL) {
        free(g_sys->pos);
        free(g_sys->vel);
        free(g_sys->force);
        for (int i = 0; i < g_sys->M; i++) {
            free(g_sys->head[i]);
        }
        free(g_sys->head);
        free(g_sys->list);
        free(g_sys);
        g_sys = NULL;
    }
}

EMSCRIPTEN_KEEPALIVE
void init_simulation(int N, double density, double target_temp, double dt, int mode_nvt, double sigma, double epsilon, double mass) {
    free_simulation();
    
    g_sys = (System*)malloc(sizeof(System));
    g_sys->N = N;
    g_sys->density = density;
    g_sys->L = sqrt((double)N / density);
    g_sys->T_target = target_temp;
    g_sys->dt = dt;
    g_sys->mode_nvt = mode_nvt;
    
    /* Configurazione parametri della specie atomica (fallback ad Argon se non specificati) */
    g_sys->sigma = (sigma > 0.1) ? sigma : DEFAULT_SIGMA;
    g_sys->epsilon = (epsilon > 0.1) ? epsilon : DEFAULT_EPSILON;
    g_sys->mass = (mass > 0.1) ? mass : DEFAULT_MASS;
    
    g_sys->pos = (Vec2D*)malloc(N * sizeof(Vec2D));
    g_sys->vel = (Vec2D*)malloc(N * sizeof(Vec2D));
    g_sys->force = (Vec2D*)malloc(N * sizeof(Vec2D));
    
    double r_cutoff = 2.5 * g_sys->sigma;
    if (r_cutoff > g_sys->L / 2.0) r_cutoff = g_sys->L / 2.0;
    
    g_sys->M = (int)(g_sys->L / r_cutoff);
    if (g_sys->M < 3) g_sys->M = 3;
    g_sys->cell_size = g_sys->L / g_sys->M;
    
    g_sys->head = (int**)malloc(g_sys->M * sizeof(int*));
    for (int i = 0; i < g_sys->M; i++) {
        g_sys->head[i] = (int*)malloc(g_sys->M * sizeof(int));
    }
    g_sys->list = (int*)malloc(N * sizeof(int));
    
    g_sys->time = 0.0;
    g_sys->step = 0;
    
    int n_side = (int)ceil(sqrt((double)N));
    double spacing = g_sys->L / n_side;
    int idx = 0;
    for (int i = 0; i < n_side && idx < N; i++) {
        for (int j = 0; j < n_side && idx < N; j++) {
            g_sys->pos[idx].x = (i + 0.5) * spacing;
            g_sys->pos[idx].y = (j + 0.5) * spacing;
            idx++;
        }
    }
    
    double stddev = sqrt(KB * target_temp / g_sys->mass);
    Vec2D v_com = {0.0, 0.0};
    for (int i = 0; i < N; i++) {
        g_sys->vel[i].x = rand_gaussian(0.0, stddev);
        g_sys->vel[i].y = rand_gaussian(0.0, stddev);
        v_com.x += g_sys->vel[i].x;
        v_com.y += g_sys->vel[i].y;
    }
    v_com.x /= N;
    v_com.y /= N;
    for (int i = 0; i < N; i++) {
        g_sys->vel[i].x -= v_com.x;
        g_sys->vel[i].y -= v_com.y;
    }
    
    compute_forces_cell(g_sys);
    compute_energies(g_sys);
    
    /* Pre-equilibrazione fluida con velocity rescaling */
    double dt_m = g_sys->dt / g_sys->mass;
    for (int eq = 0; eq < 200; eq++) {
        for (int i = 0; i < g_sys->N; i++) {
            g_sys->pos[i].x += g_sys->vel[i].x * g_sys->dt + 0.5 * g_sys->force[i].x * dt_m * g_sys->dt;
            g_sys->pos[i].y += g_sys->vel[i].y * g_sys->dt + 0.5 * g_sys->force[i].y * dt_m * g_sys->dt;
        }
        wrap_positions(g_sys);
        
        for (int i = 0; i < g_sys->N; i++) {
            g_sys->vel[i].x += 0.5 * g_sys->force[i].x * dt_m;
            g_sys->vel[i].y += 0.5 * g_sys->force[i].y * dt_m;
        }
        
        compute_forces_cell(g_sys);
        
        for (int i = 0; i < g_sys->N; i++) {
            g_sys->vel[i].x += 0.5 * g_sys->force[i].x * dt_m;
            g_sys->vel[i].y += 0.5 * g_sys->force[i].y * dt_m;
        }
        
        compute_energies(g_sys);
        rescale_velocities(g_sys);
    }
    g_gr_active = 0;
    g_gr_target_steps = 0;
    g_gr_collected_steps = 0;

    compute_energies(g_sys);
}

EMSCRIPTEN_KEEPALIVE
void start_gr_sampling(int target_steps) {
    if (g_sys == NULL) return;
    for (int k = 0; k < GR_NBINS; k++) {
        g_r_histogram[k] = 0.0;
        g_r_result[k] = 0.0;
    }
    g_gr_dr = (g_sys->L * 0.5) / GR_NBINS;
    for (int k = 0; k < GR_NBINS; k++) {
        /* Distanza ridotta r* = r / sigma della specie corrente */
        g_r_axis[k] = ((k + 0.5) * g_gr_dr) / g_sys->sigma;
    }
    g_gr_target_steps = target_steps > 0 ? target_steps : 10000;
    g_gr_collected_steps = 0;
    g_gr_active = 1;
}

EMSCRIPTEN_KEEPALIVE
void step_simulation(int num_steps) {
    if (g_sys == NULL) return;
    
    double dt_m = g_sys->dt / g_sys->mass;
    for (int s = 0; s < num_steps; s++) {
        for (int i = 0; i < g_sys->N; i++) {
            g_sys->pos[i].x += g_sys->vel[i].x * g_sys->dt + 0.5 * g_sys->force[i].x * dt_m * g_sys->dt;
            g_sys->pos[i].y += g_sys->vel[i].y * g_sys->dt + 0.5 * g_sys->force[i].y * dt_m * g_sys->dt;
        }
        
        wrap_positions(g_sys);
        
        for (int i = 0; i < g_sys->N; i++) {
            g_sys->vel[i].x += 0.5 * g_sys->force[i].x * dt_m;
            g_sys->vel[i].y += 0.5 * g_sys->force[i].y * dt_m;
        }
        
        compute_forces_cell(g_sys);
        
        for (int i = 0; i < g_sys->N; i++) {
            g_sys->vel[i].x += 0.5 * g_sys->force[i].x * dt_m;
            g_sys->vel[i].y += 0.5 * g_sys->force[i].y * dt_m;
        }
        
        compute_energies(g_sys);
        if (g_sys->mode_nvt && (g_sys->step % 5 == 0)) {
            rescale_velocities(g_sys);
        }
        
        if (g_gr_active == 1) {
            double r_max = g_sys->L * 0.5;
            double r2_max = r_max * r_max;
            double inv_dr = 1.0 / g_gr_dr;
            
            for (int i = 0; i < g_sys->N; i++) {
                for (int j = i + 1; j < g_sys->N; j++) {
                    double rx = g_sys->pos[i].x - g_sys->pos[j].x;
                    double ry = g_sys->pos[i].y - g_sys->pos[j].y;
                    rx -= g_sys->L * round(rx / g_sys->L);
                    ry -= g_sys->L * round(ry / g_sys->L);
                    double r2 = rx * rx + ry * ry;
                    if (r2 < r2_max) {
                        double r = sqrt(r2);
                        int bin = (int)(r * inv_dr);
                        if (bin >= 0 && bin < GR_NBINS) {
                            g_r_histogram[bin] += 1.0;
                        }
                    }
                }
            }
            g_gr_collected_steps++;

            if (g_gr_collected_steps >= g_gr_target_steps) {
                double rho = (double)g_sys->N / (g_sys->L * g_sys->L);
                double frames = (double)g_gr_collected_steps;
                double dr2 = g_gr_dr * g_gr_dr;
                
                for (int k = 0; k < GR_NBINS; k++) {
                    double shell_area = M_PI * (2.0 * k + 1.0) * dr2;
                    double ideal_count = rho * shell_area;
                    if (ideal_count > 1e-12) {
                        g_r_result[k] = (2.0 * g_r_histogram[k] / (frames * (double)g_sys->N)) / ideal_count;
                    } else {
                        g_r_result[k] = 0.0;
                    }
                }
                g_gr_active = 2;
            }
        }
        
        g_sys->time += g_sys->dt;
        g_sys->step++;
    }
    compute_energies(g_sys);
}

EMSCRIPTEN_KEEPALIVE
int get_gr_status() { return g_gr_active; }

EMSCRIPTEN_KEEPALIVE
int get_gr_progress() { return g_gr_collected_steps; }

EMSCRIPTEN_KEEPALIVE
int get_gr_target() { return g_gr_target_steps; }

EMSCRIPTEN_KEEPALIVE
int get_gr_nbins() { return GR_NBINS; }

EMSCRIPTEN_KEEPALIVE
double* get_gr_r_ptr() { return g_r_axis; }

EMSCRIPTEN_KEEPALIVE
double* get_gr_result_ptr() { return g_r_result; }

EMSCRIPTEN_KEEPALIVE
double* get_positions_ptr() {
    return (double*)g_sys->pos;
}

EMSCRIPTEN_KEEPALIVE
double get_box_length() {
    return g_sys ? g_sys->L : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_reduced_box_length() {
    return g_sys ? (g_sys->L / g_sys->sigma) : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_temperature() {
    return g_sys ? g_sys->T_inst : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_reduced_temperature() {
    return g_sys ? (g_sys->T_inst * KB / g_sys->epsilon) : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_energy_kin() {
    return g_sys ? g_sys->E_kin : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_energy_pot() {
    return g_sys ? g_sys->E_pot : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_energy_tot() {
    return g_sys ? g_sys->E_tot : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_reduced_energy_tot_per_particle() {
    return (g_sys && g_sys->N > 0) ? (g_sys->E_tot / (g_sys->N * g_sys->epsilon)) : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_reduced_density() {
    return g_sys ? (g_sys->density * g_sys->sigma * g_sys->sigma) : 0.0;
}

EMSCRIPTEN_KEEPALIVE
int get_step() {
    return g_sys ? g_sys->step : 0;
}

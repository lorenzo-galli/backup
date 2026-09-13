#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <emscripten/emscripten.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Costanti fisiche Argon (Angstrom, ps, u, K) */
#define EPSILON 99.365319
#define SIGMA 3.4
#define MASS 39.9
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
    
    double r2_cut = (2.5 * SIGMA) * (2.5 * SIGMA);
    if (r2_cut > (sys->L * 0.5) * (sys->L * 0.5)) {
        r2_cut = (sys->L * 0.5) * (sys->L * 0.5);
    }
    
    double r2_min = (0.65 * SIGMA) * (0.65 * SIGMA);
    double sig2 = SIGMA * SIGMA;
    
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
                        
                        sys->E_pot += 4.0 * EPSILON * (sr12 - sr6);
                        
                        double f_mag = 24.0 * EPSILON / r2_safe * (2.0 * sr12 - sr6);
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
    sys->E_kin = 0.5 * MASS * v2_sum;
    sys->E_tot = sys->E_kin + sys->E_pot;
    double dof = 2.0 * (sys->N - 1);
    sys->T_inst = 2.0 * sys->E_kin / (dof * KB);
}

static void rescale_velocities(System *sys) {
    if (sys->T_inst > 0.0 && !isnan(sys->T_inst)) {
        double factor = sqrt(sys->T_target / sys->T_inst);
        if (factor < 0.2) factor = 0.2;
        if (factor > 5.0) factor = 5.0;
        for (int i = 0; i < sys->N; i++) {
            sys->vel[i].x *= factor;
            sys->vel[i].y *= factor;
        }
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
void init_simulation(int N, double density, double target_temp, double dt, int mode_nvt) {
    free_simulation();
    
    g_sys = (System*)malloc(sizeof(System));
    g_sys->N = N;
    g_sys->density = density;
    g_sys->L = sqrt((double)N / density);
    g_sys->T_target = target_temp;
    g_sys->dt = dt;
    g_sys->mode_nvt = mode_nvt;
    
    g_sys->pos = (Vec2D*)malloc(N * sizeof(Vec2D));
    g_sys->vel = (Vec2D*)malloc(N * sizeof(Vec2D));
    g_sys->force = (Vec2D*)malloc(N * sizeof(Vec2D));
    
    double r_cutoff = 2.5 * SIGMA;
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
    
    double stddev = sqrt(KB * target_temp / MASS);
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
    
    /* Pre-equilibrazione ultra-veloce di 150 step */
    double dt_m = g_sys->dt / MASS;
    for (int eq = 0; eq < 150; eq++) {
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
        
        if (eq % 5 == 0) {
            compute_energies(g_sys);
            rescale_velocities(g_sys);
        }
    }
    compute_energies(g_sys);
}

EMSCRIPTEN_KEEPALIVE
void step_simulation(int num_steps) {
    if (g_sys == NULL) return;
    
    double dt_m = g_sys->dt / MASS;
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
        
        if (g_sys->mode_nvt && (g_sys->step % 5 == 0)) {
            rescale_velocities(g_sys);
        }
        
        g_sys->time += g_sys->dt;
        g_sys->step++;
    }
    compute_energies(g_sys);
}

EMSCRIPTEN_KEEPALIVE
double* get_positions_ptr() {
    return (double*)g_sys->pos;
}

EMSCRIPTEN_KEEPALIVE
double get_box_length() {
    return g_sys ? g_sys->L : 0.0;
}

EMSCRIPTEN_KEEPALIVE
double get_temperature() {
    return g_sys ? g_sys->T_inst : 0.0;
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
int get_step() {
    return g_sys ? g_sys->step : 0;
}

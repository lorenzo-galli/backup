#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================================
 * CONFIGURAZIONE E COSTANTI FISICHE (Unità: Angstrom, ps, u, K)
 * ============================================================================== */
#define DEFAULT_DENSITY 0.1
#define TARGET_TEMP 150.0
#define N_ATOMS 1000
#define DEFAULT_EQ_STEPS 5000
#define DEFAULT_PROD_STEPS 10000
#define DT 0.01

#define EPSILON 99.365319   /* Parametro L-J in unità E_0 */
#define SIGMA 3.4           /* Diametro L-J in Angstrom */
#define MASS 39.9           /* Massa Argon in u */
#define KB 0.831446         /* Costante di Boltzmann (E_0/K) */
#define MAX_FORCE 5000.0    /* Force capping anti-NaN per alte densità */

/* Struttura per vettori 2D */
typedef struct {
    double x;
    double y;
} Vec2D;

/* Struttura del Sistema di Dinamica Molecolare */
typedef struct {
    int N;
    double density;
    double L;
    
    Vec2D *pos;
    Vec2D *vel;
    Vec2D *force;
    
    /* Griglia per il Metodo delle Celle (Linked-Cell List) */
    int M;             /* Numero di celle per lato (M x M) */
    double cell_size;  /* Dimensione cella in Angstrom */
    int **head;        /* Matrice M x M con indice prima particella */
    int *list;         /* Array N collegato alla particella successiva */
    
    double time;
    int step;
    
    double E_kin;
    double E_pot;
    double E_tot;
    double T_inst;
} System;

/* Algoritmo di Box-Muller per genera numeri casuali gaussiani */
double rand_gaussian(double mean, double stddev) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    if (u1 < 1e-10) u1 = 1e-10;
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return mean + z0 * stddev;
}

/* Inizializzazione del Sistema */
System* create_system(int N, double density) {
    System *sys = (System*)malloc(sizeof(System));
    sys->N = N;
    sys->density = density;
    sys->L = sqrt((double)N / density);
    
    sys->pos = (Vec2D*)malloc(N * sizeof(Vec2D));
    sys->vel = (Vec2D*)malloc(N * sizeof(Vec2D));
    sys->force = (Vec2D*)malloc(N * sizeof(Vec2D));
    
    /* Configurazione Griglia per Metodo delle Celle */
    double r_cutoff = 2.5 * SIGMA;
    if (r_cutoff > sys->L / 2.0) {
        r_cutoff = sys->L / 2.0;
    }
    
    sys->M = (int)(sys->L / r_cutoff);
    if (sys->M < 3) sys->M = 3; /* Almeno 3x3 celle per periodicità 2D */
    sys->cell_size = sys->L / sys->M;
    
    /* Allocazione griglia celle */
    sys->head = (int**)malloc(sys->M * sizeof(int*));
    for (int i = 0; i < sys->M; i++) {
        sys->head[i] = (int*)malloc(sys->M * sizeof(int));
    }
    sys->list = (int*)malloc(N * sizeof(int));
    
    sys->time = 0.0;
    sys->step = 0;
    
    /* Inizializzazione Reticolo Quadrato */
    int n_side = (int)ceil(sqrt((double)N));
    double spacing = sys->L / n_side;
    int idx = 0;
    for (int i = 0; i < n_side && idx < N; i++) {
        for (int j = 0; j < n_side && idx < N; j++) {
            sys->pos[idx].x = (i + 0.5) * spacing;
            sys->pos[idx].y = (j + 0.5) * spacing;
            idx++;
        }
    }
    
    /* Inizializzazione Velocità Maxwell-Boltzmann */
    double stddev = sqrt(KB * TARGET_TEMP / MASS);
    Vec2D v_com = {0.0, 0.0};
    for (int i = 0; i < N; i++) {
        sys->vel[i].x = rand_gaussian(0.0, stddev);
        sys->vel[i].y = rand_gaussian(0.0, stddev);
        v_com.x += sys->vel[i].x;
        v_com.y += sys->vel[i].y;
    }
    
    /* Rimozione del moto del centro di massa (V_COM = 0) */
    v_com.x /= N;
    v_com.y /= N;
    for (int i = 0; i < N; i++) {
        sys->vel[i].x -= v_com.x;
        sys->vel[i].y -= v_com.y;
    }
    
    return sys;
}

/* Aggiornamento della struttura Linked-Cell List */
void build_cell_list(System *sys) {
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
        
        if (cx >= sys->M) cx = sys->M - 1;
        if (cy >= sys->M) cy = sys->M - 1;
        
        sys->list[i] = sys->head[cx][cy];
        sys->head[cx][cy] = i;
    }
}

/* Calcolo Forze Lennard-Jones e Potenziale tramite Metodo delle Celle */
void compute_forces_cell_method(System *sys) {
    for (int i = 0; i < sys->N; i++) {
        sys->force[i].x = 0.0;
        sys->force[i].y = 0.0;
    }
    sys->E_pot = 0.0;
    
    build_cell_list(sys);
    
    double r2_min = (0.7 * SIGMA) * (0.7 * SIGMA);
    double sig2 = SIGMA * SIGMA;
    
    for (int cx = 0; cx < sys->M; cx++) {
        for (int cy = 0; cy < sys->M; cy++) {
            
            for (int dx = -1; dx <= 1; dx++) {
                for (int dy = -1; dy <= 1; dy++) {
                    
                    int ncx = (cx + dx + sys->M) % sys->M;
                    int ncy = (cy + dy + sys->M) % sys->M;
                    
                    int i = sys->head[cx][cy];
                    while (i != -1) {
                        int j = sys->head[ncx][ncy];
                        while (j != -1) {
                            
                            /* Evita auto-interazione e doppio conteggio */
                            if (cx == ncx && cy == ncy && j <= i) {
                                j = sys->list[j];
                                continue;
                            }
                            
                            if (ncx < cx || (ncx == cx && ncy < cy)) {
                                j = sys->list[j];
                                continue;
                            }
                            
                            /* Vettore r_i - r_j */
                            double rx = sys->pos[i].x - sys->pos[j].x;
                            double ry = sys->pos[i].y - sys->pos[j].y;
                            
                            /* Minimum Image Convention */
                            rx -= sys->L * round(rx / sys->L);
                            ry -= sys->L * round(ry / sys->L);
                            
                            double r2 = rx * rx + ry * ry;
                            double r2_safe = (r2 < r2_min) ? r2_min : r2;
                            
                            double sr2 = sig2 / r2_safe;
                            double sr6 = sr2 * sr2 * sr2;
                            double sr12 = sr6 * sr6;
                            
                            /* Potenziale V_ij */
                            sys->E_pot += 4.0 * EPSILON * (sr12 - sr6);
                            
                            /* Magnitudine Forza */
                            double f_mag = 24.0 * EPSILON / r2_safe * (2.0 * sr12 - sr6);
                            
                            /* Force capping anti-NaN */
                            if (f_mag > MAX_FORCE) f_mag = MAX_FORCE;
                            if (f_mag < -MAX_FORCE) f_mag = -MAX_FORCE;
                            
                            sys->force[i].x += f_mag * rx;
                            sys->force[i].y += f_mag * ry;
                            sys->force[j].x -= f_mag * rx;
                            sys->force[j].y -= f_mag * ry;
                            
                            j = sys->list[j];
                        }
                        i = sys->list[i];
                    }
                }
            }
        }
    }
}

/* Calcolo Energie Macroscopiche e Temperatura */
void compute_energies(System *sys) {
    double v2_sum = 0.0;
    for (int i = 0; i < sys->N; i++) {
        v2_sum += sys->vel[i].x * sys->vel[i].x + sys->vel[i].y * sys->vel[i].y;
    }
    sys->E_kin = 0.5 * MASS * v2_sum;
    sys->E_tot = sys->E_kin + sys->E_pot;
    
    double dof = 2.0 * (sys->N - 1);
    sys->T_inst = 2.0 * sys->E_kin / (dof * KB);
}

/* Velocity Rescaling per NVT */
void rescale_velocities(System *sys) {
    if (sys->T_inst > 0.0 && !isnan(sys->T_inst)) {
        double factor = sqrt(TARGET_TEMP / sys->T_inst);
        if (factor < 0.5) factor = 0.5;
        if (factor > 2.0) factor = 2.0;
        
        for (int i = 0; i < sys->N; i++) {
            sys->vel[i].x *= factor;
            sys->vel[i].y *= factor;
        }
    }
}

/* Integratore Velocity Verlet */
void step_verlet(System *sys) {
    double dt_m = DT / MASS;
    
    /* 1. Aggiorna posizioni */
    for (int i = 0; i < sys->N; i++) {
        sys->pos[i].x += sys->vel[i].x * DT + 0.5 * sys->force[i].x * dt_m * DT;
        sys->pos[i].y += sys->vel[i].y * DT + 0.5 * sys->force[i].y * dt_m * DT;
    }
    
    /* 2. Mezzo step velocita' */
    for (int i = 0; i < sys->N; i++) {
        sys->vel[i].x += 0.5 * sys->force[i].x * dt_m;
        sys->vel[i].y += 0.5 * sys->force[i].y * dt_m;
    }
    
    /* 3. Calcolo Forze con Metodo delle Celle */
    compute_forces_cell_method(sys);
    
    /* 4. Secondo mezzo step velocita' */
    for (int i = 0; i < sys->N; i++) {
        sys->vel[i].x += 0.5 * sys->force[i].x * dt_m;
        sys->vel[i].y += 0.5 * sys->force[i].y * dt_m;
    }
    
    sys->time += DT;
    sys->step++;
    compute_energies(sys);
}

void free_system(System *sys) {
    free(sys->pos);
    free(sys->vel);
    free(sys->force);
    for (int i = 0; i < sys->M; i++) {
        free(sys->head[i]);
    }
    free(sys->head);
    free(sys->list);
    free(sys);
}

/* Main */
int main(int argc, char *argv[]) {
    srand(42); /* Seed fisso per riproducibilità */
    
    double density = DEFAULT_DENSITY;
    int eq_steps = DEFAULT_EQ_STEPS;
    int prod_steps = DEFAULT_PROD_STEPS;
    
    if (argc > 1) {
        density = atof(argv[1]);
    }
    if (argc > 2) {
        eq_steps = atoi(argv[2]);
    }
    if (argc > 3) {
        prod_steps = atoi(argv[3]);
    }
    
    char prefix[256];
    snprintf(prefix, sizeof(prefix), "run_rho_%g", density);
    
    printf("==================================================\n");
    printf(" Simulazione MD in C (Metodo delle Celle) - rho = %g\n", density);
    printf("==================================================\n");
    
    System *sys = create_system(N_ATOMS, density);
    compute_forces_cell_method(sys);
    compute_energies(sys);
    
    char file_energies[256], file_traj[256], file_vel[256];
    snprintf(file_energies, sizeof(file_energies), "%s_energies.csv", prefix);
    snprintf(file_traj, sizeof(file_traj), "%s_traj.xyz", prefix);
    snprintf(file_vel, sizeof(file_vel), "%s_vel.xyz", prefix);
    
    FILE *f_eng = fopen(file_energies, "w");
    FILE *f_trj = fopen(file_traj, "w");
    FILE *f_vel = fopen(file_vel, "w");
    
    if (!f_eng || !f_trj || !f_vel) {
        fprintf(stderr, "Errore nell'apertura dei file di output!\n");
        return 1;
    }
    
    fprintf(f_eng, "step,time,E_kin,E_pot,E_tot,T\n");
    
    printf("  --> Equilibrazione NVT (%d step)...\n", eq_steps);
    for (int i = 0; i < eq_steps; i++) {
        step_verlet(sys);
        if (i % 10 == 0) {
            rescale_velocities(sys);
        }
        
        fprintf(f_eng, "%d,%.6f,%.4f,%.4f,%.4f,%.4f\n", sys->step, sys->time, sys->E_kin, sys->E_pot, sys->E_tot, sys->T_inst);
        
        if (i % 1000 == 0) {
            printf("    [Step %d/%d] T = %.2f K, E_tot = %.2f\n", i, eq_steps, sys->T_inst, sys->E_tot);
        }
    }
    
    printf("  --> Produzione NVE (%d step)...\n", prod_steps);
    for (int i = 0; i < prod_steps; i++) {
        step_verlet(sys);
        
        fprintf(f_eng, "%d,%.6f,%.4f,%.4f,%.4f,%.4f\n", sys->step, sys->time, sys->E_kin, sys->E_pot, sys->E_tot, sys->T_inst);
        
        if (i % 10 == 0) {
            fprintf(f_trj, "%d\nFrame %d\n", sys->N, sys->step);
            fprintf(f_vel, "%d\nFrame %d\n", sys->N, sys->step);
            for (int k = 0; k < sys->N; k++) {
                fprintf(f_trj, "Ar %.6f %.6f 0.0\n", sys->pos[k].x, sys->pos[k].y);
                fprintf(f_vel, "Ar %.6f %.6f 0.0\n", sys->vel[k].x, sys->vel[k].y);
            }
        }
        
        if (i % 1000 == 0) {
            printf("    [Step %d/%d] T = %.2f K, E_tot = %.2f\n", sys->step, sys->step, sys->T_inst, sys->E_tot);
        }
    }
    
    fclose(f_eng);
    fclose(f_trj);
    fclose(f_vel);
    
    printf("Simulazione C completata con successo! Output salvati in '%s'.\n\n", prefix);
    
    free_system(sys);
    return 0;
}

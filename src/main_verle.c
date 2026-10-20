/*
 * main_verle.c
 * Молекулярная динамика: потенциал Леннарда-Джонса, алгоритм Верле
 * Вариант 1: без wrap в method_Verle, только min_image для расстояний
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mt19937-64.c"

// ========== ГЛОБАЛЬНЫЕ НАСТРОЙКИ ==========
#define h 0.0005
#define STEPS 10000
#define SAVE_EVERY 100

// ========== НАСТРОЙКИ СКОРОСТЕЙ ==========
#define GAS_SPEED 0.15
#define LIQUID_SPEED 0.06
#define SOLID_SPEED 0.0008

// ========== ФИЗИЧЕСКИЕ КОНСТАНТЫ ==========
#define sigma 1.0
#define eps 1.0
#define m 1.0
#define rcut (3.5 * sigma)
#define rcut2 (rcut * rcut)

#define MAX_N 400

double rx[MAX_N], ry[MAX_N];
double rx_before[MAX_N], ry_before[MAX_N];
double rx_next[MAX_N], ry_next[MAX_N];
double vx[MAX_N], vy[MAX_N];
double ax[MAX_N], ay[MAX_N];

int N;
double L;

double get_random(double min_value, double max_value)
{
    return genrand64_real2() * (max_value - min_value) + min_value;
}

// Периодические границы для позиций (НЕ ИСПОЛЬЗУЕТСЯ В method_Verle, только для начальной инициализации)
double wrap(double x)
{
    if (x >= L) x -= L;
    if (x < 0.0) x += L;
    return x;
}

// Минимальное изображение для расстояний (единственное место, где учитываются периодические границы)
double min_image(double dx)
{
    if (dx > L / 2.0) dx -= L;
    if (dx < -L / 2.0) dx += L;
    return dx;
}

// ========== ИНИЦИАЛИЗАЦИЯ ПОЗИЦИЙ ==========

void set_random_positions()
{
    for (int i = 0; i < N; i++) {
        rx[i] = get_random(0, L);
        ry[i] = get_random(0, L);
    }
}

void set_square_positions()
{
    int rows = (int)sqrt(N);
    if (rows * rows < N) rows++;
    
    double cell_size = L / rows;
    
    int k = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            if (k >= N) return;
            rx[k] = (j + 0.5) * cell_size;
            ry[k] = (i + 0.5) * cell_size;
            k++;
        }
    }
}

void set_triangular_positions()
{
    int rows = (int)sqrt(N);
    if (rows * rows < N) rows++;
    
    double a = L / rows;
    if (a < sigma) a = sigma;
    
    int k = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            if (k >= N) return;
            rx[k] = j * a;
            ry[k] = i * a * sqrt(3) / 2.0;
            if (i % 2 == 1) rx[k] += a / 2.0;
            k++;
        }
    }
}

// ========== ИНИЦИАЛИЗАЦИЯ СКОРОСТЕЙ ==========

void set_random_velocities(double max_v)
{
    for (int i = 0; i < N; i++) {
        vx[i] = get_random(-max_v, max_v);
        vy[i] = get_random(-max_v, max_v);
    }
}

void set_small_velocities()
{
    for (int i = 0; i < N; i++) {
        vx[i] = get_random(-SOLID_SPEED, SOLID_SPEED);
        vy[i] = get_random(-SOLID_SPEED, SOLID_SPEED);
    }
}

// ========== РАСЧЕТ СИЛ ==========

void particle_forces()
{
    const double sigma6 = pow(sigma, 6);
    const double sigma12 = pow(sigma, 12);
    
    for (int i = 0; i < N; i++) {
        ax[i] = 0.0;
        ay[i] = 0.0;
        
        for (int j = 0; j < N; j++) {
            if (j == i) continue;
            
            double dx = min_image(rx[j] - rx[i]);
            double dy = min_image(ry[j] - ry[i]);
            double r2 = dx*dx + dy*dy;
            
            if (r2 > rcut2) continue;
            if (r2 < 0.01) continue;
            
            double r2_inv = 1.0 / r2;
            double r6_inv = r2_inv * r2_inv * r2_inv;
            double r12_inv = r6_inv * r6_inv;
            
            double f = 24.0 * eps * (2.0 * sigma12 * r12_inv - sigma6 * r6_inv) * r2_inv;
            
            ax[i] += f * dx;
            ay[i] += f * dy;
        }
    }
    
    for (int i = 0; i < N; i++) {
        ax[i] /= m;
        ay[i] /= m;
    }
}

// ========== РАСЧЕТ ЭНЕРГИИ ==========

void calculate_energies(double *Ek, double *Ep)
{
    const double sigma6 = pow(sigma, 6);
    const double sigma12 = pow(sigma, 12);
    
    *Ek = 0.0;
    *Ep = 0.0;
    
    for (int i = 0; i < N; i++) {
        *Ek += 0.5 * m * (vx[i]*vx[i] + vy[i]*vy[i]);
    }
    
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = min_image(rx[j] - rx[i]);
            double dy = min_image(ry[j] - ry[i]);
            double r2 = dx*dx + dy*dy;
            
            if (r2 > rcut2) continue;
            if (r2 < 0.01) continue;
            
            double r2_inv = 1.0 / r2;
            double r6_inv = r2_inv * r2_inv * r2_inv;
            double r12_inv = r6_inv * r6_inv;
            *Ep += 4 * eps * (r12_inv - r6_inv);
        }
    }
}

// ========== АЛГОРИТМ ВЕРЛЕ (без wrap) ==========

void method_Verle()
{
    particle_forces();
    
    for (int i = 0; i < N; i++) {
        rx_next[i] = 2*rx[i] - rx_before[i] + ax[i]*h*h;
        ry_next[i] = 2*ry[i] - ry_before[i] + ay[i]*h*h;
        
        vx[i] = (rx_next[i] - rx_before[i]) / (2*h);
        vy[i] = (ry_next[i] - ry_before[i]) / (2*h);
    }
    
    for (int i = 0; i < N; i++) {
        rx_before[i] = rx[i];
        ry_before[i] = ry[i];
        rx[i] = rx_next[i];
        ry[i] = ry_next[i];
    }
}

// ========== ИНИЦИАЛИЗАЦИЯ ==========

void init(int pos_type, int vel_type, double vel_max)
{
    if (pos_type == 0) {
        set_random_positions();
    } else if (pos_type == 1) {
        set_square_positions();
    } else if (pos_type == 2) {
        set_triangular_positions();
    }
    
    if (vel_type == 0) {
        set_random_velocities(vel_max);
    } else {
        set_small_velocities();
    }
    
    particle_forces();
    
    for (int i = 0; i < N; i++) {
        rx_before[i] = rx[i];
        ry_before[i] = ry[i];
        rx[i] = rx_before[i] + vx[i]*h + 0.5*ax[i]*h*h;
        ry[i] = ry_before[i] + vy[i]*h + 0.5*ay[i]*h*h;
    }
}

// ========== ВЫВОД В ФАЙЛЫ ==========

void write_to_file(int step, const char* folder)
{
    char fname[200];
    sprintf(fname, "%s/conf_%d.csv", folder, step);
    FILE *conffile = fopen(fname, "w");
    if (!conffile) return;
    fprintf(conffile, "x,y,vx,vy\n");
    for (int i = 0; i < N; i++)
        fprintf(conffile, "%f,%f,%f,%f\n", rx[i], ry[i], vx[i], vy[i]);
    fclose(conffile);
}

void write_energy_to_file(int step, double Ek, double Ep, const char* folder)
{
    static int first = 1;
    static FILE *efile = NULL;
    if (first) {
        char fname[200];
        sprintf(fname, "%s/energy.csv", folder);
        efile = fopen(fname, "w");
        if (efile) fprintf(efile, "step,total,kinetic,potential\n");
        first = 0;
    }
    if (efile) {
        fprintf(efile, "%d,%f,%f,%f\n", step, Ek+Ep, Ek, Ep);
        fflush(efile);
    }
}

// ========== РЕЖИМЫ СИМУЛЯЦИИ ==========

#ifdef STATE_GAS
int main()
{
    init_genrand64(time(NULL));
    N = 35; L = 30.0;
    const char* output_folder = "data/gas";
    char cmd[100];
    sprintf(cmd, "mkdir -p %s", output_folder);
    system(cmd);
    init(0, 0, GAS_SPEED);
    printf("\n=== ГАЗ (случайные положения) ===\nN = %d, L = %.1f\n\n", N, L);
    for (int step = 0; step <= STEPS; step++) {
        if (step > 0) method_Verle();
        if (step % SAVE_EVERY == 0) {
            write_to_file(step, output_folder);
            double Ek, Ep; calculate_energies(&Ek, &Ep);
            write_energy_to_file(step, Ek, Ep, output_folder);
            printf("Шаг %d: Ek = %8.4f, Ep = %8.4f, E = %8.4f\n", step, Ek, Ep, Ek+Ep);
        }
    }
    return 0;
}
#endif

#ifdef STATE_GAS_REGULAR
int main()
{
    init_genrand64(time(NULL));
    N = 10; L = 12.0;
    const char* output_folder = "data/gas_regular";
    char cmd[100];
    sprintf(cmd, "mkdir -p %s", output_folder);
    system(cmd);
    init(1, 0, GAS_SPEED);
    printf("\n=== ГАЗ (регулярная решетка) ===\nN = %d, L = %.1f\n\n", N, L);
    for (int step = 0; step <= STEPS; step++) {
        if (step > 0) method_Verle();
        if (step % SAVE_EVERY == 0) {
            write_to_file(step, output_folder);
            double Ek, Ep; calculate_energies(&Ek, &Ep);
            write_energy_to_file(step, Ek, Ep, output_folder);
            printf("Шаг %d: Ek = %8.4f, Ep = %8.4f, E = %8.4f\n", step, Ek, Ep, Ek+Ep);
        }
    }
    return 0;
}
#endif

#ifdef STATE_LIQUID
int main()
{
    init_genrand64(time(NULL));
    N = 36; L = 12.0;
    const char* output_folder = "data/liquid";
    char cmd[100];
    sprintf(cmd, "mkdir -p %s", output_folder);
    system(cmd);
    init(1, 0, LIQUID_SPEED);
    printf("\n=== ЖИДКОСТЬ (малая) ===\nN = %d, L = %.1f\n\n", N, L);
    for (int step = 0; step <= STEPS; step++) {
        if (step > 0) method_Verle();
        if (step % SAVE_EVERY == 0) {
            write_to_file(step, output_folder);
            double Ek, Ep; calculate_energies(&Ek, &Ep);
            write_energy_to_file(step, Ek, Ep, output_folder);
            printf("Шаг %d: Ek = %8.4f, Ep = %8.4f, E = %8.4f\n", step, Ek, Ep, Ek+Ep);
        }
    }
    return 0;
}
#endif

#ifdef STATE_LIQUID_LARGE
int main()
{
    init_genrand64(time(NULL));
    N = 100; L = 20.0;
    const char* output_folder = "data/liquid_large";
    char cmd[100];
    sprintf(cmd, "mkdir -p %s", output_folder);
    system(cmd);
    init(1, 0, LIQUID_SPEED);
    printf("\n=== ЖИДКОСТЬ (большая) ===\nN = %d, L = %.1f\n\n", N, L);
    for (int step = 0; step <= STEPS; step++) {
        if (step > 0) method_Verle();
        if (step % SAVE_EVERY == 0) {
            write_to_file(step, output_folder);
            double Ek, Ep; calculate_energies(&Ek, &Ep);
            write_energy_to_file(step, Ek, Ep, output_folder);
            printf("Шаг %d: Ek = %8.4f, Ep = %8.4f, E = %8.4f\n", step, Ek, Ep, Ek+Ep);
        }
    }
    return 0;
}
#endif

#ifdef STATE_SOLID
int main()
{
    init_genrand64(time(NULL));
    N = 25; L = 12.0;
    const char* output_folder = "data/solid";
    char cmd[100];
    sprintf(cmd, "mkdir -p %s", output_folder);
    system(cmd);
    init(1, 1, 0.0);
    printf("\n=== ТВЕРДОЕ ТЕЛО (квадратная решетка) ===\nN = %d, L = %.1f\n\n", N, L);
    for (int step = 0; step <= STEPS; step++) {
        if (step > 0) method_Verle();
        if (step % SAVE_EVERY == 0) {
            write_to_file(step, output_folder);
            double Ek, Ep; calculate_energies(&Ek, &Ep);
            write_energy_to_file(step, Ek, Ep, output_folder);
            printf("Шаг %d: Ek = %8.4f, Ep = %8.4f, E = %8.4f\n", step, Ek, Ep, Ek+Ep);
        }
    }
    return 0;
}
#endif

#ifdef STATE_SOLID_TRI
int main()
{
    init_genrand64(time(NULL));
    N = 25; L = 12.0;
    const char* output_folder = "data/solid_tri";
    char cmd[100];
    sprintf(cmd, "mkdir -p %s", output_folder);
    system(cmd);
    init(2, 1, 0.0);
    printf("\n=== ТВЕРДОЕ ТЕЛО (треугольная решетка) ===\nN = %d, L = %.1f\n\n", N, L);
    for (int step = 0; step <= STEPS; step++) {
        if (step > 0) method_Verle();
        if (step % SAVE_EVERY == 0) {
            write_to_file(step, output_folder);
            double Ek, Ep; calculate_energies(&Ek, &Ep);
            write_energy_to_file(step, Ek, Ep, output_folder);
            printf("Шаг %d: Ek = %8.4f, Ep = %8.4f, E = %8.4f\n", step, Ek, Ep, Ek+Ep);
        }
    }
    return 0;
}
#endif

#ifndef STATE_GAS
#ifndef STATE_GAS_REGULAR
#ifndef STATE_LIQUID
#ifndef STATE_LIQUID_LARGE
#ifndef STATE_SOLID
#ifndef STATE_SOLID_TRI
int main()
{
    printf("Используйте флаги: -DSTATE_GAS, -DSTATE_GAS_REGULAR, -DSTATE_LIQUID, -DSTATE_LIQUID_LARGE, -DSTATE_SOLID, -DSTATE_SOLID_TRI\n");
    return 1;
}
#endif
#endif
#endif
#endif
#endif
#endif

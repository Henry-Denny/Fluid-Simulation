#ifndef FLUID_H
#define FLUID_H

#include <stdlib.h>
#include <math.h>

#define ITER 4
#define IX(x, y) ((x) + (y) * N)

typedef struct FluidSquare {
    int size;
    float dt;
    float diff;
    float visc;
    
    float *s;
    float *density;
    
    float *Vx;
    float *Vy;

    float *Vx0;
    float *Vy0;
} FluidSquare;

FluidSquare *FluidSquareCreate(int size, float diffusion, float viscosity, float dt);
void FluidSquareFree(FluidSquare *square);
static void set_bnd(int b, float *x, int N);
static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N);
static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N);
static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float dt, int N);
static void project(float *velocX, float *velocY, float *p, float *div, int iter, int N);
void FluidSquareStep(FluidSquare *square);
void FluidSquareAddDensity(FluidSquare *square, int x, int y, float amount);
void FluidSquareAddVelocity(FluidSquare *square, int x, int y, float amountX, float amountY);

#endif
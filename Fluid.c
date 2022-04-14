#include "Fluid.h"

FluidSquare *FluidSquareCreate(int size, float diffusion, float viscosity, float dt)
{
    FluidSquare *square = malloc(sizeof(*square));
    int N = size;
    
    square->size = size;
    square->dt = dt;
    square->diff = diffusion;
    square->visc = viscosity;
    
    square->s = calloc(N * N, sizeof(float));
    square->density = calloc(N * N, sizeof(float));
    
    square->Vx = calloc(N * N, sizeof(float));
    square->Vy = calloc(N * N, sizeof(float));
    
    square->Vx0 = calloc(N * N, sizeof(float));
    square->Vy0 = calloc(N * N, sizeof(float));
    
    return square;
}

void FluidSquareFree(FluidSquare *square)
{
    free(square->s);
    free(square->density);
    
    free(square->Vx);
    free(square->Vy);
    
    free(square->Vx0);
    free(square->Vy0);
    
    free(square);
}


static void set_bnd(int b, float *x, int N)
{
    for(int i = 1; i < N - 1; i++)
    {
        x[IX(i, 0   )] = b == 2 ? -x[IX(i, 1    )] : x[IX(i, 1      )];
        x[IX(i, N-1 )] = b == 2 ? -x[IX(i, N-2  )] : x[IX(i, N-2    )];
    }
    for(int j = 1; j < N - 1; j++)
    {
        x[IX(0  , j )] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
        x[IX(N-1, j )] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
    }
    
    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N-1)] = 0.5f * (x[IX(1, N-1)] + x[IX(0, N-2)]);
    x[IX(N-1, 0)] = 0.5f * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
    x[IX(N-1, N-1)] = 0.5f * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);
}

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N)
{
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) 
    {
        for (int j = 1; j < N - 1; j++)
        {
            for (int i = 1; i < N - 1; i++)
            {
                x[IX(i, j)] =
                    (x0[IX(i, j)]
                        + a*(    x[IX(i+1, j)]
                                +x[IX(i-1, j)]
                                +x[IX(i  , j+1)]
                                +x[IX(i  , j-1)]
                        )) * cRecip;
            }
        }
        set_bnd(b, x, N);
    }
}

static void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N)
{
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 4 * a, iter, N);
}

static void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float dt, int N)
{
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
    
    for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
        for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x    = ifloat - tmp1; 
            y    = jfloat - tmp2;
            
            if(x < 0.5f) x = 0.5f; 
            if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
            i0 = floorf(x); 
            i1 = i0 + 1.0f;
            if(y < 0.5f) y = 0.5f; 
            if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
            j0 = floorf(y);
            j1 = j0 + 1.0f; 
            
            s1 = x - i0; 
            s0 = 1.0f - s1; 
            t1 = y - j0; 
            t0 = 1.0f - t1;
            
            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;
            
            d[IX(i, j)] = 
                s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
        }
    }
    set_bnd(b, d, N);
}

static void project(float *velocX, float *velocY, float *p, float *div, int iter, int N)
{
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f*(
                    velocX[IX(i+1,  j)]
                    -velocX[IX(i-1, j)]
                    +velocY[IX(i,   j+1)]
                    -velocY[IX(i,   j-1)]
                )/N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1, 4, iter, N);
    
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (  p[IX(i+1, j)]
                                            -p[IX(i-1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (  p[IX(i, j+1)]
                                            -p[IX(i, j-1)]) * N;
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
}

void FluidSquareStep(FluidSquare *square)
{
    int N          = square->size;
    float visc     = square->visc;
    float diff     = square->diff;
    float dt       = square->dt;
    float *Vx      = square->Vx;
    float *Vy      = square->Vy;
    float *Vx0     = square->Vx0;
    float *Vy0     = square->Vy0;
    float *s       = square->s;
    float *density = square->density;
    
    diffuse(1, Vx0, Vx, visc, dt, ITER, N);
    diffuse(2, Vy0, Vy, visc, dt, ITER, N);
    
    project(Vx0, Vy0, Vx, Vy, ITER, N);
    
    advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N);
    
    project(Vx, Vy, Vx0, Vy0, ITER, N);
    
    diffuse(0, s, density, diff, dt, ITER, N);
    advect(0, density, s, Vx, Vy, dt, N);
}

void FluidSquareAddDensity(FluidSquare *square, int x, int y, float amount)
{
    int N = square->size;
    square->density[IX(x, y)] += amount;
}

void FluidSquareAddVelocity(FluidSquare *square, int x, int y, float amountX, float amountY)
{
    int N = square->size;
    int index = IX(x, y);
    
    square->Vx[index] += amountX;
    square->Vy[index] += amountY;
}
/*******************************
 type definition of configuration
********************************/

#ifndef CONFIG_H
#define CONFIG_H

#include"sys.h"

typedef struct{
    double x, vx, fx;
    double y, vy, fy;
   // double z, vz, fz;
    double r;
    double theta;
}tpatom;

// particle's properties
extern tpatom *atom;

typedef struct{
    double x, y, r;
}tpboundary;

// boundary particle' properties
extern tpboundary *b_par;

void alloc_atom(void);
void alloc_boundary_particle(int n);
void gen_rand_con(int iseed);
void add_dispersity_boundary_particle( int iseed );

#endif

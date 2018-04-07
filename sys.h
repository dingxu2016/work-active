#ifndef SYS_H
#define SYS_H

#define ratio 1.0

typedef struct {        
    int     natom;
    int     b_natom;
    double  phi;
    double  potential;
    double  stress;
    double  pressure;
    double  shear_ratio;
}tpsys;

extern tpsys sys;

typedef struct {
    double x, xinv;
    double y, yinv;
//    double z, zinv;
}tpbox;

extern tpbox box;

#endif

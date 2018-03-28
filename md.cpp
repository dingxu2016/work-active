#include<cstdio>
#include<cstdlib>
#include<cmath>
#include"sys.h"
#include"config.h"
#include"list.h"
#include"mathe.h"
#include"md.h"

void cal_force(void){
    double xij, yij, rij, dij, fr;
    double rx, ry, rr;
    double xx1, xx2, yy1, yy2, fix_x, fix_y;
    int j;

    for(int i=0; i<sys.natom; i++){
        atom[i].fx = 0.0;
        atom[i].fy = 0.0;
    }

    for(int i=0; i<sys.natom - 1; i++) {
        rx = atom[i].x;
        ry = atom[i].y;
        rr = atom[i].r;
        for(int k=0; k<=countn[i]; k++) {
            j = nl[i][k];
            xij = (rx - atom[j].x) * box.xinv ;
            xij = ( xij - (double)round(xij) ) * box.x ;
            yij = ry - atom[j].y;
            rij = sqrt(xij * xij + yij * yij);
            dij = rr + atom[j].r;
            if(rij < dij){
                fr = kspring * (dij - rij) / rij;
                //fr = pow( 1 - rij / dij , alpha - 1.0) / dij / rij; 

                atom[i].fx += fr * xij;
                atom[j].fx -= fr * xij;
                atom[i].fy += fr * yij;
                atom[j].fy -= fr * yij;
            }
        }
    }
    
    fix_x = box.x / 2.0;
    fix_y = box.y / 2.0;
    
    for(int i=0; i<sys.natom; i++)
    {
        yy1 = fabs(atom[i].y - fix_y);
        yy2 = fabs(atom[i].y + fix_y);
        if( yy1 < atom[i].r )
            atom[i].fy -= kspring * (0.5 - yy1);
        if( yy2 < atom[i].r )
            atom[i].fy += kspring * (0.5 - yy2);
    }

}

void move(double v0, double dt, double x1, double x2, int iseed){
    for(int i=0; i<sys.natom; i++){
        atom[i].vx = v0 * cos(atom[i].theta) + mob * atom[i].fx;
        atom[i].vy = v0 * sin(atom[i].theta) + mob * atom[i].fy;

        atom[i].x += atom[i].vx * dt;
        atom[i].y += atom[i].vy * dt;

        atom[i].theta += Gaussian_noise(x1, x2, iseed) * dt;
    }
}

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
    
    //体系粒子与边界粒子的相互作用
    for(int i=0; i<sys.natom; i++)
    {
        for(int j=0; j<sys.b_natom; j++)
        {
            xij = ( atom[i].x - b_par[j].x ) * box.xinv ;
            xij = ( xij - (double)round(xij) ) * box.x ;
            yij = atom[i].y - b_par[j].y ;
            rij = sqrt(xij * xij + yij * yij) ;
            dij = atom[i].r + b_par[j].r ;
            if(rij < dij)
            {
                fr = kspring * (dij - rij) / rij;
                atom[i].fx += fr * xij;
                atom[i].fy += fr * yij;
            }
        }
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

void boundary_move( double shear_velocity )
{
    double s_v;

    s_v = shear_velocity ;
    //添加上下边界流动条件
    //上边界向右流动，下边界向左流动
    //流动性：位置移动 |s_v| 的大小
    //调整边界粒子坐标 以 达到periodic条件
    for(int i=0; i<sys.b_natom; i++)
    {
        if( b_par[i].y > 0 ){
            b_par[i].x += s_v ;
            b_par[i].x -= (double)round( b_par[i].x / box.x ) * box.x;
        }
        else{
            b_par[i].x -= s_v ;
            b_par[i].x -= (double)round( b_par[i].x / box.x ) * box.x;
        }
    }

}

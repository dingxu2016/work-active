#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<random>
#include"sys.h"
#include"config.h"
#include"mathe.h"

tpatom *atom;
tpboundary *b_par;

void alloc_atom(void){
    atom = (tpatom *)malloc(sys.natom * sizeof(tpatom));
}

void alloc_boundary_particle( int n )
{
    b_par = (tpboundary *)malloc( 2 * n * sizeof(tpboundary) );
}

void gen_rand_con(int iseed){
    //initiate atom's radius
    for(int i=0; i<sys.natom; i++){
        if(i < sys.natom / 2)
            atom[i].r = 0.5;
        else
            atom[i].r = 0.5 * ratio;
    }

    //calculate box's length
    double sdisk = 0.0;
    for(int i=0; i<sys.natom; i++)
        sdisk += PI * atom[i].r * atom[i].r;
    double volsquare = sdisk / sys.phi;
    double temp = sqrt(volsquare);
    box.x = temp;
    box.y = temp;
//    box.z = temp;
    box.xinv = 1.0 / temp;
    box.yinv = 1.0 / temp;
//    box.zinv = 1.0 / temp;

    //initiate randomly atom's position
    /******** C的随机数库写的，不过这个随机数库性能太差**********
    srand(iseed);           
    for(int i=0; i<sys.natom; i++){
        atom[i].x = ( (double)rand()/RAND_MAX - 0.5) * box.x;
        atom[i].y = ( (double)rand()/RAND_MAX - 0.5) * box.y;
    }
    */

    for (int i=0; i<sys.natom; i++)
    {
        atom[i].x = ( random_number(iseed) - 0.5 ) * box.x ;
        atom[i].y = ( random_number(iseed) - 0.5 ) * (box.y - 2.0);
    }

    //initiate randomly atom's angle
    for(int i=0; i<sys.natom; i++){
        if(i<sys.natom/2)
            atom[i].theta = random_number(iseed) * PI;
        else
            atom[i].theta = atom[i-sys.natom/2].theta + PI; 
    }

    add_dispersity_boundary_particle( iseed );
}

void add_dispersity_boundary_particle( int iseed )
{
    int n;
    double upper_len = 0.0;
    double down_len  = 0.0;

    n = floor( box.x / atom[1].r / 2.0 ); //单个边界上排列的粒子的个数
    sys.b_natom = 2 * n ;   
    alloc_boundary_particle( n );
    for (int i=0; i<2*n; i++)
    {
        b_par[i].r = dispersity_number(iseed); //分布边界粒子半径, 分散度参看dispersity_number()函数实现
        if ( i < n )
            upper_len += b_par[i].r * 2.0;
        else
            down_len  += b_par[i].r * 2.0;
    }

    //rescale 粒子的半径
    for (int i=0; i<n; i++)
        b_par[i].r *= box.x / upper_len;
    for (int i=n; i<2*n; i++)
        b_par[i].r *= box.x / down_len;

    //给出边界粒子坐标
    for (int i=0; i<n; i++)
        b_par[i].y = box.y / 2.0;
    for (int i=n; i<2*n; i++)
        b_par[i].y = - box.y / 2.0;
    b_par[0].x = - box.x / 2.0 + b_par[0].r;
    b_par[n].x = - box.x / 2.0 + b_par[n].r;
    for (int i=1; i<n; i++)
        b_par[i].x = b_par[i-1].x + b_par[i-1].r + b_par[i].r;
    for (int i=n+1; i<2*n; i++)
        b_par[i].x = b_par[i-1].x + b_par[i-1].r + b_par[i].r;

}

#include<cstdio>
#include<cstdlib>
#include<cmath>
#include"sys.h"
#include"config.h"
#include"list.h"
#include"md.h"
#include"mathe.h"

int main(void){
    double dt, v0;
    int iseed;
    long int step;
    char str[80], str1[80], str2[80], str3[80], str4[80], s[80];
    double test_v;

    scanf("%le", &dt);
    scanf("%le", &v0);
    scanf("%d", &iseed);
    scanf("%d", &sys.natom);
    scanf("%le", &sys.phi);
    scanf("%ld", &step);
    scanf("%le", &sys.shear_ratio);

    check_cal_cluster = false;
    sigma = sqrt( 2.0 * Dr );

    alloc_atom();
    gen_rand_con(iseed);
    sprintf(str, "initiation_%.2f_%d_%d.txt", sys.phi, sys.natom, iseed);
    FILE *fp = fopen(str, "w+");
    for(int i=0; i<sys.b_natom; i++)
        fprintf(fp, "%26.16e\t%26.16e\t%26.16e\n", b_par[i].x, b_par[i].y, b_par[i].r);
    for(int i=0; i<sys.natom; i++)
        fprintf(fp, "%26.16e\t%26.16e\t%26.16e\n", atom[i].x, atom[i].y, atom[i].r);
    fclose(fp);

    int tj = 0; //count the  times about call meke_list

    alloc_list();
    make_list();
    cal_force();
    for (int i = 0; i < step; i++){
        boundary_move( sys.shear_ratio );
        move(v0, dt, mean, sigma, iseed);
        //test for bug
/*        test_v = 0.0;
        for(int csv = 0; csv < sys.natom; csv++)
        {
            test_v = max( test_v, sqrt( atom[csv].vx * atom[csv].vx + atom[csv].vy * atom[csv].vy ) );
        }
        if(i%100 == 0)
            printf("max_v = %26.16e\n", test_v);    */
        //test for bug
        if(check_list()){
            make_list();
            tj += 1;
            printf("tj = %d\ti = %d\n", tj, i);
        }
        cal_force();

        if( i > 19999 ){
            if( (i - 20000) % 2000 == 0){
                printf("i=%d\n", i);
                cal_cluster_size();
                sprintf(str1, "snapshot_config_%.2f_%d_%d_%d_data.txt", sys.phi, sys.natom, iseed, i);
                //sprintf(str2, "cluster_size_%.2f_%d_%d_%d_data.txt", sys.phi, sys.natom, iseed, i);
                FILE *fp1 = fopen(str1, "w+");
                //FILE *fp2 = fopen(str2, "w+");
                for (int j = 0; j < sys.b_natom; j++)
                    fprintf(fp1, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", b_par[j].x, b_par[j].y, b_par[j].r, 0.0, 0.0);
                for (int j = 0; j < sys.natom; j++){
                    atom[j].x -= round(atom[j].x * box.xinv) * box.x;
                    //atom[j].y -= round(atom[j].y * box.yinv) * box.y;
                    fprintf(fp1, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r, atom[j].vx, atom[j].vy);
                    //fprintf(fp2, "%d\t%d\n", j, cluster[j]);
                }
                fclose(fp1);
                //fclose(fp2);
            }   
        }
        //if( i >= 2 && i<=100 ){
        //    sprintf(str1, "snapshot_config_%.2f_%d_%d_%d_data.txt", sys.phi, sys.natom, iseed, i);
        //    FILE *fp1 = fopen(str1, "w+");
        //    for (int j = 0; j < sys.b_natom; j++)
        //        fprintf(fp1, "%26.16e\t%26.16e\t%26.16e\n", b_par[j].x, b_par[j].y, b_par[j].r);
        //    for (int j = 0; j < sys.natom; j++){
        //        atom[j].x -= round(atom[j].x * box.xinv) * box.x;
        //        //atom[j].y -= round(atom[j].y * box.yinv) * box.y;
        //        fprintf(fp1, "%26.16e\t%26.16e\t%26.16e\n", atom[j].x, atom[j].y, atom[j].r);
        //    }
        //    fclose(fp1);
        //}

    }
    
    sprintf(str3, "final_%.2f_%d_%d.txt", sys.phi, sys.natom, iseed);
    FILE *fp3 = fopen(str3, "w+");
    for(int i=0; i<sys.b_natom; i++)
        fprintf(fp3, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", b_par[i].x, b_par[i].y, b_par[i].r, 0.0, 0.0);
    for(int i=0; i<sys.natom; i++){
        atom[i].x -= round(atom[i].x * box.xinv) * box.x;
        //atom[i].y -= round(atom[i].y * box.yinv) * box.y;
        fprintf(fp3, "%26.16e\t%26.16e\t%26.16e\t%26.16e\t%26.16e\n", atom[i].x, atom[i].y, atom[i].r, atom[i].vx, atom[i].vy); 
    }
    fclose(fp3);
    //printf("L=%26.16e", box.x);

    sprintf(str4,"system_%.2f_%d_%d.txt", sys.phi, sys.natom, iseed);
    FILE *fp4 = fopen(str4, "w+");
    fprintf(fp4, "%26.16e\t%26.16e\n", box.x, box.y);
    fclose(fp4);

    free(atom);
    free(b_par);
    free(old_pos);
    free(countn);
    free(nl);
    if ( check_cal_cluster == true)
        free(cluster);

    return 0;
}

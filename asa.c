#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "main.h"

#define STEPS 30        //# tentativi ad una temperatura
#define ATTEMPTS 30     //# temperature

int param;
double * asa (double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList) {


    double * conf_best = calloc(dim, sizeof(double));
    
//    printa(conf_old,dim);
  
    int steps, attempts;
    
    printf("\n\nPlease insert the number of temperature to try: \n");
  //  scanf("%d",&steps);
    steps = 1;
    printf("Please insert the number of attempts for each temperature: \n");
 //   scanf("%d",&attempts);
    attempts = 1;
    printf("Select the numbers of changeable parameters (max: %d):\n",dim);
  //  scanf("%d", &param);
    param = 9;
    
    double cost_old = pow((*function)(argList, conf_old, dim),-1.);
    double cost_best = cost_old;
    double coeff = 1.;
    double  random, delta_cost, cost_new = 0., rapp = 1.;
    double c , cf, prob, x, c0;
    double a = 0.95;
    int i = 0, j=0;
    int count, win;
    c0 = cost_old/30;
    cf = c0/60;//da input

    printf("start cost = %f\n",cost_old*100);

  //  printa(conf_old,dim);

    double sens = 5;

    for (j = 0; j<dim; j++) {
        conf_best[j] = conf_old[j];
    }
    cost_best = cost_old;
    c = c0;
    time_t timer;
    (void) time (&timer);
    srand((long) timer);
    x = 0.;
    do {
        win = 0.;
        count = 0.;
        //	rapp = 1.;



        do {

//            printf("\n------------------------------------------\n");

            perturb(conf_old, conf_new, r_min, r_max, dim, rapp);

            cost_new = 1./((*function)(argList ,conf_new, dim));
            delta_cost = cost_new - cost_old;


//            printf("\ncost_new = %f cost_old = %f\n", cost_new, cost_old);


            prob = exp((-delta_cost)/c);//
            random = (double) rand()/(double)RAND_MAX;


            if (cost_new > 0.0 ) {
                if (delta_cost<0. || (prob < 1. && prob > random)) { //Così è da algoritmo!
                    win++;
                    //	printf("\nWinner!");
                    for (j = 0; j<dim; j++) {
                        conf_old[j] = conf_new[j];
                    }
                    cost_old = cost_new;
                    if (cost_new<cost_best) {
                        for (j = 0; j<dim; j++) {
                            conf_best[j] = conf_new[j];
                        }
                        cost_best = cost_new;
                    }
                }


            }

            //	printf("\nperc win = %f",rapp);
            count ++;

        } while (count<steps);//rapp!=1
        rapp = ((double) win/ (double) steps);
        x++;
        c = (c0+1) - pow(c0+1-cf,x/attempts);
	printf("\n**************************************\nx = %f c = %f \ncost best %f ", x,c, cost_best*100);

    } while (c>cf);
    printf("\n%f ", cost_best*100);
    printa(conf_best, dim);
    return conf_best;
}

void perturb (double * conf, double * conf_out, double * r_min, double * r_max, int dim, double rapp) {

    double window[dim];
    double random, new_val,r, sens;
    int i, index;

       
   
    
    index = rint(param*((double) rand()/ (double) RAND_MAX));

    if (index == param) index = param - 1;

    //printf("\nindex = %d\n",index);
    
    
    for (i = 0; i < dim; i++) {

        window[i] = (1.1)-pow(0.1,rapp);
        r = 2*(((double) rand()/(double)RAND_MAX) - 0.5);

       // random = r * window[i] * r_max[i];
	random = r * window[i] * (fabs(r_max[i]-r_min[i]));
        if (i == index) {
            new_val = conf[i]+random;
 //           printf("\n\nnew_val = %f  random = %f rmax[%d] = %f rmin = %f\n",new_val, random, i, r_max[i], r_min[i]);
        }
        else new_val = conf[i];



        if (new_val>r_max[i]) {

            conf_out[i] = 2*r_max[i] - new_val;
        }
        else {
            if (new_val<r_min[i]) {

                conf_out[i] = 2*r_min[i] - new_val;
            }
            else conf_out[i] = new_val;
        }

    }
//    printf("\nwindow = %f\n",window[index]);

 //printa(r_min,12);

}


void printa(double * conf,int dim) {
    int i=0;
    for (i=0; i<dim; i++) {
        printf("\nconf [%d] = %f", i, conf[i]);
    }
    printf("\n");
}






#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "main.h"

#define STEPS 30        //# tentativi ad una temperatura
#define ATTEMPTS 30     //# temperature

int parametri, FFindex[6];
double * FF_asa (double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList) {

    parametri = 6;
    double * conf_best = calloc(dim, sizeof(double));
    int steps, cicli, attempts[parametri];

    printa(conf_old,dim);
    printf("\n\nPlease insert the number of temperature to try: \n");
  //  scanf("%d",&steps);
    steps =100;
    printf("60\n");
    printf("Please insert the number of attempts for each temperature: \n");
 //   scanf("%d",&attempts);
    cicli = 30;
    printf("Select the numbers of changeable parameters (max: %d):\n",dim);
  //  scanf("%d", &param);
    
    double cost_old = 1./((*function)(argList ,conf_old, dim));//   pow((*function)(argList, conf_old, dim),-1.);
    double cost_best = cost_old;
    double vinte = 0.;
    double  random, delta_cost, cost_new = 0.;
    double  rapp[parametri], passo[parametri];
    double temp , cf, prob, c0, dum1_temp, dum2_temp;
    double vdum3, vdum4, pushup, pushdn, cutrap = 1.e-4, vi = 1.2, vf = 1.05; 
    int ii = 0, jj=0, xx;
    int count, quanti, win[parametri];

    for (ii=0; ii< parametri; ii++) rapp[ii] = 1.;

    c0 = 0.0005;
    cf = 5.e-6  ;//da input

    printf("start cost = %f\n",cost_old);

  //  printa(conf_old,dim);

    double sens = 5;
    int NOT = 0;
    for (jj = 0; jj<dim; jj++) {
        conf_best[jj] = conf_old[jj];
    }
    if (NOT) return conf_best;

    cost_best = cost_old;
    temp = cf; // Expected final T*
    dum2_temp = c0; // Starting T*
    dum1_temp = pow((temp/dum2_temp),(1.0/(float)(steps-1.)));
    temp = dum2_temp/dum1_temp; 
    vdum3 = pow(vf/vi,(1.0/(float)(steps-1.)));
    vdum4 = vi/vdum3;
    time_t timer;
    (void) time (&timer);
    srand((long) timer);

    printf("Siamo dentro ASA\n");
    printa(conf_old, dim);
    printa(r_min, dim);
    printa(r_max, dim);

    for (ii=0; ii< parametri; ii++){
        rapp[ii] = 1.;
        passo[ii] = (r_max[ii] - r_min[ii])/5.;
    }
    for(xx=1; xx<=steps; xx++) {
        temp = temp * dum1_temp;
        vdum4 = vdum4 * vdum3;
        pushdn = 2.0 - vdum4;
        pushup = vdum4;
        count = 0.;
        vinte =0.0;
        for (ii=0; ii< parametri; ii++){
           win[ii] = 0.;
           attempts[ii] = 0;
        }

        do {

            //printf("\n------------------------------------------\n");
            //(void) time (&timer);
            //printf("%s\n",asctime(localtime(&timer)));

            quanti = FF_perturb(conf_old, conf_new, r_min, r_max, dim, rapp, passo);
    //printf("\n Dopo Perturb, Temp = %f\n",temp);
    //for (ii=0; ii < parametri; ii++) printf("FFindex[%d] = %d\n",ii,FFindex[ii]); 

            for(ii=0; ii < quanti; ii++)  attempts[FFindex[ii]] ++;
            cost_new = 1./((*function)(argList ,conf_new, dim));
            delta_cost = cost_new - cost_old;


            //printf("\ncost_new = %f cost_old = %f\n", cost_new, cost_old);


            prob = exp((-delta_cost)/temp);//
            random = (double) rand()/(double)RAND_MAX;

            if (cost_new > 0.0 ) {
                if (delta_cost<0. || (prob < 1. && prob > random)) { //Così è da algoritmo!
                    for(ii=0; ii < quanti; ii++)  win[FFindex[ii]] ++;
                    //	printf("\nWinner!");
                    for (jj = 0; jj<dim; jj++) {
                        conf_old[jj] = conf_new[jj];
                        cost_best = cost_new;
                        conf_best[jj] = conf_new[jj];
                    }
                    cost_old = cost_new;
                    vinte ++;
                    /*if (cost_new<cost_best) {
                        cost_best = cost_new;
                        for (jj = 0; jj<dim; jj++) {
                            conf_best[jj] = conf_new[jj];
                        }
                    }*/
                    
                }
            }

            count ++;

        } while (count< cicli);//rapp!=1
        for (ii=0; ii< parametri; ii++){
           if (attempts[ii] == 0) rapp[ii]=0.5;
           else {
               rapp[ii] = win[ii]/attempts[ii];
               if ((rapp[ii] < 0.5) && (passo[ii]*pushdn > cutrap)){
                 passo[ii] = passo[ii] * pushdn;
               }
               else if((rapp[ii] > 0.5) && (passo[ii] * pushup < (r_max[ii] - r_min[ii])/3.)){
                 passo[ii] = passo[ii] * pushup;
               }
           }
        }
    
        //rapp = ( win/ (double) cicli);
        //
      	//printf("\nperc win = %f",rapp);
        //printf("\n**************************************\nx = %f c = %f\n",x,c);
        printf("\n cost_best=%f cost_now=%f steps =%d temp=%f  perc_win=%f \n", cost_best,cost_new,xx,temp,vinte/((float) count)); 
        printa(conf_old, dim);

    } 
    printf("\nEsco da ASA -----> cost best = %f \n", cost_best);
    printa(conf_best, dim);
    return conf_best;
}

int FF_perturb (double * conf, double * conf_out, double * r_min, double * r_max, int dim, double * rapp, double * passo) {

    double random, new_val,ran_extract, sens;
    int  ii, jj, quanti=0;
    while (quanti == 0){  
         quanti = (int) ((parametri+1)*((double) rand()/ (double) RAND_MAX));
         if (quanti == parametri+1) quanti =0;
    }
   
    if (quanti == parametri) {
      for (ii=0; ii < quanti; ii++) FFindex[ii] = ii;
    }
    else { 
      for (ii=0; ii < quanti; ii++){
         FFindex[ii] = (int) ((parametri)*((double) rand()/ (double) RAND_MAX));
            if (FFindex[ii] == parametri) ii=ii-1;
            else {
              for (jj=0; jj < ii; jj ++) 
                 if (FFindex[ii] == FFindex[jj]) ii=ii-1; 
            }
      }
    }
    //for (ii=quanti; ii < parametri; ii++) FFindex[ii] = 100;
    //printf("\n Quanti = %d\n",quanti);

    for (ii = 0; ii < dim; ii++) {
       conf_out[ii]=conf[ii];
    }
    for (ii=0; ii < quanti; ii++) {
        ran_extract = 2*(((double) rand()/(double)RAND_MAX) - 0.5);

        new_val =  conf[FFindex[ii]] + ran_extract * passo[FFindex[ii]];
        conf_out[FFindex[ii]] = new_val;

        if (new_val>r_max[FFindex[ii]]) {

            conf_out[FFindex[ii]] = r_min[FFindex[ii]] 
                                  + (new_val-r_max[FFindex[ii]]);
        }
        else if (new_val<r_min[FFindex[ii]]) {

            conf_out[FFindex[ii]] = r_max[FFindex[ii]] 
                                  - (r_min[FFindex[ii]]-new_val);  
        }
    }
  //  for (ii=0; ii < parametri; ii++)  printf("ii = %d conf_in= %f r_min = %f conf_out = %f r_max = %f passo = %f \n"
  //                                           ,ii,conf[ii],r_min[ii],conf_out[ii],r_max[ii], passo[ii]); 


// printa(conf_out,12);

}


void FF_printa(double * conf,int dim) {
    int i=0;
    for (i=0; i<dim; i++) {
        printf("\nconf [%d] = %f", i, conf[i]);
    }
    printf("\n");
}






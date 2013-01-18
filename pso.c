#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "main.h"


#define PARTICLES 60


double * pso (double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList){
	      
	_swarm swarm;
	
	int particles, attempts,  param;

	printf("\n\nPlease insert the number of swarm particles (min: 50) \n");
    	//scanf("%d",&particles);
	particles = 80;
    	printf("Please insert the number of swarm movements: \n");
    	//scanf("%d",&attempts);
	attempts = 100;
    	printf("Select the numbers of changeable parameters (max: %d):\n",dim);
   	 //scanf("%d", &param);
	param = 6;



	double * rP = calloc(dim, sizeof(double));
	double * rG = calloc(dim, sizeof(double));
	
	double alpha, cost_new;
	
	double * phyP = calloc(1, sizeof(double));
	double * phyG = calloc(1, sizeof(double));
	
	int i;

	printa(conf_old, dim);
	
	allocSwarm(&swarm, particles, dim);
	
	initSwarm(&swarm, conf_old, r_min, r_max, param, dim);
	
	initCosts(&swarm,  (*function), argList);
	
	int count = 0, num;
	
	double mean = 0.;
	
	while(count < attempts){
	 
	  mean = 0.;
	  num = 0.;
	  
	  for (i = 0; i < swarm.nOfParticles; i++){
	    
	    updateParticlePosition(&(swarm.swarmParticles[i]), r_min, r_max, swarm.bestSwarmPosition, param);
	    
	  //  if (i == 1) {
//	      printa(swarm.swarmParticles[i].position, dim);
//	      printa(swarm.swarmParticles[i].velocity, dim);
	    
	    cost_new = 
		pow((*function)(argList, swarm.swarmParticles[i].position, swarm.swarmParticles[i].nOfParameters),-1);
	    
		printf("----------------------------------------");
		printf("\nswarm %d, cost %lf\n", i, cost_new);

	    if (cost_new < swarm.bestSwarmCost) {
		  swarm.bestSwarmCost = cost_new;
		  updateBestPostion(swarm.swarmParticles[i].position, swarm.bestSwarmPosition, 
				    swarm.swarmParticles[i].nOfParameters);
		  
	    }
	    
	    if (cost_new < swarm.swarmParticles[i].bestParticleCost){
		  swarm.swarmParticles[i].bestParticleCost = cost_new;
		  updateBestPostion(swarm.swarmParticles[i].position, swarm.swarmParticles[i].bestPosition, 
				    swarm.swarmParticles[i].nOfParameters);
	    
	    }
	    
	    if (cost_new < INFINITY) {
		mean += cost_new; num++;}
	    
	  }
	  
	//  printa(swarm.bestSwarmPosition, dim);
	  printf("\nmean = %lf\n", mean/num);
	  printf("\nbest = %lf\n", swarm.bestSwarmCost);
	  
	  count++;
	  
	}

      return swarm.bestSwarmPosition;

}

void initSwarm(_swarm * swarm, double * conf_in, double * r_min, double * r_max, int param, int dim){
  
	int n, i, index;
	
	double r, num, pos, window;
	
	double vmax, vmin, vel, sens, new_pos;
	
	time_t timer;
	(void) time (&timer);
	srand((long) timer);
	
	
	for (n = 0; n < swarm->nOfParticles; n++){
	  
	
	  
	  /*First param parameters where initiated randomly*/
	  for (i = 0; i < param; i++){
	      
	      if (n == 0) swarm->swarmParticles[n].position[i] = conf_in[i];
	      
	      else{
	    
		
		/*Setting up the speed of the particle*/
		r = 2*(((double) rand()/(double)RAND_MAX) - 0.5);
	      
		vmax = r_max[i] - r_min[i];
		vmin = - vmax;
	      
	      
	      
		swarm->swarmParticles[n].velocity[i] = vmax * r;
	    
	      
	      
		/*Setting up the position of the particle*/
		index = rint(param*((double) rand()/ (double) RAND_MAX));

		if (index == param) index = param - 1;
		
				
		r = 2*(((double) rand()/(double)RAND_MAX) - 0.5);
		window = abs(r_max[i] - r_min[i]) * 0.25;	      
		
		num = r * window;	      
		
		pos = num + (r_max[i] - window);
		if (index != i) pos = 0;
		  
		  printf("pos = %lf\n", pos);
		  new_pos = swarm->swarmParticles[0].position[i] + pos;
	      
		if (new_pos>r_max[i]) {

		    swarm->swarmParticles[n].position[i] = 2*r_max[i] - new_pos;
		}
		else {
		    if (new_pos<r_min[i]) {

		      swarm->swarmParticles[n].position[i] = 2*r_min[i] - new_pos;
		    }
		    else swarm->swarmParticles[n].position[i] = new_pos;
		
		
		}
	      }   
	    //  if (n == 0) swarm->swarmParticles[n].position[i] = conf_in[i];
	    
	  }
	  for (i = param; i < dim; i++){
	  
	      swarm->swarmParticles[n].position[i] = conf_in[i];
	  
	  }
	
	    	
	}
	
    
}

void initCosts(_swarm * swarm, double (* function)(void * param, double * conf_old, int dim),
                void * argList){

	  int i;
	  
	  for (i = 0; i < swarm->nOfParticles; i++){
	      
	      swarm->swarmParticles[i].bestParticleCost =  
		pow((*function)(argList, swarm->swarmParticles[i].position, swarm->swarmParticles[i].nOfParameters),-1);

	      printf("\nswarm %d, init_cost %lf\n", i, swarm->swarmParticles[i].bestParticleCost);
//	      printa(swarm->swarmParticles[i].position, swarm->swarmParticles[i].nOfParameters);

	      updateBestPostion(swarm->swarmParticles[i].position, swarm->swarmParticles[i].bestPosition, 
				    swarm->swarmParticles[i].nOfParameters);
	//      printa(swarm->swarmParticles[i].bestPosition, 12);
		
	  }
	
	  swarm->bestSwarmCost = swarm->swarmParticles[0].bestParticleCost;
	  updateBestPostion( swarm->swarmParticles[0].position, swarm->bestSwarmPosition,
				swarm->swarmParticles[0].nOfParameters);
	
	
	for (i = 1; i < swarm->nOfParticles; i++){
	  if (swarm->swarmParticles[i].bestParticleCost < swarm->bestSwarmCost){
	      swarm->bestSwarmCost = swarm->swarmParticles[i].bestParticleCost;
	      updateBestPostion( swarm->swarmParticles[0].position, swarm->bestSwarmPosition,
				swarm->swarmParticles[0].nOfParameters);
	    
	  }
	  
	}
//	printf("swarm best position\n");
	//printa(swarm->bestSwarmPosition, 12);
}

void updateBestPostion (double * position_in, double * position_out, int dim){

    int i;
    for (i = 0; i < dim; i++){
      position_out[i] = position_in[i];
    
    }

}


void updateParticlePosition(_particle * particle, double * r_min, double * r_max, double * swarmBest, int param){

  int i;
  //int dim = particle->nOfParameters;
  
  
  double alpha, new_val;
  double rP, phyP;
  double rG, phyG;
  
  
  rP = rG = 2.;
  
  alpha = .9;
  
   int index = rint(param*((double) rand()/ (double) RAND_MAX));

   if (index == param) index = param - 1;
    
  
  
  for (i = 0; i < param; i++){
    
    phyG = (((double) rand()/(double)RAND_MAX));	     
    phyP = (((double) rand()/(double)RAND_MAX));	     
    
    particle->velocity[i] = 
	  alpha * particle->velocity[i] +
	  rP * phyP * (particle->bestPosition[i] - particle->position[i]) +
	  rG * phyG * (swarmBest[i] - particle->position[i]);
    
	  
	if (i == index) 
	  new_val = particle->position[i] + particle->velocity[i];
	
	else 
	  new_val = particle->position[i];
	
	if (new_val > r_max[i]) {

            particle->position[i] = 2*r_max[i] - new_val;
        }
        else {
            if (new_val < r_min[i]) {

                particle->position[i] = 2*r_min[i] - new_val;
            }
            else particle->position[i] = new_val;
        }
	
	
  
  }
 // printa(particle->position, particle->nOfParameters);
  
}

void allocSwarm (_swarm * swarm, int particles, int dim){

	swarm->nOfParticles = particles;
	
	swarm->swarmParticles = calloc(particles, sizeof(_particle));
	
	swarm->bestSwarmPosition = calloc(dim, sizeof(double));
		
	int i;
	
	for (i = 0; i < particles; i++){
	    swarm->swarmParticles[i].nOfParameters = dim;
	    swarm->swarmParticles[i].position = calloc(dim, sizeof(double));
	    swarm->swarmParticles[i].bestPosition = calloc(dim, sizeof(double));
	    swarm->swarmParticles[i].velocity = calloc(dim, sizeof(double));
	
	}
}

void allocParticles (_particle * particle, int dim){
  
    particle->nOfParameters = dim;
    
    particle->position = calloc(dim, sizeof(double));
    particle->bestPosition = calloc(dim, sizeof(double));
    particle->velocity = calloc(dim, sizeof(double));


}

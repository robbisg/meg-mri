typedef struct {
  double * position;
  double * bestPosition;
  double * velocity;
  int nOfParameters;
  double bestParticleCost;
} _particle;

typedef struct {
  _particle * swarmParticles;
  double * bestSwarmPosition;
  int nOfParticles;
  double bestSwarmCost;
  
} _swarm;



double * pso (double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList);
void initSwarm(_swarm* swarm, double* conf_in, double* r_min, double* r_max, int param, int dim);
void initCosts(_swarm * swarm, double (* function)(void * param, double * conf_old, int dim),
                void * argList);
void updateBestPostion (double * position_in, double * position_out, int dim);
void updateParticlePosition(_particle* particle, double* r_min, double* r_max, double* swarmBest, int param);
void allocSwarm (_swarm * swarm, int particles, int dim);

void allocParticles (_particle * particle, int dim);







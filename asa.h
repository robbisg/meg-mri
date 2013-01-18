/*functions implemented in asa.c*/
//double * asa (double * conf_in, double * r_min, double * r_max, int dim, double (* function)(double * conf_in, int dim), double * conf_out);
double * asa (double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
	      double * conf_new, void * argList);

void perturb (double * conf, double * conf_out, double * r_min, double * r_max, int dim, double step);
void printa(double * conf,int dim);

/*functions implemented in asa.c*/
//double * asa (double * conf_in, double * r_min, double * r_max, int dim, double (* function)(double * conf_in, int dim), double * conf_out);
double * FF_asa (double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
	      double * conf_new, void * argList);

int FF_perturb (double * conf, double * conf_out, double * r_min, double * r_max, int dim, double * rapp, double * passo);
void FF_printa(double * conf,int dim);

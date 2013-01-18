/*functions implemented in preProcessing.c*/
void centerMass (p_size** image, int dimx, int dimy, double* res, double* coord);
void centerMassVolume (double** centerMass, int slices, double* res, double* cdm_volume);
int * foundEdges (p_size ** image, int * dim);
void inertiaMatrix (int slice, double * cdmv, int flag, double * matrix, p_size *** image, double* ,int ** edgesVec);
void calCentMassVolume (metaImage * image);
void jacobi (double *aa, int n, double d[],  int *nrot, double *autovet);
void groupImage (metaImage * imageHigh, metaImage * imageLow, int * ratio, int* step, int*offset);
metaImage isotropize (metaImage * image, double * compareRes);
p_size filter (metaImage * image, int * coord, int * factor);

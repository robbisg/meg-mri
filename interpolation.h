/*Function implemented in interpolation.c*/
int trilinearInterp(double * coord, metaImage * image);
int nearestNeighborInterp (double * coord, metaImage * image);
int triCubicInterp(double * coord, metaImage * image);
double hcub (double arg);
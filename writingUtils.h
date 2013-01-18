/*functions implemented in writingUtils.c*/
void writeStart (double * conf, int dim);
void writeFinal (double * conf, int dim);
void writeMatrix (int* coordPre, int* coordPost, p_size*** imagePre, p_size*** imagePost);
void diffMatrix(p_size ** imageR, p_size ** imageF, p_size ** diff);
void diffImage (metaImage * mImagePre, metaImage * mImagePost, metaImage * mImageDiff);
void buildCubeVert (double** cubeVert, int** edgesMax, int slices);

void inverseMappingNew (metaImage * imageFloat, metaImage * imageTarget, double * conf, 
			  void (* action)(double * coordPre, double * coordPost, metaImage * imageFloat, metaImage * imageTarget, 
					  int (* interpolation)(double * coord, metaImage * image)),
			int (* interpolation)(double * coord, metaImage * image));		 
void findMaxEdges (int** edges, int dim, int** edgesMax);
void findTraslMaxEdges (double ** cubeVertTrasl, int dim, double ** newEdges);
void testPixel(int k, int i, int j, int z, int x, int y, p_size *** imageOrig, p_size *** imageRotata, double test);
double * calRoiCenter(int** roi, double* resolution);
void fullfil (metaImage * imageBig, metaImage * imageSmall);
int ** changeROIunit (int ** roi, double * resolution);
unsigned long int * contaPixel (metaImage* mimagePre, metaImage* mimageRot);
double ** traslCubeVert(int ** edgesVec, double * cdmVolume, int slices, double * angoli, double * res );
void inverseMappingNewLight (metaImage * imageFloat, metaImage * imageTarget, double * conf, 
			  void (* action)(double * coordPre, double * coordPost, metaImage * imageFloat, metaImage * imageTarget, 
					  int (* interpolation)(double * coord, metaImage * image)),
			int (* interpolation)(double * coord, metaImage * image));
			
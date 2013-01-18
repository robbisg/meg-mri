/*functions implemented in entropy.c*/

double mutualInformation (metaImage * lista, double* conf, int dim);
double logarithm (int f, int r);
double probability (int f, int r);
double entropy (unsigned long int ** ext_histogram);
void buildHistogram (int* coordPre, int* coordPost, p_size*** imagePre, p_size*** imagePost);
void metaBuildHistogram (double* coordOrig, double* coordRot, metaImage* imageOrig, metaImage* imageRot,
			 int (* interpolation)(double * coord, metaImage * image));
double normalizedEntropy (unsigned long int ** ext_histogram);
void metaWriteMatrix (double * coordOrig, double * coordRot, metaImage * imageOrig, metaImage * imageRot,
		      int (* interpolation)(double * coord, metaImage * image));
void normalizeHistogram(unsigned long int ** histogram, p_size ** output, unsigned long int count, int bit);

void initStruct (int slices, metaImage * image);
void initImage (int slices, int dimx, int dimy, metaImage * mImage, int background);
void setResolution (metaImage * image, double * resolution);
void setROI(metaImage * image, int Xmin, int Ymin, int Zmin, int Xmax, int Ymax, int Zmax);
void testResolutionProcedures(metaImage * imagePre, metaImage * imagePost);
void loadImages(metaImage * image, char * fileList, int slices);

void setConfiguration(double * angoli, double * trasl, double * zoom, double * shear, double * sens,
		      double * r_min, double * r_max, double * conf_in, int dim);
void cpVertex (metaImage * image, int ** vertex);

void setOrientation (metaImage * image, char orient);
metaImage * cor2trans (metaImage * image); 
metaImage * cor2sagit (metaImage * image); 
void sag2trans (metaImage* image);
metaImage * sag2coron (metaImage * imageS);
void tra2coron (metaImage * image); 
void tra2sagit (metaImage * image); 
void printInfo (metaImage * image);

void deInit (metaImage * image);
void deInit2 (metaImage image);



typedef struct{
  double * value;
  double * rMin;
  double * rMax;
} configuration;

typedef struct{
  configuration * type;
} conf;
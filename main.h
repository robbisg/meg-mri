#define MAX_file_len 800
#define MAX_file_to_read 900
#define BIT 256 //Numero di bit per codificare il livello di grigio.
#define THRESHOLD 12 //Soglia che mi discrimina il bianco e il nero (MIN is Black!)
#define MAX_trovati 1.5e6
#define X_max 980
#define X_min  70
#define Y_max 985 
#define Y_min  55
#define POST 1
#define WRITE 1
#define PRE 0
#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

#define RES(r, p) (r)*((p)+0.5) //(r)*((p)+0.5)
#define iRES(r, p) ((p)/(r))-0.5 //((p)/(r))-0.5

#define _MAX(a, b) (a > b ? a : b)

#define X 0
#define Y 1
#define Z 2

#define MIN 0
#define MAX 1

#define BACKGROUND 0
#define IMAGE 255

#define PARAM 12

#define p_size unsigned char


typedef struct {
  
  char orientation;
  char  * name;
  int dimx;
  int dimy;
  int slices;
  double * cdmVolume;
  int ** vertex;
  int ** roi;
  double * angoli;
  double * resolution;
  long int * histogram;
  p_size *** image;
  
} metaImage;



#include "debug.h"
#include "metaImage.h"
#include "asa.h"
#include "entropy.h"
#include "interpolation.h"
#include "preProcessing.h"
#include "fileReader.h"
#include "transformation.h"
#include "tiffReader.h"
#include "writingUtils.h"
#include "segmentation.h"

#include "imageUtils.h"
#include "FF_asa.h"
#include "pso.h"



/*functions implemented in cambio_coord.c*/
void rotX (int x, int y, int z, double theta, int * vec);
void rotY (int x, int y, int z, double theta, int * vec);
void rotZ (int x, int y, int z, double theta, int * vec);

/*functions implemented in main.c*/
double *  calcoloAngoli(metaImage* image, int** edges, int flag);
void scriviParam (double * conf, double * r_min, double * r_max, int dim, int flag);
void caricaDati (int flag, metaImage* image);
void scriviDati (int flag, metaImage* image);
void advanceProgressPercentage(int num);
double * minimize (double * (*minFunction)(double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList), double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList);

int isEqual(double * var_1, double * var_2, int length);

void testResolutionProcedures(metaImage * imagePre, metaImage * imagePost);

void testPlaneChanges();
void testInverseMapping ();
void testSegmentation();

extern int ** dimPre;
extern int ** dimPost;

extern int ** edgesVecPre;
extern int ** edgesVecPost;

extern char *files_pre;
extern char *files_post;
extern p_size *** imagePre;
extern p_size *** imagePost;
extern p_size *** imagePreRot;
extern p_size *** imagePostRot;
extern p_size *** imagePreRotFR;
extern p_size *** imagePostRotFR;
extern double * cdmVolumePre;
extern double * cdmVolumePost;
extern int slicesPre;
extern int slicesPost;
extern double * angoliPost;
extern double * angoliPre;
extern double * pixelSize;

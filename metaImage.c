#include <stdio.h>
#include <stdlib.h>
#include <tiffio.h>
#include <math.h>
#include "main.h"


void initStruct (int slices, metaImage * image){

    	int i;
	image->image = calloc (slices, sizeof(p_size **));
	image->name = calloc(50, sizeof(char));
	image->cdmVolume = calloc (3, sizeof(double));
	
	image->angoli = calloc(3, sizeof(double));
	
	image->resolution = calloc (3, sizeof(double));
	
	image->histogram = calloc (BIT, sizeof(long int));
	
	image->vertex = calloc (2, sizeof (int *));
	  image->vertex[0] = calloc (3, sizeof(int));
	  image->vertex[1] = calloc (3, sizeof(int));
	  
	  image->roi = calloc(2, sizeof(int *));
	  for (i = 0; i<2; i++){
	    image->roi[i] = calloc(3, sizeof(int));
	  }

}

void initImage (int slices, int dimx, int dimy, metaImage * mImage, int background){
  
      int i,j,k;
      
      initStruct(slices, mImage);
      
      mImage->dimx = dimx;
      mImage->dimy = dimy;
      mImage->slices = slices;
      
      p_size *** image = mImage->image;
      
 //     image = calloc (slices, sizeof(p_size **));
	for ( k = 0; k<slices; k++){
		image[k] = calloc(dimx, sizeof (p_size *));
		for (i = 0; i < dimx; i++){
			image[k][i] = calloc (dimy, sizeof(p_size ));
			for (j = 0; j<dimy; j++){
				image[k][i][j] = background;//Rendo lo spazio allocato bianco!! MIN is Black
				
			}
			
		}
		
	}
	
//	mImage->image = image;
	
}
void setResolution (metaImage * image, double * resolution){
  
      image->resolution[0] = resolution [0];
      image->resolution[1] = resolution [1];
      image->resolution[2] = resolution [2];
      
      
      
}

void setROI(metaImage * image, int Xmin, int Ymin, int Zmin, int Xmax, int Ymax, int Zmax){  
	 
	  
	  image->roi[0][0] = Xmin;
	  image->roi[0][1] = Ymin;
	  image->roi[0][2] = Zmin;
	  image->roi[1][0] = Xmax;
	  image->roi[1][1] = Ymax;
	  image->roi[1][2] = Zmax;
	  
	 
}

void deInit (metaImage * image){

      
    p_size ** image2D;
  
  
    free(image->cdmVolume);
     
    free(image->angoli);
         
    int i,j;
    free(image->name);

    for (i = 0; i < 2; i++){
      free(image->vertex[i]);
    }
    free(image->vertex);
#ifdef MEMORYDBG      
    system("vmstat 1 1");
#endif
    for (i = 0; i < image->slices; i++){
      
      image2D = image->image[i];
      
      for (j = 0; j < image->dimx; j++){
	free(image2D[j]);
      }
      free(image2D);
    }
    
    free(image->image);
	  
    free(image->resolution);
	   //free(&(image->resolution));
	  
    free(image->histogram);
	  
    for (i = 0; i<2; i++){
      free(image->roi[i]);
    }
    free(image->roi);
	   
#ifdef MEMORYDBG      
    system("vmstat 1 1");
#endif
}

void setOrientation (metaImage * image, char orient){
  
  image->orientation = orient;
  
  return;

}

metaImage * cor2trans (metaImage * imageS) {
  double conf[12];
  double dim[3], dimNew[3];
  int dimInt[3];
  metaImage imageTemp;
  int i,j,k;
  
  dim[0] = imageS->dimx;
  dim[1] = imageS->dimy;
  dim[2] = imageS->slices;
  
  double matrixOut[9];
  
  for(i = 0; i < 12; i++){
    if (i<9) matrixOut[i] = 0;
    conf[i] = 0;
  
  }
  
  //Rotation around Y of 90 degrees (positive is counterclockwise)
  conf[4] = -M_PI*0.5;
  
  //Zooming should be 1
  conf[6] =  conf[7] = 1;
  conf[8] = -1;
 
  
  forwardMatrix(matrixOut, conf);
  
  transformLight(dim, matrixOut, dimNew, conf);
  
 for(i =0; i<3; i++){
   dimInt[i] = rint(abs(dimNew[i]));
 
 }
 
 
  printf("%d	%d	%d\n",dimInt[2], dimInt[0], dimInt[1]);
//  printf("%f %f %f\n",imageTemp->resolution[2], imageTemp->resolution[0], imageTemp->resolution[1]);
  
  
  initImage(dimInt[2], dimInt[0], dimInt[1], &imageTemp, BACKGROUND);
  setResolution(&imageTemp, imageS->resolution);
  
  inverseMappingNewLight(imageS, &imageTemp, conf, metaWriteMatrix, trilinearInterp);
  
  metaImage * pointer = &imageTemp;
  
//  deInit(imageS);
  
  return pointer;



}

metaImage * cor2sagit (metaImage * imageS) {
double conf[12];
  double dim[3], dimNew[3];
  int dimInt[3];
  metaImage imageTemp;
  int i,j,k;
  
  dim[0] = imageS->dimx;
  dim[1] = imageS->dimy;
  dim[2] = imageS->slices;
  
  double matrixOut[9];
  
  for(i = 0; i < 12; i++){
    if (i<9) matrixOut[i] = 0;
    conf[i] = 0;
  
  }
  
  //Rotation around X of 90 degrees (positive is counterclockwise)
  conf[5] = -M_PI*0.5;
  
  //Zooming should be 1
  conf[6] =  conf[7] = 1;
  conf[8] = 1;
 
  
  forwardMatrix(matrixOut, conf);
  
  transformLight(dim, matrixOut, dimNew, conf);
  
 for(i =0; i<3; i++){
   dimInt[i] = rint(abs(dimNew[i]));
 
 }
 
 
  printf("%d	%d	%d\n",dimInt[2], dimInt[0], dimInt[1]);
//  printf("%f %f %f\n",imageTemp->resolution[2], imageTemp->resolution[0], imageTemp->resolution[1]);
  
  
  initImage(dimInt[2], dimInt[0], dimInt[1], &imageTemp, BACKGROUND);
  setResolution(&imageTemp, imageS->resolution);
  
  inverseMappingNewLight(imageS, &imageTemp, conf, metaWriteMatrix, trilinearInterp);
  
  metaImage * pointer = &imageTemp;
  
//  deInit(imageS);
  
  return pointer;


}

metaImage * sag2coron (metaImage * imageS) {

  double conf[12];
  double dim[3], dimNew[3];
  int dimInt[3];
  metaImage imageTemp;
  int i,j,k;
  
  dim[0] = imageS->dimx;
  dim[1] = imageS->dimy;
  dim[2] = imageS->slices;
  
  double matrixOut[9];
  
  for(i = 0; i < 12; i++){
    if (i<9) matrixOut[i] = 0;
    conf[i] = 0;
  
  }
  
  //Rotation around X of -90 degrees (positive is counterclockwise)
  conf[5] = -M_PI*0.5;
  
  //Zooming should be 1
  conf[6] = conf[8] = 1;
  conf[7] = 1;
  
 
  
  forwardMatrix(matrixOut, conf);
  
  transformLight(dim, matrixOut, dimNew, conf);
  
 for(i =0; i<3; i++){
   dimInt[i] = rint(abs(dimNew[i]));
 
 }
 
 
  printf("%d	%d	%d\n",dimInt[2], dimInt[0], dimInt[1]);
//  printf("%f %f %f\n",imageTemp->resolution[2], imageTemp->resolution[0], imageTemp->resolution[1]);
  
  
  initImage(dimInt[2], dimInt[0], dimInt[1], &imageTemp, BACKGROUND);
  setResolution(&imageTemp, imageS->resolution);
  
  inverseMappingNewLight(imageS, &imageTemp, conf, metaWriteMatrix, trilinearInterp);
  
  metaImage * pointer = &imageTemp;
  
//  deInit(imageS);
  
  return pointer;
  
}

void sag2trans (metaImage * image) {
  
  
}

void tra2coron (metaImage * image) {
  
  
}

void tra2sagit (metaImage * image) {
  
  
}

void printInfo (metaImage * image){
  
  printf("Image name: %s\n", image->name);
  
  printf("dimx: %d, dimy: %d, slices: %d\n",image->dimx, image->dimy, image->slices);
  
  printf("Center of mass: x %lf y %lf z %lf\n", image->cdmVolume[0], image->cdmVolume[1], image->cdmVolume[2]);
  
  printf("Resolution: x %lf y %lf z %lf\n", image->resolution[0],image->resolution[1],image->resolution[2]);
  
}




void cpVertex (metaImage * image, int ** vertex){

	int i,j;
	
	for (i = 0; i < 2; i++){
	  for (j = 0; j < 3; j++){
	    image->vertex[i][j] = vertex[i][j];
	  }
	}
	
	return;

}

void initConf (conf * configuration, int dim){
  
      configuration->type = calloc(dim, sizeof(double));
      int i;
      for (i = 0; i < dim; i++){
	configuration->type->rMax = calloc(3, sizeof(double));
	configuration->type->rMin = calloc(3, sizeof(double));
	configuration->type->value = calloc(3, sizeof(double));
      }
}
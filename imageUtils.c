#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"



metaImage * contrast (metaImage * image, int max){
  int i,j,k;
  double pixel;
  for (k = 0; k < image->slices; k++){
    for (i = 0; i < image->dimx; i++){
      for (j = 0; j < image->dimy; j++){
	pixel = (double)image->image[k][i][j]/(double)max;
	
	if (pixel > 1) pixel = 1;
	
	image->image[k][i][j] = pixel * 255;
	
      }
    }
  }
  
  
}
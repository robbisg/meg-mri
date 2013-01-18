#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"

int trilinearInterp(double * coord, metaImage * image){

  
    int dimx = image->dimx;
    int dimy = image->dimy;
    int dimz = image->slices;
  
  
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
  
  
    double xd, yd, zd;
    double xi, yi, zi;
    
    
    xd = modf(x, &xi);
    yd = modf(y, &yi);
    zd = modf(z, &zi);
    
    
    
    
    double i1, i2, j1, j2;
    
    int fz = (int) floor(z);
    int cz = (int) ceil(z);
     int fx = (int) floor(x);
    int cx = (int) ceil(x);
     int fy = (int) floor(y);
    int cy = (int) ceil(y);
    
    if (fz >= 0 && cz < dimz && fx >= 0 && cx < dimx
		&& fy >= 0 && cy < dimy){
    
	i1 = ((double) image->image[fz][fx][fy]) * (1 - zd) +
	    ((double) image->image[cz][fx][fy]) * (zd);
	  
	i2 = ((double) image->image[fz][fx][cy]) * (1 - zd) +
	    ((double) image->image[cz][fx][cy]) * (zd);
	    
	j1 = ((double) image->image[fz][cx][fy]) * (1 - zd) +
	      ((double) image->image[cz][cx][fy]) * (zd);
	  
	j2 = ((double) image->image[fz][cx][cy]) * (1 - zd) +
	  ((double) image->image[cz][cx][cy]) * (zd);
        
	double w1, w2;
    
	w1 = i1 * (1 - yd) + i2 * yd;
	w2 = j1 * (1 - yd) + j2 * yd;
    
	double value;
    
	value = w1 * (1 - xd) + w2 * xd;
    
	
	return rint(value);
	
    
    }
    else
    {
	return BACKGROUND;
    }
}

int triCubicInterp(double * coord, metaImage * image){

  
    int dimx = image->dimx;
    int dimy = image->dimy;
    int dimz = image->slices;
  
  
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
  
  
    double xd, yd, zd;
    double xi, yi, zi;
    
    
    xd = modf(x, &xi);
    yd = modf(y, &yi);
    zd = modf(z, &zi);
    
    
    
    
    double i1, i2, j1, j2;
    
    
    int cz = (int) floor(z) + 1;
    int fz = (int) floor(z);//(int) floor(z);
    int cx = (int) floor(x) + 1;
    int fx = (int) floor(x);//(int) floor(x);
    int cy = (int) floor(y) + 1;
    int fy = (int) floor(y);//(int) floor(y);
    
    
    if (fz >= 0 && cz < dimz && fx >= 0 && cx < dimx
		&& fy >= 0 && cy < dimy){
    
	i1 = ((double) image->image[fz][fx][fy]) * hcub(z - (double)fz)  +
	    ((double) image->image[cz][fx][fy]) * hcub(z - (double)cz)  ;
	  
	i2 = ((double) image->image[fz][fx][cy]) *  hcub(z - (double)fz)+
	    ((double) image->image[cz][fx][cy]) * hcub(z - (double)cz);
	    
	j1 = ((double) image->image[fz][cx][fy]) * hcub(z - (double)fz)+
	      ((double) image->image[cz][cx][fy]) * hcub(z - (double)cz);
	  
	j2 = ((double) image->image[fz][cx][cy]) * hcub(z - (double)fz) +
	  ((double) image->image[cz][cx][cy]) * hcub(z - (double)cz);
        
	double w1, w2;
    
	w1 = i1 * hcub(y - (double)fy) + i2 * hcub(y - (double)cy);
	w2 = j1 * hcub(y - (double)fy) + j2 * hcub(y - (double)cy);
    
	double value;
    
	value = w1 * hcub(x - (double)fx) + w2 * hcub(x - (double)cx);
    
	
	return rint(value);
	
    
    }
    else
    {
	/*if (x < 0. || y < 0. || z < 0. 
	  || fx >= dimx || fy >= dimy || fz >= dimz) */return BACKGROUND;
	/*else   return image->image[fz][fx][fy];*/
    }
}

int nearestNeighborInterp (double * coord, metaImage * image){

  int x = rint(coord[0]);
  int y = rint(coord[1]);
  int z = rint(coord[2]);
  
  if (z >= 0 && z < image->slices && x >= 0 && x < image->dimx
		&& y >= 0 && y < image->dimy){
  
    int value = image->image[z][x][y];
  
    return value;
  }
  else
    return BACKGROUND;

}

double hcub (double arg){
  
  double result;
  
  arg = fabs(arg);
  
  if (arg >= 0 && arg <1) result = 2.*pow(arg,3.) - 3.*pow(arg,2.) - 1.;
  else 
    result = 0;
    
  
  
  return result;
  
}
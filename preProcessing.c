
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"


void weightedVector (double x, double y, double * weight);




void inertiaMatrix (int slice, double * cdmv, int flag, double * matrix, p_size *** image, double * res, int ** edgesVec)  {

  /* Componenti della matrice di inerzia si calcolano come sommatorie:
   Ixx = SUM(i) (Yi**2 + Zi**2)    Ixy = SUM(i) -Xi*Yi      Ixz = SUM(i) -Xi*Zi 
   Iyx = SUM(i) -Xi*Yi       Iyy = SUM(i) (Xi**2 + Zi**2)   Iyz = SUM(i) -Yi*Zi
   Izx = SUM(i) -Xi*Zi       Izy = SUM(i) -Yi*Zi       Izz = SUM(i) (Xi**2 + Yi**2) */        

  int ii, kk,  i, j;  
  int bollanumero = 1;
  int partito ;
  double Ixx, Iyy, Izz;
  int *trovati_ora;
  
	double * coord = calloc(3, sizeof(double));

	p_size ** matrixFill;
	matrixFill = *(image+slice);

	

    for (i = edgesVec[slice][0]; i < edgesVec[slice][1]; i++) {
      for (j = edgesVec[slice][2]; j < edgesVec[slice][3]; j++) {
	// printf("%d ",matrix_out[yy][zz]);
         if (matrixFill[i][j] == 0) {
             Ixx=( RES(res[0], i) - *(cdmv+0))/1000;
             Iyy=( RES(res[1], j) -  cdmv[1])/1000;
	     Izz=( RES(res[2], slice) -cdmv[2])/1000;
             matrix[0] +=  Iyy*Iyy+Izz*Izz; // Ixx 
             matrix[1] +=  -Ixx*Iyy;        // Ixy
             matrix[2] +=  -Ixx*Izz;        // Ixz 
             matrix[4] +=  Ixx*Ixx+Izz*Izz; // Iyy         
             matrix[5] +=  -Iyy*Izz;        // Iyz
             matrix[8] +=  Ixx*Ixx+Iyy*Iyy; // Izz
             matrix[3]=matrix[1];
             matrix[6]=matrix[2];
             matrix[7]=matrix[5];
         }
      }
    }
  

  free(coord);
 

  return;
}

void centerMassVolume (double ** centerMass, int slices, double * resolution, double * cdm_volume){

	double totx =0. , toty = 0., totz = 0., volume = 0.; 
	double resz = resolution[2];
	
	printf("resz = %f\n", resz);
	
	int i = 0;
	for (i = 0; i < slices; i++){
		volume += *(*(centerMass+i)+2)/1000.;

	}
	
	for (i = 0; i < slices; i++){

		totx += (*(*(centerMass+i)+0) * *(*(centerMass+i)+2))/1000.;
		toty += (*(*(centerMass+i)+1) * *(*(centerMass+i)+2))/1000.;
		totz += (RES(resz,i)) * *(*(centerMass+i)+2)/1000.;
	}
	

	*(cdm_volume+0) = totx/volume;
	*(cdm_volume+1) = toty/volume;
	*(cdm_volume+2) = totz/volume;

}

void centerMass (p_size ** image, int dimx, int dimy, double * resolution, double * coord){
  /*Calcola il centro di massa sulla singola slice                             		 *
   *image = pointer to image values                                            		 *
   *dimx e dimy = dimensions of image                                  			 *
   *resolution = voxel resolution of a single element of image in millimeters   	 *
   *@return : coord = array vector that returns X and Y coordinates and MASS of the image*/
  
  
	/*Coordinate su X,Y,Z*/
	int i = 0, totx = 0, toty=0, j=0, M=0.;
	
	double resx = resolution[0] ;
	double resy = resolution[1] ;
	
	
	for (i = 0; i<dimx; i++){
		for (j = 0; j<dimy; j++){
			totx += RES(resx, i) * (image[i][j] > THRESHOLD ? 1 : 0);
			toty += RES(resy, j) * (image[i][j] > THRESHOLD ? 1 : 0);
			if (image[i][j] > THRESHOLD) M++;
			
		}
	}
	if (M > 0) {
	  *(coord+0) = (double) totx / (double) (M);
	  *(coord+1) = (double) toty / (double) (M);
	  *(coord+2) = (double) (M);
        }
        else {
          *(coord+0) = (double)(M);
          *(coord+1) = (double)(M);
          *(coord+2) = (double)(M);
        }
	//printf("\n\ncx = %f, cy = %f\n",*(coord+0),*(coord+1));
	
}

void jacobi (double *aa, int n, double d[],  int *nrot, double *autovet){

    int j,iq,ip,i;
    double tresh, theta,tau,t,sm,s,h,g,c,*b,*z;

    double v[3][3], a[3][3];


    a[0][0]=*aa;
    a[0][1]=*(aa+1);
    a[0][2]=*(aa+2);
    a[1][0]=*(aa+3);
    a[1][1]=*(aa+4);
    a[1][2]=*(aa+5);
    a[2][0]=*(aa+6);
    a[2][1]=*(aa+7);
    a[2][2]=*(aa+8);

    for (i=0; i < 9; i++)
 /*       printf("a[%d]=%f\n",i,*(aa+i));*/

    if((b = malloc(n * sizeof(double))) == NULL){
          fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
          exit(42);
    }
    if((z = malloc(n * sizeof(double))) == NULL){
          fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
          exit(42);
    }

    for (ip=0; ip < n ; ip++) {
       for (iq = 0; iq < n; iq++)
            v[ip][iq]=0.0;
       v[ip][ip]=1.0;
    }
    for (ip = 0; ip < n; ip++) {
       b[ip]=d[ip]=a[ip][ip];
       z[ip]=0.0;
    }
    *nrot=0;
    for(i = 1; i <= 50; i++) {
       sm=0.0;
       for (ip = 0; ip < n-1; ip++) {
           for (iq=ip+1; iq < n; iq++)
               sm += fabs(a[ip][iq]);
       }
       printf("sm = %f\n",sm);

       if (sm == 0.0) {
        printf("v[0][0]=%f\n",*(autovet+0)=v[0][0]);
        printf("v[0][1]=%f\n",*(autovet+1)=v[0][1]);
        printf("v[0][2]=%f\n",*(autovet+2)=v[0][2]);
        printf("v[1][0]=%f\n",*(autovet+3)=v[1][0]);
        printf("v[1][1]=%f\n",*(autovet+4)=v[1][1]);
        printf("v[1][2]=%f\n",*(autovet+5)=v[1][2]);
        printf("v[2][0]=%f\n",*(autovet+6)=v[2][0]);
        printf("v[2][1]=%f\n",*(autovet+7)=v[2][1]);
        printf("v[2][2]=%f\n",*(autovet+8)=v[2][2]);
        free(b);
        free(z);
        printf("N rot = %d and ii = %d\n",*nrot,i);
          return;
       }
       if (i < 4)
         tresh = 0.2*sm/(n*n);
       else
         tresh = 0.0  ;
       for (ip=0; ip < n-1; ip++) {
          for (iq = ip+1; iq < n; iq++){
             g = 100.0*fabs(a[ip][iq]);
       /*After four sweeps, skip the rotation if the off-diagonal element is small*/
             if (i > 4 && (double) (fabs(d[ip])+g) == (double) fabs (d[ip])
                  && (double)(fabs(d[iq])+g) == (double) fabs(d[iq]))
                a[ip][iq]=0.0;
             else if (fabs(a[ip][iq]) > tresh) {
                h = d[iq] - d[ip];
                if ((double) (fabs(h)+g) == (double)fabs(h))
                   t = (a[ip][iq])/h;           // t =1/(2 theta)
                else {
                   theta=0.5*h/(a[ip][iq]);
                   t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                   if (theta < 0.0) t= -t;
                }
                c=1.0/sqrt(1.+t*t);
                s=t*c;
                tau=s/(1.0+c);
                h=t*a[ip][iq];
                z[ip] -= h;
                z[iq] += h;
                d[ip] -= h;
                d[iq] += h;
                a[ip][iq] = 0.0;
                for (j=0;j<=ip-1;j++) {
                    ROTATE(a,j,ip,j,iq);
                }
                for (j=ip+1;j<=iq-1;j++) {
                    ROTATE(a,ip,j,j,iq);
                }
                for (j=iq+1;j<n;j++) {
                    ROTATE(a,ip,j,iq,j);
                }
                for (j=0;j<n;j++) {
                    ROTATE(v,j,ip,j,iq);     // qui modifico gli autovettori
                }
                (*nrot)++;
             }
         }
       }
       for (ip=0; ip < n; ip++) {
          b[ip] += z[ip];
          d[ip] = b[ip];
          z[ip] = 0.0;
       }
    }
	free(b); free(z);

    printf("Too many iterations in routine jacobi\n");
    return ;
}

int * foundEdges (p_size ** image, int * dim){
	/* Carina ma perfettibile*/
	int maxi = 0, maxj = 0, pix, mini = 1024, minj = 1024,i,j;

	int * edges = calloc (4, sizeof(int));
	for (i = 0; i<dim[0]; i++){
		
		for (j = 0; j<dim[1]; j++){
			
			if (image[i][j] > THRESHOLD){
								
				if (i<mini) mini = i; 
				if (i>maxi) maxi = i;
				if (j<minj) minj = j;
				if (j>maxj) maxj = j;
			}
		}
	 	   
	}
	
	*(edges+0) = mini;	
	*(edges+1) = maxi;
	*(edges+2) = minj;
	*(edges+3) = maxj;

	return edges;
	

}

void findMaxEdges (int ** edges, int dim, int ** edgesMax){
	
	int i;
	for (i = 0; i < dim; i++){
		if (edges[i][0]<edgesMax[0][0]) edgesMax[0][0] = edges[i][0];
		if (edges[i][1]>edgesMax[1][0]) edgesMax[1][0] = edges[i][1];
		if (edges[i][2]<edgesMax[0][1]) edgesMax[0][1] = edges[i][2];
		if (edges[i][3]>edgesMax[1][1]) edgesMax[1][1] = edges[i][3];
		
	}
//	printf("%d %d %d %d\n", edgesMax[0],edgesMax[1],edgesMax[2],edgesMax[3]);
	
}




/*This procedure calculates the whole image mass center*/
void calCentMassVolume(metaImage * mImage){

  
    int dimx1, dimy1, slices1;

    dimx1 = mImage->dimx;
    dimy1 = mImage->dimy;
    slices1 = mImage->slices;
    
    double * resol = mImage->resolution;
    
    p_size *** imageStack = mImage->image;
    
    int i,j,k;
    
    double ** cdMass = calloc (slices1, sizeof(double *));
    for (i = 0; i < slices1; i++){
      cdMass[i] = calloc (3, sizeof(double));
    }
    
   
    
    for (i = 0; i < slices1; i++){
      centerMass(imageStack[i], dimx1, dimy1, resol, cdMass[i]);
     
    }    
    
    centerMassVolume(cdMass, slices1, resol, mImage->cdmVolume);
    
    for (i = 0; i < slices1; i++){
      free(cdMass[i]);
    }
    free(cdMass);
    
}



/*Return the image of the same resolution of compareRes*/
metaImage isotropize (metaImage * image, double * compareRes){
    int i,j,k;
    static int flag = 0;
    
    double * imageRes = image->resolution;
    
    double count;
    
    double * ratio = calloc(3,sizeof(double));
    double * newResolution = calloc(3, sizeof(double));
    
    
    metaImage imageNew;
    
    if (imageRes[2] < compareRes[2]) {
      
      count = (double)image->slices * (double)imageRes[2] / (double)compareRes[2];
      printf("%d * %f / %f = %lf\n", image->slices, imageRes[2], compareRes[2], count);
      imageNew.slices = ceil(count);
      ratio[2] = compareRes[2] / imageRes[2];
      newResolution[2] = compareRes[2];
    }
    else {
      imageNew.slices = image->slices;
      ratio[2] = 1;
      newResolution[2] = imageRes[2];
    }
    
    if (imageRes[1] < compareRes[1]) {
      imageNew.dimy = ceil((double)image->dimy * (double)imageRes[1] / (double)compareRes[1]);
      ratio[1] = compareRes[1] / imageRes[1];
      newResolution[1] = compareRes[1];
    }
    else {
      imageNew.dimy = image->dimy;
      ratio[1] = 1;
      newResolution[1] = imageRes[1];
    }
    
    if (imageRes[0] < compareRes[0]) {
      imageNew.dimx = ceil((double)image->dimx * (double)imageRes[0] / (double)compareRes[0]);
      ratio[0] = compareRes[0] / imageRes[0];
      newResolution[0] = compareRes[0];
    }
    else {
      imageNew.dimx = image->dimx;
      ratio[0] = 1;
      newResolution[0] = imageRes[0];
    }
    
    printf("z %d x %d y %d\n", imageNew.slices, imageNew.dimx, imageNew.dimy);
    
    
    initImage(imageNew.slices, imageNew.dimx, imageNew.dimy, &imageNew, BACKGROUND);
    
    /*
    char * filename = calloc (50, sizeof(char));
    flag++;
    for (i = 0; i < imageNew.slices; i++){
	sprintf(filename,"result/%d_subsampled_%03d.tif",flag, i);
	writeTiffcustom(imageNew.image[i], filename, imageNew.dimx, imageNew.dimy);
    }
*/
    setResolution(&imageNew, newResolution);
    
    free(ratio);
//    free(filename);
    return imageNew;
    
}



void groupImage (metaImage * imageHigh, metaImage * imageLow, int * ratio, int * step, int * offset){
      
    int i, j, k;
    
    int * coord = calloc(3,sizeof(int));
    
    p_size *** image = imageLow->image;
    p_size *** imageH = imageHigh->image;
      
    for (i = 0; i < 3; i++){
	if (imageHigh->resolution[i] < imageLow->resolution[i]) {
	  ratio[i] = ceil(imageLow->resolution[i] / imageHigh->resolution[i]);
	  step[i] = floor(imageLow->resolution[i] / imageHigh->resolution[i]);
	}
	else {
	  ratio[i] = 1;
	  step[i] = 1;
	}
	
	if (offset[i] >= ratio[i]) offset[i] = ratio[i] - 1;
    }

    for (k = 0; k < imageLow->slices; k++){
	for (i = 0; i < imageLow->dimx; i++){
	  for (j = 0; j < imageLow->dimy; j++){
	    coord[0] = floor( i * step[0]) + offset[0];
	    coord[1] = floor( j * step[1]) + offset[1];
	    coord[2] = floor( k * step[2]) + offset[2];
//	    if (i < 10 && j < 10 && k < 10) printf("%d %d %d ---> %d %d %d\n",k,i,j,coord[2],coord[0],coord[1]);
	    
	    image[k][i][j] = filter(imageHigh, coord, ratio);
	//    printf("%d ", image[k][i][j]);
	  }
	}
    }
    
  free(coord);
	
}

p_size filter (metaImage * mImage, int * coord, int * ratio){
  
    int i,j,k,c;
    int factor[3];
    
    int value = 0;
    
    for (i = 0; i < 3; i++){
      factor[i] = ratio[i];
    }
    
    p_size *** image = mImage->image;
    
  //  printf("%d %d %d\n",factor[2],factor[0],factor[1]);
    
    c = 0;
    for (k = 0; k < factor[2]; k++){
      for (i = 0; i < factor[0]; i++){
	for (j = 0; j < factor[1]; j++){
	  if ((k + coord[2]) > mImage->slices-1 || (coord[0]+i) > mImage->dimx-1 ||(coord[1]+j) > mImage->dimy-1 ) break;
	  value += image[coord[2]+k][coord[0]+i][coord[1]+j];
	  c++;
	}
      }
    }
    value = rint((double)value / (double)(c));

    return (p_size)value;
}
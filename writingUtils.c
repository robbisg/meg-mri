#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"

#define nMAX(n1, n2) (n1 > n2 ? n1 : n2)

void writeFinal (/*metaImage * imagePre, metaImage * imagePost,*/ double * conf, int dim){
	
	int i,c,j,k;
	p_size ** diffM = calloc (dimPre[0][0], sizeof (p_size *));
   	for (c=0; c<dimPre[0][0]; c++){
		diffM[c] = calloc (dimPre[0][1], sizeof (p_size));
		
	}

	
	int *** tiffRead;
	char * argR = calloc (100, sizeof(char));//Da disallocare!!!
	char * argF = calloc (100, sizeof(char));
	int fettine = (slicesPost<slicesPre ? slicesPost : slicesPre);	

	for (i = 0; i<2*fettine; i++){
				
		diffMatrix(*(imagePostRot+i), *(imagePreRot+i),diffM);
		sprintf(argF, "result/slice_%03d.tif", i);
		writeTiff(imagePostRot[i], argF);
		sprintf(argF,"result/diff/diff-slice_%03d.tif",i);		
		writeTiff(diffM,argF);
	}
	/*Dealloco la memoria istanziata in changeImage poichè
	mi è servita in altre routine*/
	/*	
	for (i = 0; i<dimx; i++){
		free(*(diffM+i));
	}
	free (diffM);	
	*/

	for (i = 0; i<dimPre[0][0]; i++){
		free(*(diffM+i));
	}
	free (diffM);
	free(argR);free(argF);
	
}

void diffMatrix(p_size ** imageR, p_size ** imageF, p_size ** diff){
	int i,j;
	p_size pix;
	for (i = 0; i < 1024; i++){
		for (j = 0; j < 1024; j++){
			if (imageR[i][j] == imageF[i][j] && imageF[i][j]== BACKGROUND) pix = 0;
			if (imageR[i][j] == imageF[i][j] && imageF[i][j] == IMAGE) pix = 255;
			if (imageR[i][j] > imageF[i][j]) pix = 20;
			if (imageR[i][j] < imageF[i][j])  pix = 150;
			diff[i][j] = pix;
		}
	}
}

void diffImage (metaImage * mImagePre, metaImage * mImagePost, metaImage * mImageDiff){
      p_size ** imageR;
      p_size ** imageF;
      p_size ** imageDiff;
      
      int dimx = (mImagePost->dimx > mImagePre->dimx ? mImagePre->dimx : mImagePost->dimx);
      int dimy = (mImagePost->dimy > mImagePre->dimy ? mImagePre->dimy : mImagePost->dimy);
      
      int slices = (mImagePost->slices > mImagePre->slices ? mImagePre->slices : mImagePost->slices);
     
      
      int i,j,k;
      int perc;
      p_size pix;
      
      char * argF = malloc(100*sizeof(char));
      for (k = 0; k < slices; k++){
	      imageR = mImagePost->image[k];
	      imageF = mImagePre->image[k];
	      imageDiff = mImageDiff->image[k];
	      for (i = 0; i < dimx; i++){
		for (j = 0; j < dimy; j++){
			if (imageR[i][j] == imageF[i][j] && imageF[i][j]== BACKGROUND) pix = 0;
			if (imageR[i][j] == imageF[i][j] && imageF[i][j] == IMAGE) pix = 255;
			if (imageR[i][j] > imageF[i][j]) pix = 80;
			if (imageR[i][j] < imageF[i][j])  pix = 150;
			imageDiff[i][j] = pix;
			//imageDiff[i][j] = abs(imageR[i][j] - imageF[i][j]);
		}
	}
	
	      perc = ceil(100*((double) k / (double) mImageDiff->slices));
	      advanceProgressPercentage(perc);
	      sprintf(argF, "result/diff/diff-%03d.tif", k);
	      writeTiffcustom(mImageDiff->image[k], argF, mImageDiff->dimx, mImageDiff->dimy);
	  }
	  
	  free(argF);

}

void inverseMappingNew (metaImage * mImageFloat, metaImage * mImageTarget, double * conf,
	  void (* action)(double * coordPre, double * coordPost, metaImage *  imageFloat, metaImage * imageTarget, int (* interpolation)(double * coord, metaImage * image)),
			int (* interpolation)(double * coord, metaImage * image)){
  /*La funzione inverseMapping effettua il reverseMapping in due direzioni Forward e Rewind (leggere @rotationForward e @rotationRewind) per delucidazioni. 
   *Alla funzione va passata la chiamata alla funzione che effettua la rotazione: (rotationForward() o rotationRewind()).
   *Per come è stato concepito il programma inverseMappingForward usa prima una rotationForward e poi una rotationRewind, mentre inverseMappingRewind le usa invertite
   *per cui nel passaggio della funzione, una chiamata a rotationForward è equivalente a inverseMappingForward e di conseguenza una chiamata a rotationRewind corrispon-
   *de a un inverseMappingRewind*/
  
  /*imageFloat is the image where the transformation have to be applied*/
  /*imageTarget is the target image, so it is a void image if the function is metaWriteMatrix,
			   or it is the reference image if the function is metaBuildHistogram*/
  
  int i,j,k;
  double x,y,z;
  double ii,jj,kk;
  
 /*********************Variable declaration*********************************/
  
  int slicesTarget = mImageTarget->slices;
  int slicesFloat = mImageFloat->slices;
  
  int dimFloatX = mImageFloat->dimx; 
  int dimFloatY = mImageFloat->dimy;
  
  int dimTargetX = mImageTarget->dimx;
  int dimTargetY = mImageTarget->dimy;
  
  double * resFloat = mImageFloat->resolution;
  double * resTarget = mImageTarget->resolution;
  
  double * vettFloat = calloc (3, sizeof(double));
  double * vettTemp = calloc (3, sizeof(double));
  double * vettRoto = calloc (3, sizeof(double));
  
  double value;    					/*Foo value for edge detection*/
  double * rotPoint;					/*Rotation point for the transformation*/
  
  double * cdmVolumeFloat = mImageFloat->cdmVolume;	
  double * cdmVolumeTarget  = mImageTarget->cdmVolume;
  

	
   double ** cubeVert = calloc(8, sizeof(double *));	/*Shape containing the structure to match*/
  for (i = 0; i < 8; i++){
	*(cubeVert + i) = calloc(3, sizeof(double));	
  }
	
  static int flag = 0;
   /*********************************************************************************************/
   
    int ** roi;

      
   /*Region of interest setted on */
   if (!hasRoi(mImageFloat)) {
	  roi = mImageFloat->vertex; /*Vertex are in pixels*/
	  rotPoint = mImageFloat->cdmVolume;	/*cdm in real unit*/
	  
	  
      }
      else  {	
	
	
	  printf("Region of Interest setted\n");
	  roi = changeROIunit(mImageFloat->roi, mImageFloat->resolution); /*ROI is setted in millimeters or real unit*/
	  rotPoint = calRoiCenter(roi, mImageFloat->resolution);
      }
   //   printf("%f %f %f\n", mImageFloat->resolution[2], mImageFloat->resolution[0], mImageFloat->resolution[1]);

    /**The traslation of image vertex give us other coordinates where the image go**/ 	
   cubeVert = traslCubeVert (roi, rotPoint, mImageFloat->slices, conf, mImageFloat->resolution);
     
   double ** newEdges = calloc(2, sizeof(double *));
	
   for (i = 0; i < 2; i++){
	*(newEdges + i) = calloc(3, sizeof(double));
	for (j = 0; j < 3; j++){
		if ( i == MIN) value = 5000.;
		if ( i == MAX) value = -5000.;
			
		*(*(newEdges + i)+j) = value;

	}	
  }

  /*********************************************************************************************/
  

   findTraslMaxEdges(cubeVert, 8, newEdges);
  
  
  if (action == &metaWriteMatrix){
    for (i = 0; i < 3; i++){
      
      /**Setting up the mass center to other new images*/
      mImageTarget->cdmVolume[i] = mImageFloat->cdmVolume[i];
      mImageTarget->angoli[i] = mImageFloat->angoli[i];
      
    }
    
   
 }
 
 /**************************************************************************/
 
 int ** imageEdges = calloc(2, sizeof(int *));
 for (i = 0; i < 2; i++){
    imageEdges[i] = calloc(3, sizeof(int));
  }
  

 if(action == &metaWriteMatrix){  
 /**New vertex of the image were passed to new image (mImageTarget)**/
    for (i = 0; i < 2; i++){
      for (j = 0; j < 3; j++){
	
	imageEdges[i][j] = ceil(iRES(mImageTarget->resolution[j],
				     RES(mImageFloat->resolution[j],newEdges[i][j])));
	if (imageEdges[0][j] < 0) imageEdges[0][j] = 0;
      }
    }
  
    if (imageEdges[MAX][X] > mImageTarget->dimx) imageEdges[MAX][X] = mImageTarget->dimx;
    if (imageEdges[MAX][Y] > mImageTarget->dimy) imageEdges[MAX][Y] = mImageTarget->dimy;
    if (imageEdges[MAX][Z] > mImageTarget->slices) imageEdges[MAX][Z] = mImageTarget->slices;
   
  
    cpVertex(mImageTarget,imageEdges);
  }
 
 if (action == &metaBuildHistogram){
 
    if (hasRoi(mImageTarget)){
      roi = changeROIunit(mImageTarget->roi, mImageTarget->resolution);
      
      
    }
    else
    {
     roi = mImageTarget->vertex;
    
   //   rotPoint = mImageFloat->cdmVolume;
    
    }
     for (i = 0; i < 2; i++){
      for (j = 0; j < 3; j++){
	imageEdges[i][j] =  ceil(iRES(mImageTarget->resolution[j],
				     RES(mImageFloat->resolution[j],roi[i][j])));
	
      }
    }
 
 }
  
  rotPoint = mImageFloat->cdmVolume;
  /*************************************************************************/
  
  double * coordFloat = calloc(3, sizeof(double));
  double * coordTarget = calloc(3, sizeof(double));
  
  
  double * cdmTrasf = calloc(3, sizeof(double));
  double * cdmToTransform = calloc(3, sizeof(double));
  
  for (i=0; i<3; i++){
    
    cdmToTransform[i] = rotPoint[i] - rotPoint[i];
  }
 
  transform(cdmToTransform, cdmTrasf, conf);
  
  for (i=0; i<3; i++){
    cdmTrasf[i] += rotPoint[i];
  }
  
#ifdef DEBUG
  for (i = 0; i < 3; i++){
    printf("vertx	m %d M %d\n", mImageFloat->vertex[0][i],mImageFloat->vertex[1][i]);
    printf("roi		m %d M %d\n", roi[0][i], roi[1][i]);
    printf("newEdges	m %f M %f\n", newEdges[0][i], newEdges[1][i]);
    printf("imageEdges	m %d M %d\n", imageEdges[0][i], imageEdges[1][i]);
    printf("rotPoint	  %f\n", rotPoint[i]);
  }
#endif

  
  
 
   for (k = imageEdges[MIN][Z]; k < imageEdges[MAX][Z]; k++){
 	for (i = imageEdges[MIN][X]; i < imageEdges[MAX][X]; i++){
 		for (j = imageEdges[MIN][Y]; j < imageEdges[MAX][Y]; j++){

/*		 
 for (k = 0; k < slicesTarget; k++){
    for (i = 0; i < dimTargetX; i++){
      for (j = 0; j < dimTargetY; j++){
			 	*/		  
 	    *(vettFloat+0) = (double) RES(resTarget[0],i) - cdmTrasf[0];
 	    *(vettFloat+1) = (double) RES(resTarget[1],j) - cdmTrasf[1];
 	    *(vettFloat+2) = (double) RES(resTarget[2],k) - cdmTrasf[2];
	
	
	    transformInv(vettFloat,vettRoto, conf);
				
	    x = iRES(resTarget[0], *(vettFloat+0) + cdmTrasf[0]);
	    y = iRES(resTarget[1], *(vettFloat+1) + cdmTrasf[1]);
	    z = iRES(resTarget[2], *(vettFloat+2) + cdmTrasf[2]);

				
	    ii = iRES(resFloat[0], vettRoto[0] + cdmTrasf[0]);
	    jj = iRES(resFloat[1], vettRoto[1] + cdmTrasf[1]);
	    kk = iRES(resFloat[2], vettRoto[2] + cdmTrasf[2]);
		      
	    
	    coordFloat[0] = ii;
	    coordFloat[1] = jj;
	    coordFloat[2] = kk;
				
	    coordTarget[0] = x;
	    coordTarget[1] = y;
	    coordTarget[2] = z;
	    	
		
	    (*action)(coordFloat, coordTarget, mImageFloat, mImageTarget, interpolation);
			
      }
    } 
		
  }
	

	
	/***************************File writing*******************************/
	
	char * argF = calloc (100, sizeof(char));
	
	if (action == &metaWriteMatrix){ 
	  flag++;

	  calCentMassVolume(mImageTarget);
	  
	  int perc;
	  printf("\nWriting images to disk ...   \n");
	   for (i = 0; i < slicesTarget ; i++){
	
	     perc = ceil((double)i/(double)slicesTarget)*100;
	     advanceProgressPercentage(perc);
	     
	      sprintf(argF, "result/%d_imTransformed-%03d.tif", flag, i);
	      writeTiffcustom(mImageTarget->image[i], argF, mImageTarget->dimx, mImageTarget->dimy);
	  }
	  printf("\n");
	}
	
	
	
	/****************************************************************/
	
	
	/******************Free calls**********************************/
	if (hasRoi(mImageTarget) && action == &metaBuildHistogram){
	
	  for (i = 0; i < 2; i++){
	    free(roi[i]);
	  }
	  free(roi);
	}
	for (i = 0; i < 8; i++){
		if (i < 2) free(*(newEdges+i));
		free(*(cubeVert + i));	
	}
	
	for (i = 0; i < 2; i++){
	    
	    free(imageEdges[i]);
	}
	
	free(imageEdges);
	
		
	free(cubeVert);
	free(newEdges);
	
	free(vettFloat);
	free(vettRoto);
	free(vettTemp);
	free(coordFloat);
	free(coordTarget);
	free(cdmTrasf);
	free(cdmToTransform);	
	free(argF);
	
	return;
}
void inverseMappingNewLight (metaImage * mImageFloat, metaImage * mImageTarget, double * conf,
	 void (* action)(double * coordPre, double * coordPost, metaImage *  imageFloat, metaImage * imageTarget, int (* interpolation)(double * coord, metaImage * image)),
			     int (* interpolation)(double * coord, metaImage * image)){
  /*La funzione inverseMapping effettua il reverseMapping in due direzioni Forward e Rewind (leggere @rotationForward e @rotationRewind) per delucidazioni. 
   *Alla funzione va passata la chiamata alla funzione che effettua la rotazione: (rotationForward() o rotationRewind()).
   *Per come è stato concepito il programma inverseMappingForward usa prima una rotationForward e poi una rotationRewind, mentre inverseMappingRewind le usa invertite
   *per cui nel passaggio della funzione, una chiamata a rotationForward è equivalente a inverseMappingForward e di conseguenza una chiamata a rotationRewind corrispon-
   *de a un inverseMappingRewind*/
  
  /*imageFloat is the image where the transformation have to be applied*/
  /*imageTarget is the target image, so it is a void image if the function is metaWriteMatrix,
			   or it is the reference image if the function is metaBuildHistogram*/
  
  int i,j,k;
  double x,y,z;
  double ii,jj,kk;
  
 /*********************Variable declaration*********************************/
  
  int slicesTarget = mImageTarget->slices;
  int slicesFloat = mImageFloat->slices;
  
  int dimFloatX = mImageFloat->dimx; 
  int dimFloatY = mImageFloat->dimy;
  
  int dimTargetX = mImageTarget->dimx;
  int dimTargetY = mImageTarget->dimy;
  
  double * resFloat = mImageFloat->resolution;
  double * resTarget = mImageTarget->resolution;
  
  double * vettFloat = calloc (3, sizeof(double));
  double * vettTemp = calloc (3, sizeof(double));
  double * vettRoto = calloc (3, sizeof(double));
  
  double value;    					/*Foo value for edge detection*/
  double * rotPoint;					/*Rotation point for the transformation*/
  
  double * cdmVolumeFloat = mImageFloat->cdmVolume;	
  double * cdmVolumeTarget  = mImageTarget->cdmVolume;
  
  double matrixForward[9];
  double matrixInverse[9];
	
  for (i = 0; i < 9; i++){
    matrixForward[i] = 0.;matrixInverse[i] = 0.;	/*Matrix cleaning*/
  }
  
  
  double ** cubeVert = calloc(8, sizeof(double *));	/*Shape containing the structure to match*/
  for (i = 0; i < 8; i++){
	*(cubeVert + i) = calloc(3, sizeof(double));	
  }
	
  static int flag = 0;
   /*********************************************************************************************/
   
    int ** roi;
  
  /******************************/

 if (action == &metaWriteMatrix){ 
    conf[0] += (RES(mImageTarget->resolution[0],mImageTarget->dimx) 
	      - RES(mImageFloat->resolution[0],mImageFloat->dimx)) * 0.5;
    conf[1] += (RES(mImageTarget->resolution[1],mImageTarget->dimy)
	      - RES(mImageFloat->resolution[1],mImageFloat->dimy)) * 0.5;
    conf[2] += (RES(mImageTarget->resolution[2],mImageTarget->slices)
	      - RES(mImageFloat->resolution[2],mImageFloat->slices)) * 0.5;
 } 
 
  /**********************/
//  printa(conf, 12);
      
   /*Region of interest setted on */
   if (!hasRoi(mImageFloat)) {
	  roi = mImageFloat->vertex;
	  rotPoint = mImageFloat->cdmVolume;	
	  
	  
      }
      else  {	
	
	
	  printf("Region of Interest setted\n");
	  roi = changeROIunit(mImageFloat->roi, mImageFloat->resolution);
	  rotPoint = calRoiCenter(roi, mImageFloat->resolution);
      }
 //     printf("%f %f %f\n", mImageFloat->resolution[2], mImageFloat->resolution[0], mImageFloat->resolution[1]);

    /**The traslation of image vertex give us other coordinates where the image go**/ 	
   cubeVert = traslCubeVert (roi, rotPoint, mImageFloat->slices, conf, mImageFloat->resolution);
     
   double ** newEdges = calloc(2, sizeof(double *));
	
   for (i = 0; i < 2; i++){
	*(newEdges + i) = calloc(3, sizeof(double));
	for (j = 0; j < 3; j++){
		if ( i == MIN) value = 5000.;
		if ( i == MAX) value = -5000.;
			
		*(*(newEdges + i)+j) = value;

	}	
  }

  /*********************************************************************************************/
  

   findTraslMaxEdges(cubeVert, 8, newEdges);
  
  
 
 
 /**************************************************************************/
 
 int ** imageEdges = calloc(2, sizeof(int *));
 for (i = 0; i < 2; i++){
    imageEdges[i] = calloc(3, sizeof(int));
  }
  

 if(action == &metaWriteMatrix){  
 /**New vertex of the image were passed to new image (mImageTarget)**/
    for (i = 0; i < 2; i++){
      for (j = 0; j < 3; j++){
	imageEdges[i][j] = rint(newEdges[i][j]);
//	printf("imageEdges M %d\n",  imageEdges[1][j]);
//	printf("newEdges M %f\n",  newEdges[1][j]);
	if (imageEdges[0][j] < 0) imageEdges[0][j] = 0;
      }
    }
  
    printf("\n");
  
    if (imageEdges[MAX][X] > mImageTarget->dimx) 
      imageEdges[MAX][X] = mImageTarget->dimx;
    if (imageEdges[MAX][Y] > mImageTarget->dimy) 
      imageEdges[MAX][Y] = mImageTarget->dimy;
    if (imageEdges[MAX][Z] > mImageTarget->slices) 
      imageEdges[MAX][Z] = mImageTarget->slices;
   
  
    cpVertex(mImageTarget,imageEdges);
  }
  
 /*******************************************************************************/
 if (action == &metaBuildHistogram){
 
    if (hasRoi(mImageTarget)){
	roi = changeROIunit(mImageTarget->roi, mImageTarget->resolution);
        //roi = mImageTarget->roi;
      
    }
    else
    {
      roi = mImageTarget->vertex;
    
   //   rotPoint = mImageFloat->cdmVolume;
    
    }
     for (i = 0; i < 2; i++){
      for (j = 0; j < 3; j++){
	imageEdges[i][j] = roi[i][j];
	//imageEdges[i][j] = 2*roi[i][j];
      }
    }
 
 }
  
  rotPoint = mImageFloat->cdmVolume;
  /*************************************************************************/
  
  double * coordFloat = calloc(3, sizeof(double));
  double * coordTarget = calloc(3, sizeof(double));
  
  double * cdmTrasf = calloc(3, sizeof(double));
  double * cdmToTransform = calloc(3, sizeof(double));

  for (i=0; i<3; i++){
    
    cdmToTransform[i] = rotPoint[i] - rotPoint[i];
  }
 
 forwardMatrix(matrixForward, conf);

 transformLight(cdmToTransform, matrixInverse, cdmTrasf, conf);
//  transform(cdmToTransform, cdmTrasf, conf);
  
  for (i=0; i<3; i++){
    cdmTrasf[i] += rotPoint[i];
  }
  
#ifdef DEBUG
  for (i = 0; i < 3; i++){
    printf("resolution t	m %f\n",mImageTarget->resolution[i]);
    printf("vertx float	m %d 	M %d\n", mImageTarget->vertex[0][i],mImageTarget->vertex[1][i]);
    printf("roi	targ	m %d 	M %d\n", roi[0][i], roi[1][i]);
    printf("newEdges t	m %f 	M %f\n", newEdges[0][i], newEdges[1][i]);
    printf("imageEdges t	m %d	 M %d\n", imageEdges[0][i], imageEdges[1][i]);
    printf("rotPoint t	  %f\n\n", rotPoint[i]);
  }
#endif

  /************************************************************************/
   if (action == &metaWriteMatrix){
    for (i = 0; i < 3; i++){
      
      /**Setting up the mass center to other new images*/
      mImageTarget->cdmVolume[i] = cdmTrasf[i];
      mImageTarget->angoli[i] = mImageFloat->angoli[i];
      
    }
    
   
 }
  
  /**************************************************************************/
  inverseMatrix(matrixInverse, conf);
  
  
   for (k = imageEdges[MIN][Z]; k < imageEdges[MAX][Z]; k++){
 	for (i = imageEdges[MIN][X]; i < imageEdges[MAX][X]; i++){
 		for (j = imageEdges[MIN][Y]; j < imageEdges[MAX][Y]; j++){
 
		 
/* for (k = 0; k < slicesTarget; k++){
    for (i = 0; i < dimTargetX; i++){
      for (j = 0; j < dimTargetY; j++){*/
			 			  
 	    *(vettFloat+0) = (double) RES(resTarget[0],i) - cdmTrasf[0];
 	    *(vettFloat+1) = (double) RES(resTarget[1],j) - cdmTrasf[1];
 	    *(vettFloat+2) = (double) RES(resTarget[2],k) - cdmTrasf[2];
	
	
	    invTransformLight(vettFloat, matrixForward, vettRoto, conf);
	   // transformInv(vettFloat,vettRoto, conf);
				
	    x = iRES(resTarget[0], *(vettFloat+0) + cdmTrasf[0]);
	    y = iRES(resTarget[1], *(vettFloat+1) + cdmTrasf[1]);
	    z = iRES(resTarget[2], *(vettFloat+2) + cdmTrasf[2]);

				
	    ii = iRES(resFloat[0], vettRoto[0] + cdmTrasf[0]);
	    jj = iRES(resFloat[1], vettRoto[1] + cdmTrasf[1]);
	    kk = iRES(resFloat[2], vettRoto[2] + cdmTrasf[2]);
		      
	    
	    coordFloat[0] = ii;
	    coordFloat[1] = jj;
	    coordFloat[2] = kk;
				
	    coordTarget[0] = x;
	    coordTarget[1] = y;
	    coordTarget[2] = z;
	    	
		
	    (*action)(coordFloat, coordTarget, mImageFloat, mImageTarget,interpolation);
			
      }
    } 
		
  }
	

	
	/***************************File writing*******************************/
	
	char * argF = calloc (100, sizeof(char));
	
	if (action == &metaWriteMatrix){ 
	  flag++;

	  calCentMassVolume(mImageTarget);
	  
	  int perc;
	  
	  char * fileName = calloc( 50, sizeof(char));

	  if (flag > 0  ) {
	    for (i = 0; i < slicesTarget ; i++){
	
	      perc = ceil((double)i/(double)slicesTarget)*100;
	      advanceProgressPercentage(perc);
	     
	      sprintf(argF, "result/%d_imTransformed-%03d.tif", flag, i);
	      writeTiffcustom(mImageTarget->image[i], argF, mImageTarget->dimx, mImageTarget->dimy);
	    }
	    printf("\n");
	  }
	}
	
	
	/****************************************************************/
	
	
	/******************Free calls**********************************/
	if (hasRoi(mImageTarget) && action == &metaBuildHistogram){
	
	  for (i = 0; i < 2; i++){
	    free(roi[i]);
	  }
	  free(roi);
	}
	for (i = 0; i < 8; i++){
		if (i < 2) free(*(newEdges+i));
		free(*(cubeVert + i));	
	}
	
	for (i = 0; i < 2; i++){
	    
	    free(imageEdges[i]);
	}
	
	free(imageEdges);
	
		
	free(cubeVert);
	free(newEdges);
	
	free(vettFloat);
	free(vettRoto);
	free(vettTemp);
	free(coordFloat);
	free(coordTarget);
	free(cdmTrasf);
	free(cdmToTransform);	
	free(argF);
	
	return;
}
void collectDataWrite(){



}
double ** traslCubeVertLight(int ** edgesMax, double * cdmVolume, int slices, double * conf, double * res){

  
	int i,j;
	
	
	double * matrixF = calloc(9, sizeof(double));
	/*Alloco la memoria per memorizzare le coordinate
	del parallelepipedo da traslare*/
	double ** cubeVert = calloc(8, sizeof(double *));
	for (i = 0; i < 8; i++){
		*(cubeVert + i) = calloc(3, sizeof(double));	
	}

	double ** cubeVertTrasl = calloc(8, sizeof(double *));
	for (i = 0; i < 8; i++){
		*(cubeVertTrasl + i) = calloc(3, sizeof(double));	
	}

	/*Trovo le coordinate dei vertici del parallelepipedo*/
	buildCubeVert(cubeVert, edgesMax, slices);


	/*Sottraggo il centro di massa per far sì che ruoti attorno a quel punto*/
		
	for (i = 0; i < 8; i++){
		for (j = 0; j < 3; j++){		
			cubeVert[i][j] = RES(res[j],(cubeVert[i][j])) - cdmVolume[j];
		}
	}

	/*Effettuo la rotazione dei vertici*/
	
	forwardMatrix(matrixF, conf);
	for (i = 0; i < 8; i++){

	    transformLight(cubeVert[i], cubeVertTrasl[i], matrixF, conf);
		// transform(cubeVert[i], cubeVertTrasl[i], conf);
	
	}
	
	for (i = 0; i < 8; i++){
		for (j = 0; j < 3; j++){
			cubeVertTrasl[i][j] = iRES(res[j],(cubeVertTrasl[i][j])) + iRES(res[j],(cdmVolume[j]));
		
		}
		

	}
	
	
	for(i = 0; i < 8; i++){
	
	    free(*(cubeVert+i));
	}
	free(cubeVert);
	free(matrixF);
	
	return cubeVertTrasl;
}

double ** traslCubeVert(int ** edgesMax, double * cdmVolume, int slices, double * conf, double * res){

  
	int i,j;
	
	
	
	/*Alloco la memoria per memorizzare le coordinate
	del parallelepipedo da traslare*/
	double ** cubeVert = calloc(8, sizeof(double *));
	for (i = 0; i < 8; i++){
		*(cubeVert + i) = calloc(3, sizeof(double));	
	}

	double ** cubeVertTrasl = calloc(8, sizeof(double *));
	for (i = 0; i < 8; i++){
		*(cubeVertTrasl + i) = calloc(3, sizeof(double));	
	}

	/*Trovo le coordinate dei vertici del parallelepipedo*/
	buildCubeVert(cubeVert, edgesMax, slices);


	/*Sottraggo il centro di massa per far sì che ruoti attorno a quel punto*/
		
	for (i = 0; i < 8; i++){
		for (j = 0; j < 3; j++){		
			cubeVert[i][j] = RES(res[j],(cubeVert[i][j])) - cdmVolume[j];
		}
	}

	/*Effettuo la rotazione dei vertici*/
	for (i = 0; i < 8; i++){

		 transform(cubeVert[i], cubeVertTrasl[i], conf);
	
	}
	
	for (i = 0; i < 8; i++){
		for (j = 0; j < 3; j++){
			cubeVertTrasl[i][j] = iRES(res[j],(cubeVertTrasl[i][j])) + iRES(res[j],(cdmVolume[j]));
		
		}
		

	}
	
	
	for(i = 0; i < 8; i++){
	
	    free(*(cubeVert+i));
	}
	free(cubeVert);
	
	return cubeVertTrasl;
}

void buildCubeVert (double ** cubeVert, int ** edgesMax, int slices){

	int i,j,k;
	int x,y,z;

	for (i = 0; i < 8; i++){
		if ((i % 2) == 0){
			
			z = edgesMax[0][2];
			
			if (( i / 2) < 2) x = edgesMax[0][0];
			else x = edgesMax[1][0];
			if (( i / 2) % 2 == 0) y = edgesMax[0][1];
			else y = edgesMax[1][1];
			
		}
		else{
			z = edgesMax[1][2];
			
			if (( (i-1) / 2) < 2) x = edgesMax[0][0];
			else x = edgesMax[1][0];
			if (( (i-1) / 2) % 2 == 0) y = edgesMax[0][1];
			else y = edgesMax[1][1];
			
		}

		cubeVert[i][0] = x;
		cubeVert[i][1] = y;
		cubeVert[i][2] = z;
		
// 		for (j = 0; j < 3; j++){
// 			printf("cube vert = %f ",cubeVert[i][j]);
// 		}
// 		printf("\n");
	}

}

void findTraslMaxEdges (double ** cubeVertTrasl, int dim, double ** newEdges){
	int i,j;
	
	for (i = 0; i < dim; i++){
		for (j = 0; j < 3; j++){

			if (cubeVertTrasl[i][j] < newEdges[0][j]) newEdges[0][j] = cubeVertTrasl[i][j];
			if (cubeVertTrasl[i][j] > newEdges[1][j]) newEdges[1][j] = cubeVertTrasl[i][j];

		}
	}
}

void fullfil (metaImage * imageBig, metaImage * imageSmall){
 
    int i,j,k;
    
    int xOff, yOff, zOff;
    
    
    int offset[3];
    
    xOff = rint((double)(imageBig->dimx - imageSmall->dimx)*0.5);
    yOff = rint((double)(imageBig->dimy - imageSmall->dimy)*0.5);
    zOff = rint((double)(imageBig->slices - imageSmall->slices)*0.5);
    
    offset[0] = xOff;
    offset[1] = yOff;
    offset[2] = zOff;
    
    for (k = 0; k < imageSmall->slices; k++){
      for (i = 0; i < imageSmall->dimx; i++){
	for (j = 0; j < imageSmall->dimy; j++){
	  
	  imageBig->image[k+zOff][i+xOff][j+yOff] = imageSmall->image[k][i][j];
	  
	}
      }
      
    }
    
    for (i = 0; i < 2; i ++){
      for (j = 0; j < 3; j++){
	imageBig->vertex[i][j] = imageSmall->vertex[i][j]+offset[j];
	imageBig->cdmVolume[j] = imageSmall->cdmVolume[j]+offset[j];
      }
    }
#ifdef WRITEIMG  
    char * argF = calloc(80, sizeof(char));
     for (k = 0; k < imageBig->slices; k++){
     
      sprintf(argF, "result/fulfilled/post-fulfilled-%03d.tif", k);
      writeTiffcustom(imageBig->image[k], argF, imageBig->dimx, imageBig->dimy);
     }
    
    free(argF);
#endif    
}

int hasRoi (metaImage * image){
  
  /*Returns 1 if it's setted, and 0 otherwise*/
  
  int i,j;
  
  for (i = 0; i < 2; i++){
    for (j = 0; j<3; j++){
    
      if (image->roi[i][j] != 0) {
	return 1;
	break;
      }
    }
  }
  return 0;
  
}


int ** changeROIunit (int ** roi, double * resolution){

    
    int i, n;
    
    int ** newROI = calloc(2, sizeof(int *));
    for (i = 0; i < 2; i++){
    
      newROI[i] = calloc(3, sizeof(int *));
    }
    
    for (n = 0; n < 2 ; n++){
      for (i = 0; i < 3; i++){
    
	newROI[n][i] = rint(iRES(resolution[i], (double)roi[n][i]));
      }
    }
    return newROI;

}


double * calRoiCenter(int ** roi, double * resolution){
  
  double * result = calloc(3,sizeof(double));
  
  int i;
  
  for (i = 0; i < 3; i++){
    result[i] = (RES(resolution[i], (double)roi[1][i])  + RES(resolution[i],(double)roi[0][i]))*0.5;  
  }

  return result;
}


void testInverseMapping (){
	
	int i;

	double pointA [3] = {0, 0, 1};
	double pointB [3] = {0, 1, 0};
	double pointC [3] = {1, 0, 0};

	double angoli [3] = {0, 0, -M_PI_2};

	double * pointAT = calloc(3, sizeof(double));
	double * pointBT = calloc(3, sizeof(double));
	double * pointCT = calloc(3, sizeof(double));

	rotationForward(pointA, pointAT, angoli);
	rotationForward(pointB, pointBT, angoli);
	rotationForward(pointC, pointCT, angoli);
	
	for (i = 0; i < 3; i++){
		printf("%f %f \n", pointA[i], pointAT[i]);
	}	
	printf("\n");
	for (i = 0; i < 3; i++){
		printf("%f %f \n", pointB[i], pointBT[i]);
	}
	printf("\n");
	for (i = 0; i < 3; i++){
		printf("%f %f \n", pointC[i], pointCT[i]);
	}
	printf("\n");
	double * pointATR = calloc(3, sizeof(double));
	double * pointBTR = calloc(3, sizeof(double));
	double * pointCTR = calloc(3, sizeof(double));

	rotationRewind(pointAT, pointATR, angoli);
	rotationRewind(pointBT, pointBTR, angoli);
	rotationRewind(pointCT, pointCTR, angoli);

	for (i = 0; i < 3; i++){
		printf("%f %f \n", pointA[i], pointATR[i]);
	}	
	printf("\n");
	for (i = 0; i < 3; i++){
		printf("%f %f \n", pointB[i], pointBTR[i]);
	}
	printf("\n");
	for (i = 0; i < 3; i++){
		printf("%f %f \n", pointC[i], pointCTR[i]);
	}
	printf("\n");
	
	free(pointAT);
	free(pointBT);
	free(pointCT);
	
	free(pointATR);
	free(pointBTR);
	free(pointCTR);
}


void testPixel(int k, int i, int j, int z, int x, int y, p_size *** imageOrig, p_size *** imageRotata, double test){
	if (z >= 0 && z < slicesPre && x >= 0 && x < 1024 && y >= 0 && y < 1024)
	  if (imageOrig[z][x][y] == 0) printf("%d %d %d value = %d ---> %f\n", z, x,y, imageOrig[z][x][y],test);
}


unsigned long int * contaPixel (metaImage* mimagePre, metaImage* mimageRot){
    
    int i,j,k;
    unsigned long int pixelPre = 0, pixelRot = 0;
    
    unsigned long int * returnValue = calloc(2, sizeof(unsigned long int));
    
    int slicesPre = mimagePre->slices;
    int slicesRot = mimageRot->slices;
    
    p_size *** imagePre = mimagePre->image;
    p_size *** imageRot = mimageRot->image;
    
    
    for (k = 0; k < slicesPre; k++){
      for (i = 0; i < mimagePre->dimx; i++){
	for (j = 0; j < mimagePre->dimy; j++){
	  if (imagePre[k][i][j] == IMAGE) pixelPre++;
	}
      }
    }
    for (k = 0; k < slicesRot; k++){
      for (i = 0; i < mimageRot->dimx; i++){
	for (j = 0; j < mimageRot->dimy; j++){
	  if (imageRot[k][i][j] == IMAGE) pixelRot++;
	}
      }
    }
    printf("\nPre %ld Post %ld\n", pixelPre, pixelRot);
    
    returnValue[0] = pixelPre;
    returnValue[1] = pixelRot;
    
    return returnValue;
}


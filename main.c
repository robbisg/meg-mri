#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "main.h"


double * pixelSize;
int ** edgesVecPre;
int ** edgesVecPost;

/*DEPRECATED VARIABLES*/
int ** dimPre;
int ** dimPost;
double * angoliPost;
double * angoliPre;
p_size *** imagePre;
p_size *** imagePost;
p_size *** imagePreRot;
p_size *** imagePostRot;
p_size *** imagePreRotFR;
p_size *** imagePostRotFR;
double * cdmVolumePre;
double * cdmVolumePost;
p_size *** output;
int slicesPre;
int slicesPost;
/*DEPRECATED VARIABLES*/

const char progress[4] = {'|', '/', '-', '\\'};


int main(int argc, char * argv[]){ 

	
//    testPlaneChanges ();



//	printf("int = %d, short = %d, unsigned char = %d, char = %d\n\n",
//	  sizeof(void ***), sizeof(short int ***), sizeof(unsigned char ***), sizeof(char***));


//	exit(42);
//	testInverseMapping();

	FILE * file, *parPtr;

	

	int i,k; 	

	char * filesPre, *filesPost; 				/*Variables containing the list of images*/

	char scelta;						/*Choice option*/
 	char param_file[50], stringa_read[50]; 			/*Variables for loading data*/
 	
 	
	metaImage mImagePre, mImagePost;
	
   	
	if (((file = fopen("pre/datiPre.txt", "r")) == NULL) || ((file = fopen("post/datiPost.txt", "r")) == NULL) 
				|| ((file = fopen("lista_param_input.par", "r"))== NULL))
	{
		printf("Dati non presenti, si procedera con il calcolo.\n\n");
		scelta = 'c';
	}			
	else
	{
		printf("Dati presenti. \n\nDigita 'c' per calcolarli nuovamente o 'l' per caricarli da file.\n");
		//scanf("%c",&scelta);
	
	}	

	scelta = 'c';
	
	
	  /*These pointers let me to better control the flow of information*/			
	  metaImage * mImageFloat;
	  metaImage * mImageTarget;
	
	
	if (scelta == 'c'){
	  
	  if((filesPre = (char *) malloc(MAX_file_len * MAX_file_to_read)) == NULL){
	    fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
	    exit(42);
	  } 

	  if((filesPost = (char *) malloc(MAX_file_len * MAX_file_to_read)) == NULL){
	    fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
	    exit(42);
	  }


	  mImagePre.slices = readFile(filesPre, PRE);     	/*LF*/
	  mImagePost.slices = readFile(filesPost, POST);	/*HF*/

	  printf("slicesPre = %d\n", mImagePre.slices);
	  printf("slicesPost = %d\n", mImagePost.slices);
	
//	  mImagePre.slices = mImagePost.slices = 50;

	  initStruct(mImagePre.slices, &mImagePre);
	  initStruct(mImagePost.slices, &mImagePost);
	
	  pixelSize = calloc(3, sizeof(double ));
	
	
	  pixelSize[0] = 0.188;
	  pixelSize[1] = 0.25;
	  pixelSize[2] = 6.;
	  /*Setting up pixel resolution*/
	  setResolution(&mImagePre, pixelSize); 
	
	  pixelSize[0] = 1.;
	  pixelSize[1] = 1.;
	  pixelSize[2] = 6.;
	
	
	  /*Setting up pixel resolution*/
	  setResolution(&mImagePost, pixelSize);
	
	  free(pixelSize);
	  
	
	
				/*Calcolo dei dati*/
	
	  double lowRes[3] = {1  , 1 , 6   };
	  metaImage imageLowPre, imageLowPost;

	  
	  
	  printf("Image Loading...    \n");
	
	  /*Image loading*/
	  loadImages(&mImagePre, filesPre, mImagePre.slices);
	  loadImages(&mImagePost, filesPost, mImagePost.slices);		
	
#ifdef MEMORYDBG
	  printf("\nmemory before grouping!!!\n\n");
	  system("free -m");
#endif 
	  
	  
/***************************************************************************************************/

	  if (isEqual(lowRes, mImagePre.resolution, 3) && isEqual(lowRes, mImagePost.resolution, 3)){
	    /*Do not group*/
	    mImageTarget = &mImagePre;
	    mImageFloat = &mImagePost;
	    
	  }
	  
	  else {
	    /*Image Grouping*/
	    int * offset = calloc(3, sizeof(int));
	    int * window = calloc(3, sizeof(int));
	    int * overlap = calloc(3, sizeof(int));
	    
	    
	    offset[0] = 0;	 /*Pixel offset along x direction*/
	    offset[1] = 0;	 /*Pixel offset along y direction*/
	    offset[2] = 0;	 /*Pixel offset along z direction*/
	    
	    
	    overlap[0] = 0;	/*Overlap along x direction*/
	    overlap[1] = 0;	/*Overlap window step along x direction*/
	    overlap[2] = 0;	/*Overlap window step along x direction*/
	        
	    imageLowPre = isotropize(&mImagePre, lowRes);
	    groupImage(&mImagePre, &imageLowPre, window, overlap, offset);
	    
	    
	    offset[0] = 0;	 /*Pixel offset along x direction*/
	    offset[1] = 0;	 /*Pixel offset along y direction*/
	    offset[2] = 0;	 /*Pixel offset along z direction*/
	    
	    imageLowPost = isotropize(&mImagePost, lowRes);
	    groupImage(&mImagePost, &imageLowPost, window, overlap, offset);
	        
	    free(window);
	    free(overlap);
	    free(offset);

	    
#ifdef MEMORYDBG
	    printf("\nmemory after grouping!!!\n\n");
	    system("free -m");    
#endif	  
	    
	    deInit(&mImagePre);
	    deInit(&mImagePost);
	    
	     mImageTarget = &imageLowPre;
	     mImageFloat  = &imageLowPost;
	  }
	/*********************************************************/	
	/*******************Calculation**************************/
	/********************************************************/
	
	
	 
	
	  printf("slicesPre = %d\n\n\n\n\n", mImageTarget->slices);
	  printf("slicesPost = %d\n\n\n\n\n", mImageFloat->slices);
	
	  edgesVecPre = calloc(mImageTarget->slices, sizeof(int *));
	  edgesVecPost = calloc(mImageFloat->slices, sizeof(int *));
	  
	  dimPre = calloc(mImageTarget->slices, sizeof(int *));
	  dimPost = calloc(mImageFloat->slices, sizeof (int *));
	  
	  
	  /*Image shape extractione*/
	
	  printf("Select filtering...");
	  filterImage(mImageTarget, thresholdInt);
	  filterImage(mImageFloat, thresholdInt);
	
//	  filterImage(mImageTarget, noFilter);
//	  filterImage(mImageFloat, noFilter);
	  /******************Image tiff writing**********************************/
	  char * fileName = calloc(100, sizeof(char));

	  sprintf(fileName,"pre/pre_grouped");
	  writeMetaImage(mImageTarget, fileName);
		
	  sprintf(fileName,"post/post_grouped");
	  writeMetaImage(mImageFloat, fileName);	

	  free(fileName);
	  
	
	  /****************************Euler's Angles**************************************/
	  mImageTarget->angoli = calcoloAngoli(mImageTarget, edgesVecPre, PRE);
	  mImageFloat->angoli = calcoloAngoli(mImageFloat,  edgesVecPost, POST);
	
	
	
	  /*Queste assegnazioni servono di supporto per findMaxEdges*/
	  mImageFloat->vertex[0][0] = mImageTarget->vertex[0][0] = mImageTarget->dimx;
	  mImageFloat->vertex[0][1] = mImageTarget->vertex[0][1] = mImageTarget->dimy;
	  mImageFloat->vertex[0][2] = mImageTarget->vertex[0][2] = 0; 
	  mImageFloat->vertex[1][2] = mImageFloat->slices;
	  mImageTarget->vertex[1][2] = mImageTarget->slices;
	
	
	  /*********************************Find Shape edges******************************/
	  findMaxEdges(edgesVecPost,  mImageFloat->slices, mImageFloat->vertex);
	  findMaxEdges(edgesVecPre, mImageTarget->slices, mImageTarget->vertex);

	
	
	  scriviDati(POST, mImageFloat);
	  scriviDati(PRE, mImageTarget);	
	  
	  for (i = 0; i < mImageTarget->slices; i++) {
	    free(edgesVecPre[i]);
	    free(dimPre[i]);
	  }
	  free(edgesVecPre);
	  free(dimPre);
	  
	  for (i = 0; i < mImageFloat->slices; i++) {
	    free(edgesVecPost[i]);
	    free(dimPost[i]);
	  }
	  free(edgesVecPost);
	  free(dimPost);
      }	



        /*Memory allocation for the minimization algorithm*/
	double * conf_in = calloc(12, sizeof (double));		/*Start configuration*/
	double * conf_out = calloc(12, sizeof (double));	/*Final configuration*/
	double r_min [12]; 					/*Upper configuration limit*/
	double r_max [12]; 					/*Lower configuration limit*/

			/*Loading data from files*/
			/*Da rivedere*/
     if (scelta == 'l'){

	caricaDati(PRE,  &mImagePre); 
	caricaDati(POST, &mImagePost);	
	
	
	sprintf (param_file,"lista_param_win.par");
        if ( (parPtr = fopen(param_file, "r")) == NULL)
           printf("Errore apertura file\n");

        for (i = 0; i < PARAM; i++) {
           for (k = 1; k <= 10 ; k++) {
              fscanf(parPtr, "%s", stringa_read+k*20);
         //     printf("%s\n",stringa_read+k*20);
           }
           sscanf((stringa_read+4*20),"%16lf\n",&conf_in[i]);
           sscanf((stringa_read+7*20),"%16lf\n",&r_min[i]);
           sscanf((stringa_read+9*20),"%16lf\n",&r_max[i]);
           printf("%d %f %f %f\n",i,conf_in[i],r_min[i],r_max[i]);
        }
        fclose(parPtr);
	printa(conf_in, 12);
	printa(r_min, 12);
	printa(r_max, 12);

	mImageTarget = &mImagePre;
	mImageFloat = &mImagePost;
      }
	
      
/******************************************************************************************************/
/********************************PRIMA TRASFORMAZIONE AFFINE*******************************************/
/******************************************************************************************************/


      int j;  
      
    
      
      metaImage mImageDiff;     	/*Difference image*/
      metaImage mImagePostSized;	/*Image with same dimension of target*/
      
      metaImage mImageOrigOvers;	/*image oversampled*/


#ifdef DEBUG   
      for (i = 0; i < 2; i++){
	for (j = 0; j < 3; j++){
	  printf("pre vrtx = %d ", mImageTarget->vertex[i][j]);
	  printf("post vrtx = %d ", mImageFloat->vertex[i][j]);
	}
	printf("\n");
      }
        

     printf("\n");  
     printf("\n");
#endif
 #ifdef MEMORYDBG     
      system("free -m");
#endif
  /**************Codice da eliminare**************************************/  
  
      int dimXBig = _MAX(mImageFloat->dimx, mImageTarget->dimx);
      int dimYBig = _MAX(mImageFloat->dimy, mImageTarget->dimy);
      int dimZBig = _MAX(mImageFloat->slices, mImageTarget->slices);
      
      /*Uniformation of comparison image space*/
      initImage(dimZBig, dimXBig, dimYBig, &mImagePostSized, BACKGROUND);
      setResolution(&mImagePostSized, mImageFloat->resolution);
      
      fullfil(&mImagePostSized, mImageFloat);

 #ifdef MEMORYDBG     
      system("free -m");
#endif
	
 /***************************************************************************/     
      
#ifdef DEBUG 
      /*Printing debugging information*/
      for (i = 0; i < 2; i++){
	for (j = 0; j < 3; j++){
	  printf("%d ", mImagePostSized.vertex[i][j]);
	}
	printf("\n");
      }
#endif
   
      metaImage mImagePreRot;    /*Image Transformed*/
      metaImage mImagePostRot;	 /*Image Post Transformed*/
      
      
 #ifdef MEMORYDBG     
      system("free -m");
#endif
      
                                   
      /*initialization of images*/
      initImage(mImageTarget->slices, dimXBig, dimYBig, &mImagePreRot, BACKGROUND);
      setResolution(&mImagePreRot, mImageTarget->resolution);
      
 #ifdef MEMORYDBG     
      system("free -m");
#endif
      
      
      double * conf_in_rot = calloc(PARAM, sizeof(double));
      
      
      int dim = PARAM;
	double r_min_rot [12]; 					/*Upper configuration limit*/
	double r_max_rot [12];
      /*Setting parameters for the first transformation*/
      double shear[]  = {0., 0., 0.};//{0., 0., 0.};
      double zoom[]   = {1., 1., 1.};//{-1.3, 0.94, 0.92};
      double angoli[] = {0., 0., 0.};//{3.11, 0., 0.};
      double trasl[]  = {0., 0., 0.};
      
          
      /*Sensibility of the transformation parameters*/
      /*Traslation, angles, Zoom and shear*/
      double sens [] = {30., 1., 1., 0.9}; 
      
     
      setConfiguration(angoli, trasl, zoom, shear, sens, r_min_rot, r_max_rot, conf_in_rot, dim);
               
        
      printa(conf_in_rot, 12);
    
      /*First transformation of image, without the calculation of MI*/
      inverseMappingNewLight(mImageTarget, &mImagePreRot, conf_in_rot, metaWriteMatrix, trilinearInterp);
     
      free(conf_in_rot);
      
#ifdef DEBUG   
      for (i = 0; i < 2; i++){
	for (j = 0; j < 3; j++){
	  printf("pre vrtx = %d ", mImagePreRot.vertex[i][j]);
	  printf("post vrtx = %d ", mImagePostSized.vertex[i][j]);
	}
	printf("\n");
      }
        

     printf("\n");  
     printf("\n");
#endif
   
#ifdef MEMORYDBG     
      system("free -m");
#endif
	 
	  
	  deInit(mImageTarget);
	  
 #ifdef MEMORYDBG     
      system("free -m");
#endif
      
#ifdef DEBUG   
	  printf("\n");  
	  printf("\n");
#endif
        
      
     
 
/******************************************************************************************************/
/*********************************MINIMIZATION PROCESS*************************************************/
/******************************************************************************************************/

      //contaPixel(&mImagePreRot, &mImagePostSized);
       
            
      
   
         
      double shearT[]  = {0., 0., 0.};//{0., 0., 0.};
      double zoomT[]   = {1., 1., 1.};//{-1.3, 0.94, 0.92};
      double angoliT[] = {0., 0., 0.};//{-0.13, -0.20, 0.25};//{-0.,12 -0.2, 0.2};
      
      double traslT [] = {0., 0., 0.};//{-22., 9., -14.};//{-19., 10., -16.};//{-20., 8., 10.};
	
      double sensT [] = {0., 0., 0.0, 0.}; 
	
      for (i = 0; i < 3; i++){
	//traslT[i]  =  mImagePreRot.cdmVolume[i] - mImagePostSized.cdmVolume[i];
      } 
      
    //  double traslT2 [] = {0., 0., 0.};
      
     if (scelta == 'c') { 
	setConfiguration(angoliT, traslT, zoomT, shearT, sensT, r_min, r_max, conf_in, PARAM);
	
	scriviParam(conf_in, r_min, r_max, PARAM, PRE);
     }
      printa(conf_in, 12);
   
     
     /**************************MINIMIZATION***************************************/ 
      
      metaImage imageList[2];
      
      mImageTarget = &mImagePreRot;
      mImageFloat = &mImagePostSized;
      
      
      imageList[0] = (*mImageFloat);	/*imageFloat*/
      imageList[1] = (*mImageTarget);  	/*imageTarget*/
      
     
      imageHistogram2D(&mImagePostSized, mImagePostSized.histogram);
      imageHistogram2D(&mImagePreRot, mImagePreRot.histogram);
  
  //    setROI(&mImagePreRot, 9,19, 6, 139, 216, 30);
      
      
      /*********************************************************************/
      
      
      
      printf("\n");
      system("date");
      printf("\n");
     
      int sample = 1;
      
      double ** statConf = calloc(sample, sizeof(double *));
      for (i = 0; i < sample; i++){
	  statConf[i] = calloc(PARAM, sizeof (double));
      }
	
      for (i = 0; i < sample; i++){
	
	conf_in = minimize(asa, conf_in, r_min, r_max, dim, mutualInformation, conf_out, imageList);
	printf("run %d MI %lf\n",i, pow(mutualInformation(imageList, conf_in, dim),-1));
      
      }
    
      system("date");
      printf("\n");
     /******************************************************************************/
          
      initImage(mImagePostSized.slices, mImagePostSized.dimx, mImagePostSized.dimy, &mImageDiff, BACKGROUND);
      setResolution(&mImageDiff, mImagePostSized.resolution);
      
      
      initImage(mImagePostSized.slices, mImagePostSized.dimx, mImagePostSized.dimy, &mImagePostRot, BACKGROUND);
      setResolution(&mImagePostRot, mImagePostSized.resolution);
      
       inverseMappingNewLight(&mImagePostSized, &mImagePostRot, conf_in, metaWriteMatrix, trilinearInterp);
       printf("%lf\n\n",pow(mutualInformation(imageList, conf_in, dim),-1));
    
      diffImage(&mImagePostRot, &mImagePreRot, &mImageDiff);
      
      contaPixel(&mImagePreRot,&mImagePostSized);
     
      calcParameters(&mImagePreRot, &mImagePostRot);
      
      deInit(&mImageDiff);
      deInit(&mImagePostRot);
      deInit(&mImagePreRot);
      deInit(&mImagePostSized);
     
      
  
      scriviParam(conf_in, r_min, r_max, PARAM, POST);
      
      /*****************************************************************************/
      /*****************************************************************************/
      /**************Rotazione dell'immagine in scala di grigi**********************/
      /*****************************************************************************/
      /*****************************************************************************/
      
       /*Building of oversampled image*/
       metaImage mImageOriginal;
       
              
      printf("Calculation of original images parameters....\n");
      mImageOriginal.name = "mImageOriginal";
      mImageOriginal.slices = readFile(filesPost, POST);
      initStruct(mImageOriginal.slices, &mImageOriginal);
      loadImages(&mImageOriginal, filesPost, mImageOriginal.slices);
      filterImage(&mImageOriginal, noFilter);
      
      int ** edgesVecOrig = calloc(mImageOriginal.slices, sizeof(int *));
      
//      mImageOriginal.slices = 50;
      
      mImageOriginal.resolution[0] = 1;
      mImageOriginal.resolution[1] = 1;
      mImageOriginal.resolution[2] = 1;
      
      for (i = 0; i < 3; i++){
	mImageOriginal.vertex[0][i] = 0;
      }
      edgesVecPre = calloc(mImagePre.slices, sizeof(int *));
      
      mImageOriginal.vertex[1][0] = mImageOriginal.dimx;
      mImageOriginal.vertex[1][1] = mImageOriginal.dimy;
      mImageOriginal.vertex[1][2] = mImageOriginal.slices;
      
      mImageOriginal.angoli = calcoloAngoli(&mImageOriginal, edgesVecOrig, POST);
      
      printInfo(&mImageOriginal);
             
      
      
     // cpVertex(&mImageOriginal, mImagePostSized.vertex);
    //  printf("roi = %d \n", mImageOriginal.roi[0][0]);
   
      

    // setResolution(&mImageOriginal, res);
      
      
      double pixel [3] = {1., 1., 1.};
      
      mImageOrigOvers.dimx = rint(iRES(pixel[0],mImageOriginal.dimx));
      mImageOrigOvers.dimy = rint(iRES(pixel[1],mImageOriginal.dimy));
      mImageOrigOvers.slices = rint(iRES(pixel[2],mImageOriginal.slices));
       
      initImage(mImageOrigOvers.slices, mImageOrigOvers.dimx, mImageOrigOvers.dimy,
		&mImageOrigOvers, BACKGROUND);
      
      setResolution(&mImageOrigOvers, pixel);

    //  printa(conf_in,12);
     
      
      printf("Image Oversized!\n\n");
      inverseMappingNewLight(&mImageOriginal, &mImageOrigOvers, conf_in, metaWriteMatrix, trilinearInterp);

     
      free(conf_in);
      free(conf_out);

	exit(42);
}



double * calcoloAngoli (metaImage* image, int** edges, int flag){

	int  ii,c;	
	
	int dimx = image->dimx;
	int dimy = image->dimy;
	
	printf("%d %d\n\n\n",dimx,dimy);
	
	double * resolution = image->resolution;
	int righe = image->slices;
	
	double ** cMassMatrix = calloc (righe, sizeof (double));
	for (c =0 ;c<righe; c++){
		cMassMatrix[c] = calloc (3, sizeof(double));
	}

	ii = 0;
	
	int * dim = calloc(2,sizeof(int));
	dim[0] = dimx;
	dim[1] = dimy;
	
        FILE *fileInerzia;

	int * edge;

	/*Calcolo il centro di massa delle immagini e li salvo in una matrice cMassMatrixPre*/
	while (ii <  righe) {
	  
		

		edge = foundEdges(*(image->image+ii),dim);
		*(edges+ii) = edge;


		centerMass(*(image->image+ii),dimx,dimy,resolution,*(cMassMatrix+ii));
#ifdef DEBUG
		printf("%d %f %f %f\n", ii, *(*(cMassMatrix+ii)+0),
			*(*(cMassMatrix+ii)+1),*(*(cMassMatrix+ii)+2));
#endif

		ii++;
  	}


	
	/*Calcolo del centro di massa totale sul volume*/
	centerMassVolume (cMassMatrix, righe, resolution, image->cdmVolume);
	printf("%f %f %f\n",  *(image->cdmVolume+0),
			*(image->cdmVolume+1),*(image->cdmVolume+2));


//	double * inerzia = calloc (9, sizeof(double));//Disallocare
	double inerzia[9] = {0};
	ii=0;


	/*Calcolo della matrice del tensore d'inerzia*/
	while (ii <  righe) {
		
		inertiaMatrix(ii,image->cdmVolume,flag,inerzia,image->image, resolution,edges);	
		ii++;
  	}

	if (flag == PRE) fileInerzia = fopen("pre/InerziaPre.txt", "w");
	else fileInerzia = fopen("post/InerziaPost.txt", "w");
        
	fprintf(fileInerzia,"%lf %lf %lf\n",inerzia[0],inerzia[1],inerzia[2]);	
        fprintf(fileInerzia,"%lf %lf %lf\n",inerzia[3],inerzia[4],inerzia[5]);	
        fprintf(fileInerzia,"%lf %lf %lf\n",inerzia[6],inerzia[7],inerzia[8]);	
	fclose(fileInerzia);

	int ordine = 3;

	double * diag = calloc (3, sizeof(double));
	double * autovet = calloc(9,sizeof(double));
	int * nrot = calloc (1, sizeof(int));
	
	/*Applico jacobi e ottengo autovalori e autovettori del tensore di inerzia*/
	jacobi(inerzia, ordine, diag, nrot, autovet);

	double *v_orig_1 = calloc (3, sizeof(double));
	double *v_orig_2 = calloc (3, sizeof(double));
	double *v_orig_3 = calloc (3, sizeof(double));	


	*(v_orig_1+0) = *(autovet+0);
   	*(v_orig_1+1) = *(autovet+3);
  	*(v_orig_1+2) = *(autovet+6);
   	*(v_orig_2+0) = *(autovet+1);
  	*(v_orig_2+1) = *(autovet+4);
 	*(v_orig_2+2) = *(autovet+7);
	*(v_orig_3+0) = *(autovet+2);
   	*(v_orig_3+1) = *(autovet+5);
   	*(v_orig_3+2) = *(autovet+8);
	

	double * angoli = calloc(3,sizeof(double));

	ruotamento_on_X_Y_Z(v_orig_1,v_orig_2,v_orig_3,angoli);//DEPRECATED

	free(v_orig_1);	
	free(v_orig_2);
	free(v_orig_3);

	for (c = 0; c < righe; c++){
		free(*(cMassMatrix+c));
	}	
	free(cMassMatrix);

	
	free(diag);
	free(nrot);
	free(autovet);

	free(dim);
	
	
	
	return angoli;

}

void scriviDati (int flag, metaImage * image){
	
	FILE * file;
	
	int slices = image->slices;
	
	if (flag == PRE) file = fopen("pre/datiPre.txt", "w");
	else file = fopen("post/datiPost.txt", "w");
	
	fprintf(file,"%d\n",slices);
	fprintf(file,"%f %f %f\n",image->cdmVolume[0],image->cdmVolume[1],image->cdmVolume[2]);
	fprintf(file,"%f %f %f\n",image->angoli[0],image->angoli[1],image->angoli[2]);
	
	int i;
	for(i = 0; i < 2; i++){
		fprintf(file,"%d %d %d\n",image->vertex[i][0],image->vertex[i][1],image->vertex[i][2]);
	}
	for(i = 0; i < 3; i++){
		fprintf(file,"%lf ",image->resolution[i]);
	}
	printf("\n");
	
	fclose(file);

}

void caricaDati (int flag, metaImage * image){
	
	FILE * file;
	
	
	
	
	if (flag == PRE) file = fopen("pre/datiPre.txt", "r");
	else file = fopen("post/datiPost.txt", "r");

	printf("Loading data...    ");

	char line[1000];

	int i=0;

	while (fgets (line, 1000, file) != NULL){
		
		i++;
		switch (i){
			case 1: sscanf(line,"%d\n",&(image->slices));
				break;
			case 2: sscanf(line,"%lf %lf %lf",&(image->cdmVolume[0]),&(image->cdmVolume[1]),&(image->cdmVolume[2]));
				break;			
			case 3: sscanf(line,"%lf %lf %lf",image->angoli,&(image->angoli[1]),&(image->angoli[2]));
				break;
			case 4: sscanf(line,"%d %d %d",&(image->vertex[i-4][0]),&(image->vertex[i-4][1]),
			      				&(image->vertex[i-4][2]));
			case 5: sscanf(line,"%d %d %d",&(image->vertex[i-4][0]),&(image->vertex[i-4][1]),
			      				&(image->vertex[i-4][2]));				
			default: {
					/*edges[i-4] = calloc(4, sizeof(int));
					 sscanf(line,"%d %d %d %d\0",&(edges[i-4][0]),&(edges[i-4][1]),
			      				&(edges[i-4][2]),&(edges[i-4][3]));*/
					/*sscanf(line,"%d %d %d",&(image->vertex[i-4][0]),&(image->vertex[i-4][1]),
			      				&(image->vertex[i-4][2]));
*/
					sscanf(line,"%lf %lf %lf",&(image->resolution[0]),&(image->resolution[1]),&(image->resolution[2]));
					break;
				}

		}
		
	}

	i = 0.;
		
	char * filename = calloc (50, sizeof(char));
        p_size c;
//	int ** matrixFill;

	int ** dim = calloc(image->slices, sizeof(int *));
        for (c = 0; c < image->slices; c++){
  		*(dim+c) = calloc(2, sizeof(int));
  	}
	int perc;
	
	image->image = calloc(image->slices, sizeof(p_size ***));
	
	for (i = 0; i < image->slices; i++){
	  
		if (flag == PRE) {
			sprintf(filename,"pre/fill-slice_%03d.tif",i);
		}
		else{
			sprintf(filename,"post/fill-slice_%03d.tif",i);
		}
		
		perc = ceil(100*((double) i / (double) image->slices));
		advanceProgressPercentage(perc);
	
                *(image->image+i) = readTiff(*(dim+i), filename);
 
	}	
	image->dimx = dim[0][0];
	image->dimy = dim[0][1];
	
	
	for (c = 0; c < image->slices; c++){
	  free(dim[c]);
	}
	free(dim);
	
	
	fclose(file);

}

void scriviParam (double * conf, double * r_min, double * r_max, int dim, int flag){
	FILE *file;
	int ii = 0, index;
	if (flag == PRE) file = fopen("lista_param_input.par", "w");
	else file = fopen("lista_param_win.par", "w");
        /*
	char nome[50];
	
	for (ii = 0; ii < dim; ii++){
	    index = ii/3;
	    if (index == 0) nome = "trasl";
	    if (index == 1) nome = "angoli";
	    if (index == 2) nome = "shear";
	    if (index == 3) nome = "zoom";
	    
	  fprintf(file,"%3d) %s   =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii+1,nome,conf[ii],r_min[ii],r_max[ii]);
	}*/
	
	
	ii=1;
        fprintf(file,"%3d) delta_x    =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
        ii=2;
        fprintf(file,"%3d) delta_y    =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
        ii=3;
        fprintf(file,"%3d) delta_z    =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
        ii=4;
	fprintf(file,"%3d) alpha      =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
        ii=5;
        fprintf(file,"%3d) beta       =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
        ii=6;
        fprintf(file,"%3d) gamma      =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
	ii=7;
	fprintf(file,"%3d) zoom1      =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
	ii=8;
	fprintf(file,"%3d) zoom2      =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
	ii=9;
	fprintf(file,"%3d) zoom3      =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
	ii=10;
	fprintf(file,"%3d) shear1     =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
	ii=11;
	fprintf(file,"%3d) shear2     =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
	ii=12;
	fprintf(file,"%3d) shear3     =  %16lf  0.000000      [ %14lf - %14lf ]\n",ii,conf[ii-1],r_min[ii-1],r_max[ii-1]);
	

	fclose(file);

}

void setConfiguration(double * angoli, double * trasl, double * zoom, double * shear, double * sens,
		      double * r_min, double * r_max, double * conf_in, int dim){
      
      int i;
      for (i = 0; i < 3; i++){
	  conf_in[i+0] =  trasl[i];
	  conf_in[i+3] =  angoli[i];
	  conf_in[i+6] =  zoom[i];
	  conf_in[i+9] =  shear[i];
	  
	 	  
	  r_max[i+0] = trasl[i]+sens[0];
	  r_max[i+3] = angoli[i]+sens[1];
	  r_max[i+6] = zoom[i]+sens[2];
	  r_max[i+9] = shear[i]+sens[3];
	  
	  	  
	  r_min[i+0] = trasl[i]-sens[0];
	  r_min[i+3] = angoli[i]-sens[1];
	  r_min[i+6] = zoom[i]-sens[2];
	  r_min[i+9] = shear[i]-sens[3];
	  
      }

}

void loadImages(metaImage * image, char * fileList, int slices){
  
  int ii = 0;
  int dim[2];
  int perc;
  while (ii <  slices) {
    
		perc = ceil(100*((double) ii / (double) slices));
		advanceProgressPercentage(perc);
		image->image[ii] = readTiff(dim, fileList+ii*MAX_file_len);
		ii++;
  }

  image->dimx = dim[0];
  image->dimy = dim[1];
  
  printf("\n");
}

int ruotaAssi(metaImage * image, metaImage * imageRot){
  
  
  int i,j,k,x,y,z;
  
  int dimx = image->dimx;
  int dimy = image->dimy;
  
  int slices = image->slices;
  
  double * matrix = calloc(9, sizeof(double));
  
  double * vettO = calloc(3, sizeof(double));
  double * vettR = calloc(3, sizeof(double));
  
  matrix_X(-0.2, matrix);
	
	vettO[0] = (double) dimx;
	vettO[1] = (double) dimy;
	vettO[2] = (double) slices;
  
  righe_per_colonne(vettO, matrix, vettR);
  
	x = 128.;
	y = 128.;
	z = 228.;
	
	
	
  initImage(z,x,y,imageRot,BACKGROUND);
  
  for (k = 0; k < slices; k++){
    for (i = 0; i < dimx; i++){
      for (j = 0; j < dimy; j++){
      
	vettO[0] = (double) i;
	vettO[1] = (double) j;
	vettO[2] = (double) k;
	
	
	righe_per_colonne(vettO, matrix, vettR);
	
	x = abs(rint(vettR[0]));
	y = abs(rint(vettR[1]));
	z = abs(rint(vettR[2]));
	
	imageRot->image[z][x][y] = image->image[k][i][j];
	
      }
    }
  }
    
    char * argF = calloc (100, sizeof(char));//Da disallocare!!!
    for (i = 0; i<228; i++){
		sprintf(argF, "rotate/rotate_%03d.tif", i);
		writeTiffcustom(imageRot->image[i], argF, imageRot->dimx, imageRot->dimy);
	}  
      
  return 0;
}

double * minimize (double * (*minFunction)(double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList), double * conf_old, double * r_min, double * r_max, int dim, double (* function)(void * param, double * conf_old, int dim),
              double * conf_new, void * argList){

	double * conf_out = (*minFunction)(conf_old, r_min, r_max, dim, function, conf_new, argList);
	
	return conf_out;

}

void calcParameters(metaImage * mImagePre, metaImage * mImagePost){
  
  p_size *** imagePre = mImagePre->image;
  p_size *** imagePost = mImagePost->image;
  int i;
  
  short int ** borders = calloc(mImagePost->dimx, sizeof(short int *));
  for (i = 0; i < mImagePost->dimx; i++){
    borders[i] = calloc(2, sizeof(short int ));
  }
  
  int x,y,z;
  double pixMatch = 0, pixLFWhite = 0;
  
  metaImage provaBorders;
  
  initImage(mImagePost->slices, mImagePost->dimx, mImagePost->dimy, &provaBorders, BACKGROUND);
  
  for (z = 0; z < mImagePost->slices; z++){
    
    findBorders(mImagePost, borders, z);
    
    for (x = 0; x < mImagePost->dimx; x++){
      for (y = borders[x][0]; y < borders[x][1]; y++){
	  
	
	if (mImagePre->image[z][x][y] == IMAGE){
	  pixLFWhite++;
	  if (mImagePost->image[z][x][y] == mImagePre->image[z][x][y])  
	    pixMatch++;
	}
	
	
	provaBorders.image[z][x][y] = IMAGE;

	}
     }
    
    writeMetaImage(&provaBorders, "result/borders");
    
  }
  
  printf("\n\nMatching white pixels LF: %lf\n\n", pixMatch/pixLFWhite);
  
  
  
  
  
}

void findBorders (metaImage * mImage, short int ** borders, int slice){
  
  p_size ** image = mImage->image[slice];
  
  int i,j;
  int white = 0;
  p_size value = 0;
  
  for (i = 0; i < mImage->dimx; i++){
    white = 0;
    for (j = 0; j < mImage->dimy; j++){
      
      value = image[i][j];
      if (value == IMAGE && white == 0) {
	white = 1;
	borders[i][0] = j;
      }
      if (value == IMAGE && white == 1){
	borders[i][1] = j;
      }
      
    }
 //   printf("x %d ymin %d ymax %d \n", i, borders[i][0], borders[i][1]);
  }
  
 
  
  
}

void advanceProgressPercentage(int num) {
  printf("\r%d  %%",num);
  fflush(stdout);
}



void testPlaneChanges(){

	char * filesPre, *filesPost; 	/*Variables containing the list of images*/

	metaImage mImagePre, mImagePost;
	

	if((filesPre = (char *) malloc(MAX_file_len * MAX_file_to_read)) == NULL){
          fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
          exit(42);
        } 

	if((filesPost = (char *) malloc(MAX_file_len * MAX_file_to_read)) == NULL){
          fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
          exit(42);
        }

	
	mImagePre.slices = readFile(filesPre, PRE);     /*TAC*/
	mImagePost.slices = readFile(filesPost, POST);	/*MRI*/


	initStruct(mImagePre.slices, &mImagePre);
	
	initImage(192, 256, 256, &mImagePost, BACKGROUND);
	
	pixelSize = calloc(3, sizeof(double ));
	
	pixelSize[0] = 1.;
	pixelSize[1] = 1.;
	pixelSize[2] = 1.;
	setResolution(&mImagePre, pixelSize); /*Setting up pixel resolution*/
	
	pixelSize[0] = 1.;
	pixelSize[1] = 1.;
	pixelSize[2] = 1.;
	
	/*Setting up pixel resolution*/
	setResolution(&mImagePost, pixelSize);
	
	int i;
	
	loadImages(&mImagePre, filesPre, mImagePre.slices);
	mImagePre.cdmVolume[0] = mImagePre.dimx*0.5;
	mImagePre.cdmVolume[1] = mImagePre.dimy*0.5;
	mImagePre.cdmVolume[2] = mImagePre.slices*0.5;
	for (i = 0; i < 3; i++){
	  mImagePre.vertex[0][i] = 0;
	  mImagePre.vertex[1][i] = mImagePre.cdmVolume[i]*2;
	}
	
	printf("%d %d %d\n", mImagePre.dimx, mImagePre.dimy, mImagePre.slices);
	
	//initImage(mImagePre.slices,mImagePre.dimx, mImagePre.dimy, &mImagePost, BACKGROUND);
	
	
	
//	cor2sagit(&mImagePre);
//	sag2coron(&mImagePre);
//	cor2trans(&mImagePre);
	
	//printf("%d %d %d\n", new.);
}

int isEqual(double * var_1, double * var_2, int length){
  
  int i;
  
  for (i = 0; i < length; i++){
    if (var_1[i] != var_2[i]) return 0;
    
  }
  return 1;
}

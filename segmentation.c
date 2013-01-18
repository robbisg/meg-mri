#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "main.h"





/*This is the basic function for image filtering, the filter algorithm should have this simple signature*/
void filterImage (metaImage * image, void (* algorithm)(metaImage * image)){
#ifdef MEMORYDBG	    
	    printf("\nmemory before filtering!!!\n\n");
	    system("ps -aux | grep 'Mutual'");
#endif	    
	    
	    (*algorithm)(image);
#ifdef MEMORYDBG		    
	    printf("\nmemory after filtering!!!\n\n");
	    system("ps -aux | grep 'Mutual'");
#endif	   
}


void noFilter (metaImage * image){
    return;
}





void inverseThresholdInt (metaImage * image){
	  int ii;
	   for (ii = 0; ii < image->slices; ii++){
	  
	    image->image[ii] = inverseThreshold(image->image[ii], image->dimx, image->dimy, ii);
	  
	  }
}
p_size ** inverseThreshold (p_size ** image, int dimx, int dimy, int ii){

	int c = 0,i,j;
	p_size ** matrix_out = calloc (dimx, sizeof (p_size *));
   	for (c=0; c<dimx; c++){
		matrix_out[c] = calloc (dimy, sizeof (p_size));
		
	}
	
	for (i=0; i < dimx; i++) {
    		for (j=0; j < dimy; j++) {
                    
      			matrix_out[i][j] = (p_size)BACKGROUND;
                  
      		 }
        }
        
        for (i=0; i < dimx; i++) {
	    for (j=0; j < dimy; j++) {
                    
	      
		  if (BACKGROUND == 255 && image[i][j] >= THRESHOLD) matrix_out[i][j] = (p_size)IMAGE;
                  if (BACKGROUND == 0 && image[i][j] <= THRESHOLD) matrix_out[i][j] = (p_size)IMAGE;
	    }
        }
        
        
        
   	for (c=0; c<dimx; c++){
		free(image[c]);
		
	}
        free(image);
        
        return matrix_out;
}
void adaptiveThresholdInt (metaImage * image){
	  int ii;
	  
	  imageHistogram2D(image, image->histogram);
	  int thre = findThresholdII(image->histogram);
	  printf("threshold = %d\n", thre);
	  
	  thre = 26;
	   for (ii = 0; ii < image->slices; ii++){
	  
	    image->image[ii] = adaptiveThreshold(image->image[ii], image->dimx, image->dimy, ii, thre);
	  
	  }
}


/*Threshold segmentation: based on a threshold the images become binary*/
p_size ** adaptiveThreshold (p_size ** image, int dimx, int dimy, int ii, int threshold){

	int c = 0,i,j;
	p_size ** matrix_out = calloc (dimx, sizeof (p_size *));
   	for (c=0; c<dimx; c++){
		matrix_out[c] = calloc (dimy, sizeof (p_size));
		
	}
	
	for (i=0; i < dimx; i++) {
    		for (j=0; j < dimy; j++) {
                    
      			matrix_out[i][j] = (p_size)BACKGROUND;
                  
      		 }
        }
        
        for (i=0; i < dimx; i++) {
	    for (j=0; j < dimy; j++) {
                    
	      
		  if (BACKGROUND == 255 && image[i][j] < threshold) matrix_out[i][j] = (p_size)IMAGE;
                  if (BACKGROUND == 0 && image[i][j] > threshold) matrix_out[i][j] = (p_size)IMAGE;
	    }
        }
        
        for (c=0; c<dimx; c++){
	  free(image[c]);
		
	}
        free(image);
        
        
        return matrix_out;
}
void thresholdInt (metaImage * image){
	  int ii;
	  
	  imageHistogram2D(image, image->histogram);
	  int thre = 46;
	  
	   for (ii = 0; ii < image->slices; ii++){
	  
	    image->image[ii] = adaptiveThreshold(image->image[ii], image->dimx, image->dimy, ii, thre);
	  
	  }
}


/*Threshold segmentation: based on a threshold the images become binary*/
p_size ** threshold(p_size** image, int dimx, int dimy, int ii, int threshold)
{

	int c = 0,i,j;
	p_size ** matrix_out = calloc (dimx, sizeof (p_size *));
   	for (c=0; c<dimx; c++){
		matrix_out[c] = calloc (dimy, sizeof (p_size));
		
	}
	
	for (i=0; i < dimx; i++) {
    		for (j=0; j < dimy; j++) {
                    
      			matrix_out[i][j] =(p_size) BACKGROUND;
                  
      		 }
        }
        
        for (i=0; i < dimx; i++) {
	    for (j=0; j < dimy; j++) {
                    
	      
		  if (BACKGROUND == 255 && image[i][j] < (p_size)threshold) matrix_out[i][j] = (p_size)IMAGE;
                  if (BACKGROUND == 0 && image[i][j] > (p_size)threshold) matrix_out[i][j] = (p_size)IMAGE;
	    }
        }
        
       	for (c=0; c<dimx; c++){
	  free(image[c]);
		
	}
        free(image);
        
        
        return matrix_out;
}
/*This is bolla algorithm interface the algorithm pipeline should be implemented here*/
void bollaInt (metaImage * image){
      
	  int ii;
	  for (ii = 0; ii < image->slices; ii++){
	  
	    image->image[ii] = bolla(image->image[ii], image->dimx, image->dimy, ii);
	  
	  }
}

p_size ** bolla (p_size ** image, int dimx, int dimy, int ii){
	int i,j,c,partito=0;
	int * trovati_ora;
	c=0;
	p_size ** matrix_out = calloc (dimx, sizeof (p_size *));
   	for (c=0; c<dimx; c++){
		matrix_out[c] = calloc (dimy, sizeof (p_size));
		
	}
//	printf ("\nIn centerMass\n\n");
	if((trovati_ora =  calloc(MAX_trovati, sizeof(int))) == NULL){
	     fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
             exit(42);
	}
	
	

//	printf ("trovati ora = %d\n",trovati_ora[0]);
	for (i = 0; i< dimx; i++){
		for (j = 0; j<dimy; j++){	

			if(partito == 0){
			   if (image[i][j]>=THRESHOLD) {
				trovati_ora[0]=1;  
           			trovati_ora[trovati_ora[0]] = j*dimx+i;
            			bolla_connessa(i, j, &trovati_ora[0], image, dimx, dimy, matrix_out);
           			 if (trovati_ora[0] >  10000){
            			   partito = 1;
//             			   printf ("trovati ora = %d\n",trovati_ora[0]);
               				//matrix_out[i][j]=255;
           			 }
			    }		
			}
		}

	}
	for (c=0; c<dimx; c++){
		free(image[c]);
		
	}
        free(image);



	return matrix_out;
}

int  bolla_connessa(int y_vox, int z_vox, int *centro_ora, p_size ** image, int dimx, int dimy, p_size** matrix_out) {

int ii, kk, yy, xx,qq,i,j, vicinato=1;   //se vicinato=0 nella ENNESIMA corona non ho trovato piu' vicini
int bolla =1;                        // centro_ora[0]=numero centri trovati, il cui nome e' negli elementi seguenti
int voxel_already_bollati [dimx][dimy] ;

	for (i=0; i < dimx; i++) {
    		for (j=0; j < dimy; j++) {
      			voxel_already_bollati[i][j]=0;
      		 }
        }


  while(vicinato) {
    for (ii=1; ii <= *centro_ora  ; ii ++) {
      if (*centro_ora >= MAX_trovati) {
         printf ("stop MAX_trovati too small, min=%d\n",*centro_ora);
         exit(42);
      }
      vicinato=0;
      for (xx =  -1; xx <= 1; xx++) {
        for (yy =  -1; yy <= 1; yy++) {
           if ((xx!=0 || yy!=0) && ((*(centro_ora+ii))%dimy+yy>=0 && (*(centro_ora+ii))%dimy+yy<dimx
			&& (*(centro_ora+ii))/dimy+xx>=0 && (*(centro_ora+ii))/dimy+xx<dimy)) {
   //    printf("%d %d \n",(*(centro_ora+ii))%dimy+yy,(*(centro_ora+ii))/dimy+xx);
             if (image[(*(centro_ora+ii))%dimy+yy][(*(centro_ora+ii))/dimy+xx] >=  THRESHOLD) {
               if(voxel_already_bollati[(*(centro_ora+ii))%dimy+yy][(*(centro_ora+ii))/dimy+xx]==0) {
                 (*centro_ora) ++;
                 *(centro_ora + *centro_ora)=((*(centro_ora+ii))/dimy+xx)*dimy+(*(centro_ora+ii))%dimy+yy;
                 voxel_already_bollati[(*(centro_ora+ii))%dimy+yy][(*(centro_ora+ii))/dimy+xx]=bolla;
                 matrix_out[(*(centro_ora+ii))%dimy+yy][(*(centro_ora+ii))/dimy+xx]=255;
                 vicinato=1;
               }
             }
             else voxel_already_bollati[(*(centro_ora+ii))%dimy+yy][(*(centro_ora+ii))/dimy+xx]=-1;
           }
        }
      }
  //        printf ("centro_ora =%d\n",*(centro_ora));
  //        printf ("calcolato Y vox=%d and y vox=%d\n",(*(centro_ora+ii))%dimy,y_vox);
  //        printf ("calcolato Z vox=%d and z vox=%d\n",(*(centro_ora+ii))/dimy,z_vox);
    }
  }
}


/*This is passBand algorithm interface the algorithm pipeline should be implemented here*/
void passBandInt (metaImage * image){
  
      imageHistogram2D(image, image->histogram);
      int *threshold = findThreshold(image->histogram);
      
      int ii;
      
      for (ii = 0; ii < image->slices; ii++){
     
	  image->image[ii] = passBandThreshold(image->image[ii],image->dimx, image->dimy, threshold, ii);
      }
      

}

p_size ** passBandThreshold (p_size ** image, int dimx, int dimy, int * threshold, int ii){
	
	int c = 0,i,j;
	p_size ** matrix_out = calloc (dimx, sizeof (p_size *));
   	for (c=0; c<dimx; c++){
		matrix_out[c] = calloc (dimy, sizeof (p_size));
		
	}
	
	for (i=0; i < dimx; i++) {
    		for (j=0; j < dimy; j++) {
                    
      			matrix_out[i][j] = BACKGROUND;
                  
      		 }
        }
        
        
        double mean = ((double)threshold[1] - (double)threshold[0])/2.;
        
        int thresholdMin = ceil((double)threshold[1] - mean);
	int thresholdMax = ceil((double)threshold[1] + mean);
	
//	printf("thresholds = %d %d\n",thresholdMax,thresholdMin);
	
	
        
        for (i = 0; i< dimx; i++){
		for (j = 0; j<dimy; j++){	
	
			   if (image[i][j]>=thresholdMin && image[i][j]<=thresholdMax) {
				matrix_out[i][j] = IMAGE;
			}
		}

	}


	for (c=0; c<dimx; c++){
		free(image[c]);
		
	}
        free(image);

	return matrix_out;
	
	
}




/*This procedure compute the histogram of the image volume*/
/*The first variable is the input, and the algorithm fill in the second variable*/
void imageHistogram2D (metaImage * mImage, long int * histogram ){
      
      int i,j,k;
      int dimx = mImage->dimx;
      int dimy = mImage->dimy;
      int slices = mImage->slices;
      p_size *** image = mImage->image;
      
      int index;
    
      for (k = 0; k < slices; k++){
	for (i = 0; i < dimx; i++){
	  for (j = 0; j < dimy; j++){
	    
	    index = image[k][i][j];
	//    printf("%d ",index);
	    histogram[index]++ ;
	  }
	}
      }
}

/*DA VERIFICARE*/
int * findThreshold (long int * histogram){
    int i;
    
    long int max1 = 0, max2 = 0; 
    int index1 = 800, index2 = 800, flag;
    int * threshold = calloc(2, sizeof(int)); 
    
    long int gradient = 0.;
    
    max1 = histogram[0];
    index1 = 0;
  //  printf("[%d] %ld\n",0,histogram[0]);
    for (i = 1; i < 256; i++){
	//  printf("[%03d] %ld\n",i,histogram[i]);
	  
	  if (i < 254) gradient = histogram[i+2] - histogram[i];
	  else gradient = 0.;
	  
//	  printf(" max1 %d\n", max1);
	  
	  if (gradient < 0.) flag = 0; //La funzione decresce, 
	      else flag = 1;//La funzione cresce ---> massimo locale
	  
	  if (histogram[i] > max1) {//Massimo assoluto
	      max1 = histogram[i];
	      index1 = i;
	  }
	  else {
	    if (histogram[i] > max2 && flag == 1){//Massimo locale
		max2 = histogram[i];
		index2 = i;
	    }
	    
	  }
    }
    
    
    threshold[0] = index1;
    threshold[1] = index2;
    printf("histogram maxs = %ld %ld\n", max1, max2);
    printf("histogram max levels = %d %d\n", index1, index2);
    return threshold;

}

int findThresholdII(long int * histogram){
    int i;
    
    long int max1 = 0, max2 = 0; 
    int index1 = 800, index2 = 800, flag=0;
    int * threshold = calloc(2, sizeof(int)); 
    
    long int gradient = 0.;
    
    max1 = histogram[0];
    index1 = 0;
    
    for (i = 0; i < 256; i++){
    
       if (histogram[i] > max1) {//Massimo assoluto
	      max1 = histogram[i];
	      index1 = i;
      }
    }
    for (i = index1; i < 256; i++){
      if (i < 254) gradient = histogram[i+2] - histogram[i];
	  else gradient = 0.;
      
      if (gradient > 0.) flag++;
      else flag = 0;
      
     
      
      if (flag == 2) {index2 = i; break;}
      
    }
    
    printf("histogram maxs = %ld %ld\n", max1, max2);
    printf("histogram max levels = %d %d\n", index1, index2);
    return index2;
}

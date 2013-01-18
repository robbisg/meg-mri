#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "main.h"
#include "interpolation.h"




void weightedVector (double x, double y, double * weight);

unsigned long int countHist;
unsigned long int ** histogram;
//unsigned long int histogram[256][256];
static unsigned long int pre, post, miss;


double mutualInformation (metaImage * List, double * conf, int dim){
	
	
	//changeImage(conf, dim, slicesPost, 0);
	int j,c;
	
#ifdef DRAWHISTO	
	p_size ** normHistog = calloc(BIT, sizeof (p_size *));
	for (c=0; c<BIT; c++){
		normHistog[c] = calloc (BIT, sizeof (p_size));
	}
	
	;
#endif	
	
 	histogram = calloc (BIT, sizeof (unsigned long int *));
    	for (c=0; c<BIT; c++){
 		histogram[c] = calloc (BIT, sizeof (unsigned long int));
 	}
	
	metaImage imageFloat = List[0];//mImagePostSized;	/*imageFloat*/
	metaImage imageTarget = List[1];//mImagePreRot
	
	
		 
	pre = 0;
	post = 0;
	miss = 0;
 	inverseMappingNewLight(&imageFloat, &imageTarget, conf, metaBuildHistogram, nearestNeighborInterp);
	
	#ifdef ENTRDBG
	  printf("black pixels : pixPre = %d --- pixPost = %d\n",pre, post);
	#endif
	
	
	double entr; 
	
	entr = normalizedEntropy(histogram);
	//entr = entropy(histogram);
	int i, bit = BIT;
      #ifdef ENTRDBG
	
	printf("\n%ld  %ld  %ld \n",histogram[0][0],histogram[0][200], histogram[0][255]);
	printf("%ld  %ld  %ld \n",histogram[255][0],histogram[255][200],histogram[255][255]);
      
	
	
	double pixMatching = 0.0;
	double pixMatchingW = 0.0;
	double pixMatchingB = 0.0;
	
	double sum = histogram[0][0] + histogram[0][255] + histogram[255][0] + histogram[255][255];
	
	pixMatching = (double) (histogram[0][0] + histogram[255][255]) / sum ;
	
	printf ("\nPercentage of matching pixels = %.02lf\%\n", pixMatching*100);
	printf ("\nPercentage of matching pixels = %.02lf\%\n", 100.-pixMatching*100.);
	
	unsigned long int * pixelWhite = contaPixel(&imageTarget, &imageFloat);
	
	pixMatchingW = histogram[255][255] / (double)pixelWhite[0];
	pixMatchingB = histogram[0][0] / (double)pre;
	pixMatching = (double)pixelWhite[0]/ ((double) post + (double)pixelWhite[0]);
	
	printf ("\nPercentage of white pixel Image pre = %.02lf\%\n", pixMatchingW*100);
	printf ("\nPercentage of black pixel Image pre = %.02lf\%\n", pixMatchingB*100);
	printf ("\nPercentage of white in whole image = %.02lf\%\n", pixMatching*100);
	
	pixMatchingW = histogram[255][255] / (double)pixelWhite[1];
	pixMatchingB = histogram[0][0] /(double) post;
	pixMatching = (double)pixelWhite[1]/ ((double) post + (double)pixelWhite[1]);

	printf ("\nPercentage of white pixel Image post = %.02lf\%\n", pixMatchingW*100);
	printf ("\nPercentage of black pixel Image post = %.02lf\%\n", pixMatchingB*100);
	printf ("\nPercentage of white in whole image = %.02lf\%\n", pixMatching*100);
    #endif	
#ifdef DRAWHISTO	
	normalizeHistogram(histogram, normHistog, countHist, bit);

	for (i = 0; i<bit; i++){
		free(*(normHistog+i));
	}
	free(normHistog);
#endif

 	for (i = 0; i<bit; i++){
 		free(*(histogram+i));
 	}
	free(histogram);
		
	
	
	return entr;
} 


void buildHistogram ( int * coordPre, int * coordPost, p_size *** imagePre, p_size *** imagePost){

	int x, y, z;
	int i, j, k;
	int pixPre, pixPost;

	x = coordPost[X];
	y = coordPost[Y];
	z = coordPost[Z];	

	i = coordPre[X];
	j = coordPre[Y];
	k = coordPre[Z];


// 	if (z >= 0 && z < slicesPre && x >= 0 && x < dimPre[z][0] && y >= 0 && y < dimPre[z][1] )
// 				pixPre = imagePre[k][i][j];
// 				pixPre = readPixelMemory(z,x,y,PRE);
// 				
// 		else pixPre = 255;
				
		//pixPre = readPixel(slice,i,j,PRE);
		
		pixPost = imagePost[z][x][y];		

		pixPre = imagePre[k][i][j];

		if (pixPost == IMAGE) post++;
		if (pixPre == IMAGE) pre++;

		histogram[pixPre][pixPost]++;
		
}


void weightedVector (double x, double y, double * weight){		
	*(weight+3) = rint(x - floor(x)) * rint(y - floor(y));
	*(weight+2) = rint(x - floor(x)) * rint(ceil(y) - y);
	*(weight+1) = rint(ceil(x) - x) * rint(y - floor(y));
	*(weight+0) = rint(ceil(x) - x) * rint(ceil(y) - y);
	
	int i;
	

}

double entropy (unsigned long int ** ext_histogram){
	histogram = ext_histogram;
	countHist=0;
	double entr = 0;
    	int i,j;
	for (i = 0; i<BIT; i++){
		for (j=0; j<BIT;j++){
			countHist+=histogram[i][j];
		}
	}
	
	for (i = 0; i<BIT; i++){
		for (j = 0; j<BIT; j++){
		//	printf("%d    ",histogram[i][j]);
			entr+=probability(i,j)*logarithm(i,j);  
        	}
		//printf("\n");
    	}
	
    	return entr;
}


double probability (int f, int r){
       double prob = (double) histogram[f][r]/(double) countHist;
       return prob;
}


       
double logarithm(int f, int r){
       double probF = 0.0;
       double probR = 0.0;
       int i = 0;
       for (i = 0; i < BIT; i++){
           probF+=probability(f,i);
           probR+=probability(i,r);
       }
       double log2 = 0.0;
       if (probF == 0 || probR == 0 || probability (f,r) == 0) 
		log2 = 0;
       else 
		log2 = log(probability(f,r)/(probF*probR))/log(2);
      // printf ("log2 = %f\n",log2);
       return log2;
       }


void writeMatrix (int * coordOrig, int * coordRot, p_size *** imageOrig, p_size *** imageRot){
	
	/*La matrice output viene allocata (e riempita di bianco) in writeFinal, il compito di questa funzione
	è di riempire i buchi neri nella stessa*/
	
	/*Le coordinate passate alla funzione sono non resolution-influenced*/
	
	int x, y, z;
	int i, j, k;
	int pixPre, pixPost;

	i = coordRot[X];
	j = coordRot[Y];
	k = coordRot[Z];	

	x = coordOrig[X];
	y = coordOrig[Y];
	z = coordOrig[Z];

	imageRot[k][i][j] = imageOrig[z][x][y];
}

double normalizedEntropy (unsigned long int ** ext_histogram){
	histogram = ext_histogram;
	countHist=0;
	double entr = 0;
    	int i,j;
	for (i = 0; i<BIT; i++){
		for (j=0; j<BIT;j++){
			countHist+=histogram[i][j];
		}
	}
	  
	//  printf("count = %ld\n", countHist);
	
	unsigned long int probPost, probPre; double probTot = 0.0;
	double entroPost = 0.0; double entroPre = 0.0;
	
	double entroTot = 0.0;
	for (i = 0; i < BIT; i++){
	    probPost = 0;
	    probPre = 0;
	    for (j = 0 ; j < BIT; j++){
		probPost += histogram[i][j];
		probPre += histogram[j][i];
		probTot = probability(i,j);
		if (probTot == 0)  entroTot+=0;
		else entroTot += probTot * log(probTot);
	    }
	    
	    if (probPost == 0)  entroPost+=0;
	    else entroPost += (double)probPost/(double)countHist * log((double)probPost/(double)countHist);
	    if (probPre == 0)  entroPre+=0;
	    else entroPre += (double)probPre/(double)countHist * log((double)probPre/(double)countHist);
	    if (probTot == 0)  entroTot+=0;
	    else entroTot += probTot * log(probTot);
	    
	    
//	    printf("probPost = %ld probPre = %ld\n", probPost,probPre);
	    
	}
	
	//  printf("H(post)=%lf H(pre)=%lf H(pre,post)=%lf\n",entroPost,entroPre,entroTot);
	    
	//  entr = (-entroPost - entroPre) + entroTot;
	  
	//  printf("classical = %lf\n", entr);
	  
	  entr = -(entroPost + entroPre)/-entroTot;
	  
	//  printf("normalized = %lf\n", entr);
	
  	return entr;

}
void metaBuildHistogram (double * coordOrig, double * coordRot, metaImage * imageOrig, metaImage * imageRot,
			 int (* interpolation)(double * coord, metaImage * image)){
	
	int x, y, z;
	int i, j, k;
	int pixPre, pixPost;

	x = rint(coordRot[X]);
	y = rint(coordRot[Y]);
	z = rint(coordRot[Z]);	

	i = rint(coordOrig[X]);
	j = rint(coordOrig[Y]);
	k = rint(coordOrig[Z]);
	
	int flag = 0;

	if (k >= 0 && k < imageOrig->slices && i >= 0 && i < imageOrig->dimx
		&& j >= 0 && j < imageOrig->dimy){
		//pixPre = imageOrig->image[k][i][j];
		pixPre = (*interpolation)(coordOrig, imageOrig);
		
	}
	else {
	  flag = 1;
	}
		  
	if(z >= 0 && z < imageRot->slices && x >= 0 && x < imageRot->dimx
		&& y >= 0 && y < imageRot->dimy){
	
		//pixPost = imageRot->image[z][x][y];		
		pixPost = (*interpolation)(coordRot, imageRot);
	}
	else {
	  if (flag == 1) {flag = 3;}
	  else {flag = 2;}
	}
	
	
	
	if (pixPost == BACKGROUND) post++;
	if (pixPre  == BACKGROUND) pre++;
	
	
	
 	 if (flag == 0) histogram[pixPost][pixPre]++;
	 if (flag == 1) histogram[pixPost][BACKGROUND]++;
	 if (flag == 2) histogram[BACKGROUND][pixPre]++;
	 //else histogram[0][255]++;
	
}

void metaWriteMatrix (double * coordOrig, double * coordRot, metaImage * imageOrig, metaImage * imageRot,
		       int (* interpolation)(double * coord, metaImage * image)){
	double x, y, z;
	int i, j, k;
	int pixPre, pixPost;

	i = rint(coordRot[X]);
	j = rint(coordRot[Y]);
	k = rint(coordRot[Z]);	

	x = coordOrig[X];
	y = coordOrig[Y];
	z = coordOrig[Z];
	
	int value = (*interpolation)(coordOrig, imageOrig);
	
	if (k >= 0 && k < imageRot->slices && i >= 0 && i < imageRot->dimx
		&& j >= 0 && j < imageRot->dimy/* && 
	    z >= 0 && z < imageOrig->slices && x >= 0 && x < imageOrig->dimx
		&& y >= 0 && y < imageOrig->dimy*/){
	    imageRot->image[k][i][j] = value;
	}
}

void normalizeHistogram(unsigned long int ** histogram, p_size ** output, unsigned long int count, int bit){
    static int flag = 0;
    int i,j;
    unsigned long int max;
    double value;
    
    for (i = 0; i < bit; i++){
      for (j = 0; j < bit; j++ ){
	if (histogram[i][j]>max) max = histogram[i][j];
      
      }
    }
    
    
    for (i = 0; i < bit; i++){
      for (j = 0; j < bit; j++ ){
	  value = bit * ((double) histogram[i][j])/((double) max*0.5);
	  output[i][j] = rint(value);
      }
    }
    
    flag++;
    
    char file_str[150]="\0";
   
    sprintf(file_str,"result/histogram_%03d.txt",flag);
    FILE * file;
    file = fopen(file_str, "w");
    
   
    for (i = 0; i < bit; i++){
      for (j = 0; j < bit; j++ ){
	   
	 fprintf(file, "%ld ",histogram[i][j]);
      }
      fprintf(file,"\n");
    }
   
 
  
    fclose(file);
    
    
    writeTiffcustom(output, "result/histogram.tif", bit, bit);
    
}






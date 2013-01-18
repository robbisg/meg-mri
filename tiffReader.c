#include <tiffio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "main.h"

TIFF *in, *out;
uint32 imagelength, imagewidth, strip_byte_count;
uint16 bitSample, samplePixel, compress, photo, fill_order, plan_conf, color1, color2, color3, color4;
int rows_strip, ns, bufsize;
unsigned long imageOffset, result;

p_size ** readTiff(int * dim, char * file){
	
	int ii; 
	float xres, yres;
	
	
	char out_str[20]="\0";

	
	
	//Apro il primo file
	if((in = TIFFOpen(file, "r")) == NULL){
    		fprintf(stderr, "Could no open incoming image\n");
    		exit(42);
  	}

	bufsize  = TIFFStripSize(in);
	ns = TIFFNumberOfStrips(in);

	if (TIFFGetField(in, TIFFTAG_BITSPERSAMPLE , &bitSample)) ;      
	if (TIFFGetField(in, TIFFTAG_ROWSPERSTRIP, &rows_strip)) ;
	if (TIFFGetField(in, TIFFTAG_IMAGELENGTH, &imagelength)) ;
	if (TIFFGetField(in, TIFFTAG_IMAGEWIDTH , &imagewidth)) ;
	if (TIFFGetField(in, TIFFTAG_COMPRESSION, &compress)) ;
	if (TIFFGetField(in, TIFFTAG_FILLORDER, &fill_order)) ;
	if (TIFFSetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplePixel)) ;
//	printf("rows per strip = %d compression= %d fill_order=%d samplePixel = %d\n",rows_strip, compress, fill_order, samplePixel );
//	printf("imagelength = %d image width = %d bufsize = %d ns = %d\n",imagelength, imagewidth,bufsize,ns);
	char *buf;
        int  kk,s;
        p_size tempbyte;
	
	
	/*Costruisco la matrice che conterrà l'immagine*/
	p_size ** image = calloc (imagelength, sizeof (p_size *));
	int c;
   	for (c=0; c<imagelength; c++){
		image[c] = calloc (imagewidth, sizeof (p_size));
	}
	
	dim[0] = imagelength;
	dim[1] = imagewidth;
	
	/*tsize_t bufsize  = TIFFStripSize(in);
        tstrip_t s, ns = TIFFNumberOfStrips(in);*/

        uint32 *bytecounts;
        
	//printf("%x\n",output);
        //uint32 imagelength,imagewidth;

	if((buf = malloc(bufsize*ns)) == NULL){
	        fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
	        exit(42);
        }
	int i;
        imageOffset = 0;
	for (s = 0; s < ns ; s++){
            if((result = TIFFReadEncodedStrip (in, s,buf+ imageOffset ,bufsize)) == -1){
               fprintf(stderr, "Read error on input strip number %d\n", s);
                       exit(42);
            }

	       /* Legge bufsize caricato in memoria bit per bit  :-) */
	
	 
            for (kk=0; kk< bufsize/imagewidth; kk++) {

               for(ii = kk*imagewidth+imageOffset; ii < (kk+1)*imagewidth+imageOffset; ii++){
                  tempbyte = ((unsigned) buf[ii]) % 256;
		  image[ii/imagewidth][ii%imagewidth] = tempbyte;
		 //  if (tempbyte > 150) printf("t %d b %d i %d\n",buf[ii],tempbyte,image[ii/imagewidth][ii%imagewidth]);
		//  printf("%d %d %d %d\n",ii, imagewidth,ii/imagewidth,ii%imagewidth);
                 }
        
            }
	//	printf("\n");	
	    imageOffset += result;
        }
	
	free(buf);
	TIFFClose(in);
	return image;

}

int one_byte_char_into_integer(char *byte_char, int array_bin) {

   int count, tempbyte;
   unsigned c,displayMask =1 << 15;
   
   count=0;
   tempbyte=0;
   for (c = 1; c <= 8; c++) {
      count++;
      if ( (unsigned) byte_char[array_bin] & displayMask ) {
         tempbyte = tempbyte+pow(2,(8-count));
      }
      byte_char [array_bin] <<= 1;
   }

  /* printf("%4d\n", tempbyte); */

             /*   printf("byte_char=%2x tempbyte=%d \n",
                      (p_size) byte_char, tempbyte); */
   return (tempbyte);
}

void writeTiff (p_size ** matrix_out, char * name){

	char * buf_out;    

     if((out = TIFFOpen(name, "w")) == NULL){
       fprintf(stderr, "Could not open output.tif for writinge\n");
       exit(42);
     } 
       
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, rows_strip);
  TIFFSetField(out,TIFFTAG_IMAGELENGTH, imagelength);
  TIFFSetField(out, TIFFTAG_IMAGEWIDTH , imagewidth);
  TIFFSetField(out,TIFFTAG_BITSPERSAMPLE , bitSample);

  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, 1);

  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(out, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
//  TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

  TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);


  TIFFSetField(out, TIFFTAG_COMPRESSION   , compress   );
//  TIFFSetField(out,TIFFTAG_FILLORDER     , 1  );
  
	

   if((buf_out = (char *) calloc(bufsize,ns)) == NULL){
      fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
      exit(42);
    }
	int ii,s,kk;
    imageOffset = 0;

    for (s = 0; s < ns ; s++){
           for(ii = imageOffset; ii < imagewidth+imageOffset; ii++){
              buf_out[ii] = matrix_out[s][ii-imageOffset];
		//if (buf_out[ii] == -1) buf_out[ii] = 255;
		//printf("%d %d |", ii-imageOffset, ns-s-1 );
		
           }
           if((TIFFWriteEncodedStrip( out, s, buf_out+ imageOffset ,bufsize))==-1){
                fprintf(stderr, "Write error on input strip number %d\n", s);
                 exit(42);
           }
	
      imageOffset = imageOffset + bufsize;
    }
	TIFFClose(out);
    _TIFFfree(buf_out);
     
}


int readPixelMemory (int slice, int x, int y, int flag){

	if (flag == PRE) {
		return imagePre[slice][x][y];
	}
	else{
		return imagePost[slice][x][y];
	}
}



int readPixel (int slice, int x, int y, int flag){
	
	int ii; 
	FILE *meshFilePtr;
	float xres, yres;
	char * file = calloc(100,sizeof(char));
	
	char out_str[20]="\0";

	if (flag == 0) {
		sprintf(file,"pre/fill-slice_1%02d.tif",slice);
	}
	else{
		sprintf(file,"post/fill-slice_1%02d.tif",slice);
	}
	
//	printf("\n\n\n%s\n\n\n",argv);


	//Apro il primo file
	if((in = TIFFOpen(file, "r")) == NULL){
    		fprintf(stderr, "Could not open incoming image\n");
    		exit(42);
  	}

	bufsize  = TIFFStripSize(in);
	ns = TIFFNumberOfStrips(in);
	int fettine=1;
//	printf ("slices=%d\n",slices); 
//	printf ("bufsize=%d,numOfstrips=%d\n",bufsize,ns);
	if (TIFFGetField(in, TIFFTAG_BITSPERSAMPLE , &bitSample)) ;
//	printf ("sample per pixel =%d\n",samplePixel);
       
	if (TIFFGetField(in, TIFFTAG_ROWSPERSTRIP, &rows_strip)) ;
	if (TIFFGetField(in, TIFFTAG_IMAGELENGTH, &imagelength)) ;
	if (TIFFGetField(in, TIFFTAG_IMAGEWIDTH , &imagewidth)) ;
	if (TIFFGetField(in, TIFFTAG_COMPRESSION, &compress)) ;
	if (TIFFGetField(in, TIFFTAG_FILLORDER, &fill_order)) ;
	
//	printf("rows per strip = %d compression= %d fill_order=%d\n",rows_strip, compress, fill_order );
//	printf("imagelength = %d image width = %d bufsize = %d\n",imagelength, imagewidth);

	 

	char *buf, *buf_out;
        int tempbyte, kk,s;
	
	
	/*tsize_t bufsize  = TIFFStripSize(in);
        tstrip_t s, ns = TIFFNumberOfStrips(in);*/

        uint32 *bytecounts;
        

        //uint32 imagelength,imagewidth;

	if((buf = malloc(bufsize*sizeof(char))) == NULL){
	        fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
	        exit(42);
        }
	int i;
        imageOffset = 0;
	
	s = x/rows_strip;

        if((result = TIFFReadEncodedStrip (in, s,buf+ imageOffset ,bufsize)) == -1){
           fprintf(stderr, "Read error on input strip number %d\n", s);
                exit(42);
        }

	ii = (x%rows_strip*imagewidth)+y;

       /* Legge bufsize caricato in memoria bit per bit  :-) */
	tempbyte=one_byte_char_into_integer(buf, ii);
	//printf("buf/imwidth = %d ",bufsize/imagewidth);
	//printf("\noffset = %d ",imageOffset);
		
	   
       // printf("ii %d s %d %d", ii, s,tempbyte); 
            

		
	free(file);
	free(buf);
	TIFFClose(in);
	return tempbyte;
	
}

void writeMetaImage (metaImage * image, char * path){

    int i, slices = image->slices;
    
    int charlength = strlen(path);
    
    char * argF = calloc(charlength+20, sizeof(char));
    
     int perc;
      printf("\nWriting images to disk ...   \n");
	   for (i = 0; i < slices ; i++){
	
	     perc = ceil((double)i/(double)slices)*100;
	     advanceProgressPercentage(perc);
	      sprintf(argF, "%s-%03d.tif", path, i);
	      writeTiffcustom(image->image[i], argF, image->dimx, image->dimy);
	  }

}

void writeTiffcustom (p_size ** matrix_out, char * name, int dimx, int dimy){

	char * buf_out;    

    if((out = TIFFOpen(name, "w")) == NULL){
       fprintf(stderr, "Could not open output.tif for writing\n");
       exit(42);
     } 
     
     ns = dimx;
     bufsize = dimy;
      
  
  TIFFSetField(out,TIFFTAG_IMAGELENGTH, dimx);
  TIFFSetField(out, TIFFTAG_IMAGEWIDTH , dimy);
  TIFFSetField(out,TIFFTAG_BITSPERSAMPLE , 8);

  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, 1);

  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(out, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
  TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

  TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);


  TIFFSetField(out, TIFFTAG_COMPRESSION   , compress   );

  
//	printf("dimx = %d dimy = %d\n",dimx,dimy);
	
	
//	printf("rows per strip = %d compression= %d fill_order=%d\n",rows_strip, compress, fill_order );
//	printf("imagelength = %d image width = %d bufsize = %d ns = %d\n",dimx, dimy, bufsize,ns);
	

   if((buf_out = calloc(bufsize*ns, sizeof(char))) == NULL){
      fprintf(stderr, "Could not allocate enough memory for the uncompressed image\n");
      exit(42);
    }
	
	
	int ii,s,kk;
    imageOffset = 0;  
    
    
    for (s = 0; s < ns ; s++){
           for(ii = imageOffset; ii < bufsize+imageOffset; ii++){
              buf_out[ii] = matrix_out[s][ii - imageOffset];
		//if (buf_out[ii] == -1) buf_out[ii] = 255;
	      
	  //    if (matrix_out[s][ii - imageOffset] > 100) printf("b %d m %d\n", buf_out[ii], matrix_out[s][ii - imageOffset]);
				
           }
           if((TIFFWriteEncodedStrip( out, s, buf_out+ imageOffset ,bufsize))==-1){
                fprintf(stderr, "Write error on input strip number %d\n", s);
                 exit(42);
           }
   //        printf("\n");
	
      imageOffset = imageOffset + bufsize;
    }
    
    
 //   printf("\n\n\n\n\n\n\nArrivo qui\n\n\n %d %d %d\n\n\n\n",s,  ii, imageOffset);
    
    TIFFClose(out);
    _TIFFfree(buf_out);
     
}

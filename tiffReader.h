/*functions implemented in tiffReader.c*/
p_size ** readTiff(int * dim, char * file);
void writeTiff (p_size ** matrix_out, char * name);
int readPixelMemory (int slice, int x, int y, int flag);
void writeTiffcustom (p_size ** matrix_out, char * name, int dimx, int dimy);
void writeMetaImage (metaImage * image, char * name);
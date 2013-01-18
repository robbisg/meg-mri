Mutual.x : main.c asa.o interpolation.o transformation.o pso.o FF_asa.o segmentation.o entropy.o writingUtils.o tiffReader.o fileReader.o preProcessing.o metaImage.o
	gcc -g -Wall main.c interpolation.o asa.o FF_asa.o pso.o segmentation.o writingUtils.o metaImage.o transformation.o entropy.o tiffReader.o fileReader.o preProcessing.o -lm -ltiff -o Mutual.x
asa.o : asa.c
	gcc -c -g -lm asa.c
entropy.o : entropy.c
	gcc -c -g -lm entropy.c
tiffReader.o : tiffReader.c
	gcc -c -lm -ltiff -g tiffReader.c
fileReader.o : fileReader.c
	gcc -c -lm -g fileReader.c
transformation.o : transformation.c
	gcc -c -g -lm transformation.c 
segmentation.o : segmentation.c
	gcc -c -lm -g segmentation.c
preProcessing.o : preProcessing.c
	gcc -c -lm -g preProcessing.c
writingUtils.o : writingUtils.c
	gcc -c -lm -g writingUtils.c
metaImage.o : metaImage.c
	gcc -c -lm -g metaImage.c
pso.o : pso.c
	gcc -c -lm -g pso.c
interpolation.o : interpolation.c
	gcc -c -lm -g interpolation.c
FF_asa.o : FF_asa.c
	gcc -c -lm -g FF_asa.c
install :
	mkdir -p pre/
	mkdir -p post/
	mkdir -p result/diff
	mkdir -p result/fulfilled
	mkdir -p result/sub_sampled
clean : 
	rm -f *.o *~
	rm -f Mutual.x
	rm -f pre/*.tif
	rm -f post/*.tif
	rm -f result/diff/*.tif
	rm -f result/*.tif
	rm -f result/fulfilled/*.tif
	rm -f result/sub_sampled/*.tif
run_script: 
	./Mutual.x
copy :
	mkdir -p result_
	mv result/ result_
	mv lista_param_win.par result_/
	
# .SUFFIXES: .o .
# CC = gcc
# RM = rm -f
# PROG = Mutual.x
# INCL = main.h
# LIB = -ltiff -lm
# OPTS = -g
# WALL = -Wall
# SOURCE = main.c \
# 	metaImage.c \
# 	preProcessing.c \
# 	segmentation.c \
# 	transformation.c \
# 	tiffReader.c \
# 	asa.c \
# 	fileReader.c \
# 	entropy.c \
# $(OBJ) = $(SOURCE:.c = .o)
# %.o : %.c
# 	$(CC) -c $(WALL) $(LIB) $< -o $@
# $(PROG) : $(OBJ)
# 	$(CC) $(LIB) $(OPTS) $(WALL) $(OBJ) -o $(PROG)
# clean : 
# 	$(RM) *.o *~
# 	$(RM) $(PROG)

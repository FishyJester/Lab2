#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <string.h>

#define twoPI 6.28318512
#define MAX 10000000000


long d, cores, dim;
int **attractors, **convergences;
char *row_done;
pthread_mutex_t row_done_mutex;
double stepSize;
struct timespec ts;
char convPixCols[612];
char *convColors[51], *attrColors[11];
char attrFilename[40], convFilename[40]; 

void calcTrueRoots(double trueRe[], double trueIm[], const double dInv);
void calcNextZ(double *re, double *im, const double dInv);
void rotationFunc(const double *re, const double *im, double *reStar, double *imstar);
void *newton(void *compArg);
void *writeFunction(void *writeArg);


int main(int argc, char* argv[]){
  char *inp = argv[1] + 1;
  char *ptr;

  d = strtol(argv[3], &ptr, 10);

  if(*inp == 't'){
    cores = strtol(argv[1] + 2, &ptr, 10);
    dim = strtol(argv[2] + 2, &ptr, 10);
  }
  else{
    cores = strtol(argv[2] + 2, &ptr, 10);
    dim = strtol(argv[1] + 2, &ptr, 10);
  }
  
  sprintf(attrFilename, "newton_attractors_x%d.ppm", d);
  sprintf(convFilename, "newton_convergence_x%d.ppm", d);
  
  pthread_t threads[cores+1];
  
  ts.tv_sec = 0;
  if(d == 4)
  	ts.tv_nsec = 500000000;
  else
  	ts.tv_nsec = 100000000;

  convColors[0] = convPixCols + 0;
  int lastByte = 0;
  for(size_t i = 1; i < 51; ++i){
	int numBytes = log10(i*5)+1;	
	lastByte += (numBytes*3+4);
    	convColors[i] = convPixCols + lastByte;
    	sprintf(convColors[i], "%d %d %d ", i*5,i*5,i*5);
  }
  convColors[0] = "0 0 0 ";
  convColors[50] = "255 255 255 ";
  
  attrColors[0] = "0 0 0 ";
  attrColors[1] = "0 128 0 ";
  attrColors[2] = "128 0 0 ";
  attrColors[3] = "128 128 0 ";
  attrColors[4] = "128 0 128 ";
  attrColors[5] = "128 128 128 ";
  attrColors[6] = "0 128 128 ";
  attrColors[7] = "0 0 128 ";
  attrColors[8] = "100 149 237 ";
  attrColors[9] = "255 255 224 ";
  attrColors[10] = "255 215 0 ";

  

  stepSize = 4.0/(dim-1);
  attractors = malloc(sizeof *attractors * dim);
  convergences = malloc(sizeof *convergences * dim);
  row_done = calloc(dim, sizeof *row_done);

  if(d == 1){
	FILE *attrFile = fopen(attrFilename, "wb");
	FILE *convFile = fopen(convFilename, "wb");

	for(size_t i = 0; i < dim; ++i){
		for(size_t j = 0; j < dim; ++j){
			fwrite(attrColors[1], 1, 6, attrFile);
			fwrite(convColors[0], 1, 6, convFile);
		}

		fwrite("\n", 1, 1, convFile);
		fwrite("\n", 1, 1, attrFile);
	}
	fclose(attrFile);
	fclose(convFile);
  }else{


  	for(size_t i = 0; i < cores; ++i){
    		size_t *compArg = malloc(sizeof *compArg); 
    		*compArg = i;
    
    		if(pthread_create(&threads[i], NULL, newton, (void*)compArg)){
      			printf("error creating thread");
      			exit(1);
    		}
  	}
  
  	if(pthread_create(&threads[cores], NULL, writeFunction, NULL)){
    		printf("error creating thread");
    		exit(1);
    	}

  	for(size_t i = 0; i < cores; ++i)
    		pthread_join(threads[i], NULL);

  	pthread_join(threads[cores], NULL);
  }
  free(attractors);
  free(convergences);
  free(row_done);
  
  return 0;
  
}


void *newton(void *compArg){
  size_t row_index = *((size_t*)compArg);
  free(compArg);
  double re0, im0, dInv, trueRe[d], trueIm[d];
  dInv = 1.0/d;
  //printf("%f\n", dInv);
  calcTrueRoots(trueRe, trueIm, dInv); 
  //printf("%f %f %f %f %f %f\n", trueRe[0], trueIm[0], trueRe[1], trueIm[1], trueRe[2], trueIm[2]);
  
  for(row_index; row_index < dim; row_index+=cores){
    int *attractor = malloc(sizeof *attractor * dim);
    int *convergence = malloc(sizeof *convergence * dim);
    re0 = -2.0;
    im0 = 2.0 - (row_index * stepSize);
    
    for(size_t j = 0; j < dim; ++j){
      double re = re0, im = im0;
      size_t iteration = 0;
      int exitFlag = 0;

// Handles one coordinate at a time
      while(iteration < 1000){
	
	calcNextZ(&re, &im, dInv);
	
	if(re <= 1e-3 && re >= -1e-3 && im <= 1e-3 && im >= -1e-3){
	  attractor[j] = 0;
	  convergence[j] = iteration;
	  break;
	}
	
	if(re >= MAX || re <= -MAX || im >= MAX || im <= -MAX){
	  attractor[j] = 0;
	  convergence[j] = iteration;
	  break;
	}
	
	for(size_t i = 0; i < d; ++i){
	  double reDist = re - trueRe[i];
	  double imDist = im - trueIm[i];

	  if(reDist <= 1e-3 && reDist >= -1e-3 && imDist <= 1e-3 && imDist >= -1e-3){
	    attractor[j] = i+1;
	    convergence[j] = iteration;
	    exitFlag = 1;
	    break;
	  }
	}
	
	if(exitFlag == 1)
	  break;

	iteration++;
      }
      
      if(iteration > 50)
	      convergence[j] = 50;

      re0 += stepSize;
    }
    attractors[row_index] = attractor;
    convergences[row_index] = convergence;
    
    pthread_mutex_lock(&row_done_mutex);
    row_done[row_index] = 1;
    pthread_mutex_unlock(&row_done_mutex);
  }
  return NULL;
}

void *writeFunction(void *writeArg){
  char *row_done_loc = calloc(dim, sizeof *row_done_loc);
  char attrHeaderBuff[50], convHeaderBuff[50];
  FILE *attrFile, *convFile;
  int maxRGB = 255;	

  attrFile = fopen(attrFilename, "wb");
  convFile = fopen(convFilename, "wb");
  sprintf(convHeaderBuff, "P3\n%ld %ld\n%d\n", dim, dim, maxRGB);
  fwrite(&convHeaderBuff, 1, strlen(convHeaderBuff), convFile);
  sprintf(attrHeaderBuff, "P3\n%ld %ld\n255\n", dim, dim);
  fwrite(&attrHeaderBuff, 1, strlen(attrHeaderBuff), attrFile);
  
  //printf("Headers written\n");

  for(size_t k = 0; k < dim;){
    pthread_mutex_lock(&row_done_mutex);
    if(row_done[k] != 0){
      memcpy(row_done_loc, row_done, dim);
    }
    pthread_mutex_unlock(&row_done_mutex);

    if(row_done_loc[k] == 0){
	nanosleep(&ts, NULL);
      	continue;
    }
    
    for(; k < dim && row_done_loc[k] != 0; ++k){
      int *attractor_loc, *convergence_loc;
      attractor_loc = attractors[k];
      convergence_loc = convergences[k];

      for(size_t n = 0; n < dim; ++n){
	fwrite(attrColors[attractor_loc[n]], 1, strlen(attrColors[attractor_loc[n]]), attrFile);
	fwrite(convColors[convergence_loc[n]], 1, strlen(convColors[convergence_loc[n]]), convFile);
      }
      fwrite("\n", 1, 1, attrFile);
      fwrite("\n", 1, 1, convFile);
      free(attractor_loc);
      free(convergence_loc);
    }
  }
  fclose(attrFile);
  fclose(convFile);
  free(row_done_loc);
  
  return NULL;
}


void calcTrueRoots(double trueRe[], double trueIm[], const double dInv){
  double arg, reRot, imRot;
  arg = twoPI*dInv;
  reRot = cos(arg);
  imRot = sin(arg);

  trueRe[0] = 1.0;
  trueIm[0] = 0.0;

  for(int n = 1; n < d; ++n){
    trueRe[n] = trueRe[n-1] * reRot - trueIm[n-1] * imRot;
    trueIm[n] = trueRe[n-1] * imRot + trueIm[n-1] * reRot;
  }
}

void calcNextZ(double *re, double *im, const double dInv){

  double reStar, imStar;
  rotationFunc(re, im, &reStar, &imStar);

  double rInv = 1.0/(reStar*reStar + imStar*imStar);
  //printf("rInv is: %f\n", rInv);
  *re = *re - dInv * *re + dInv * reStar * rInv;
  //printf("re in Z: %f\n", *re);
  *im = *im - dInv * *im - dInv * imStar * rInv;
  //printf("im in Z: %f\n", *im);
}

void rotationFunc(const double *re, const double *im, double *reStar, double *imStar){

*reStar = *re;
*imStar = *im;

  if(d != 2){
    double tmp;

    for(time_t i = 0; i < d-2; ++i){
      tmp = *re * *reStar - *im * *imStar;
      *imStar = *im * *reStar + *re * *imStar;
      *reStar = tmp;
    }
  }
}

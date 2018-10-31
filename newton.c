#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <string.h>

#define twoPI 6.28318512
#define MAX 10000000000


long d, cores, dim;
int **attractors;
char *row_done;
pthread_mutex_t row_done_mutex;
double stepSize;
struct timespec ts;
char pixCols[44];
char *colors[11];

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
  
  
  pthread_t threads[cores+1];
  ts.tv_sec = 0;
  ts.tv_nsec = 100000;

  for(size_t i = 0, j = 0; i < 11; ++i, j+=4){
    colors[i] = pixCols + j;
    sprintf(colors[i], "%d %d %d ", i,i,i);
  }

  stepSize = 4.0/(dim-1);
  attractors = malloc(sizeof *attractors * dim);
  /* convergences = malloc(sizeof *convergences * dim);*/
  row_done = calloc(dim, sizeof *row_done);
  
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

  
  free(attractors);
  free(row_done);
  
  return 0;
  
}


void *newton(void *compArg){
  size_t row_index = *((size_t*)compArg);
  free(compArg);
  double re0, im0, dInv, trueRe[d], trueIm[d];
  dInv = 1.0/d;

  calcTrueRoots(trueRe, trueIm, dInv);
  
  for(row_index; row_index < dim; row_index+=cores){
    int *attractor = malloc(sizeof *attractor * dim);
    re0 = -2.0;
    im0 = 2.0 - (row_index * stepSize);
    //printf("im0 is: %f\n", im0);

    for(size_t j = 0; j < dim; ++j){
      double re = re0, im = im0;
      size_t iteration = 0;
      int exitFlag = 0;
    
      while(iteration < 1000){
	calcNextZ(&re, &im, dInv);
	
	if(re <= 1e-3 && re >= -1e-3 && im <= 1e-3 && im >= -1e-3){
	  attractor[j] = 0;
	  break;
	}
	
	if(re >= MAX || re <= -MAX || im >= MAX || im <= -MAX){
	  attractor[j] = 0;
	  break;
	}

	for(size_t i = 0; i < d; ++i){
	  double reDist = re - trueRe[i];
	  double imDist = im - trueIm[i];

	  if(reDist <= 1e-3 && reDist >= -1e-3 && imDist <= 1e-3 && imDist >= -1e-3){
	    attractor[j] = i+1;
	    //printf("reDist: %f imDist %f index: %d\n", reDist, imDist, i+1);
	    exitFlag = 1;
	    break;
	  }
	}
	
	if(exitFlag == 1)
	  break;

	iteration++;
      }
      re0 += stepSize;
    }
    attractors[row_index] = attractor;

    printf("row_index: %d\n", row_index);
    
    pthread_mutex_lock(&row_done_mutex);
    row_done[row_index] = 1;
    pthread_mutex_unlock(&row_done_mutex);
  }
  printf("%f %f %f %f %f %f\n", trueRe[0], trueIm[0], trueRe[1], trueIm[1], trueRe[2], trueIm[2]);
  return NULL;
}

void *writeFunction(void *writeArg){
  char *row_done_loc = calloc(dim, sizeof(char));
  char headerBuff[40];
  FILE *fp;
  fp = fopen("testfile.ppm", "wb");
  sprintf(headerBuff, "P3\n%ld %ld\n%ld\n", dim, dim, d);
  fwrite(&headerBuff, 1, strlen(headerBuff)+1, fp);
  

  for(size_t k = 0; k < dim;){
    pthread_mutex_lock(&row_done_mutex);
    if(row_done[k] != 0){
      memcpy(row_done_loc, row_done, sizeof(char) * dim);
    }
    pthread_mutex_unlock(&row_done_mutex);

    if(row_done_loc[k] == 0){
      continue;
    }
    
    for(; k < dim && row_done_loc[k] != 0; ++k){
      int *attractor_loc;
      attractor_loc = attractors[k];

      for(size_t n = 0; n < dim; ++n){
	fwrite(colors[attractor_loc[n]], 1, 4, fp);
      }
      fwrite("\n", 1, 1, fp);
      free(attractor_loc);
    }
  }
  fclose(fp);
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

  *re = *re - dInv * *re + dInv * reStar * rInv;
  *im = *im - dInv * *im - dInv * imStar * rInv;
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

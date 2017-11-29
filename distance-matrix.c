#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define BUFFER_LENGTH 40

void readFile(char *inputfile,
              unsigned int **xs,
              unsigned int **ys,
              unsigned int *locationCount) {
  char buffer[BUFFER_LENGTH];
  char *end;

  FILE *fh = fopen(inputfile, "r");
  if (fh == NULL) {
    perror("readFile");
    exit(EXIT_FAILURE);
  }
  printf("readFile: Opened %s\n", inputfile);
  fgets(buffer, BUFFER_LENGTH, fh);
  (*locationCount) = (unsigned int) strtol(buffer, &end, 10);
  printf("readFile: Allocating arrays...\n");
  (*xs) = (unsigned int *) malloc((*locationCount) * sizeof(int));
  (*ys) = (unsigned int *) malloc((*locationCount) * sizeof(int));
  printf("readFile: Reading %u locations...\n", (*locationCount));
  for (int i = 0; i < (*locationCount); ++i) {
    fgets(buffer, BUFFER_LENGTH, fh);
    (*xs)[i] = (unsigned int) strtol(buffer, &end, 10);
    (*ys)[i] = (unsigned int) strtol(end + 1, NULL, 10);
  }
  fclose(fh);
}

float distance(int x1, int y1, int x2, int y2) {
  if ((x1 == x2) && (y1 == y2)) {
    return 0.0;
  }
  int dx = x2 - x1;
  int dy = y2 - y1;
  return sqrtf(dx * dx + dy * dy);
}

void swap(float *distances, unsigned int *distanceIndices, size_t i, size_t j) {
  float tempDistance = distances[i];
  unsigned int tempIndex = distanceIndices[i];
  distances[i] = distances[j];
  distances[j] = tempDistance;
  distanceIndices[i] = distanceIndices[j];
  distanceIndices[j] = tempIndex;
}

size_t partition(float *distances,
                 unsigned int *distanceIndices,
                 size_t lo,
                 size_t hi) {
  float pivot = distances[hi];
  size_t i = lo - 1;
  for (size_t j = lo; j < hi; ++j) {
    if (distances[j] < pivot) {
      i++;
      swap(distances, distanceIndices, i, j);
    }
  }
  if (distances[hi] < distances[i + 1]) {
    swap(distances, distanceIndices, i + 1, hi);
  }
  return i + 1;
}

void quicksort(float *distances,
               unsigned int *distanceIndices,
               size_t lo,
               size_t hi) {
  if (lo < hi) {
    size_t p = partition(distances, distanceIndices, lo, hi);
    if (p > 0) {
      quicksort(distances, distanceIndices, lo, p - 1);
    }
    quicksort(distances, distanceIndices, p + 1, hi);
  }
}

#ifndef TESTING

int main(int argc, char *argv[]) {
  char *inputfile, *outputfile;
  unsigned int *xs, *ys, *distanceIndices, locationCount, offset, index;
  size_t distanceCount, k, coordinateSize, distanceSize, distanceIndexSize;
  struct timeval starttime, endtime;
  float *distances;
  FILE *outputfh;
  int i, j, xi, yi, xj, yj, rowIndex, columnIndex;
  double elapsedMs;

  if (argc != 3) {
    printf("usage: %s inputfile outputfile", argv[0]);
    exit(EXIT_FAILURE);
  }

  inputfile = argv[1];
  outputfile = argv[2];
  gettimeofday(&starttime, NULL);

  readFile(inputfile, &xs, &ys, &locationCount);
  printf("main: Read %u point coordinates\n", locationCount);
  coordinateSize = locationCount * sizeof(int);
  printf("main: Allocated 2 coordinate arrays of %lu bytes\n", coordinateSize);

  distanceCount = locationCount * locationCount;
  distanceSize = distanceCount * sizeof(float);
  distances = (float *) malloc(distanceSize);
  if (distances == NULL) {
    perror("malloc");
    exit(EXIT_FAILURE);
  }
  printf("main: Allocated distance array of %lu bytes\n", distanceSize);

  distanceIndexSize = distanceCount * sizeof(unsigned int);
  distanceIndices = (unsigned int *) malloc(distanceIndexSize);
  if (distanceIndices == NULL) {
    perror("malloc");
    exit(EXIT_FAILURE);
  }
  printf("main: Allocated distance index array of %lu bytes\n",
         distanceIndexSize);

  for (i = 0; i < locationCount; ++i) {
    xi = xs[i];
    yi = ys[i];
    offset = i * locationCount;
    for (j = 0; j < locationCount; ++j) {
      xj = xs[j];
      yj = ys[j];
      index = offset + j;
      distances[index] = distance(xi, yi, xj, yj);
      distanceIndices[index] = index;
    }
  }

  printf("main: Computed %lu distances\n", distanceCount);

  quicksort(distances, distanceIndices, 0, distanceCount - 1);

  printf("main: Sorted distances\n");

  outputfh = fopen(outputfile, "w");
  if (outputfh == NULL) {
    perror("writeSortedFile");
    exit(EXIT_FAILURE);
  }
  printf("main: Opened file %s\n", outputfile);
  fprintf(outputfh, "%lu\n", distanceCount);

  for (k = 0; k < distanceCount; ++k) {
    rowIndex = distanceIndices[k] % locationCount;
    columnIndex = distanceIndices[k] / locationCount;
    fprintf(outputfh, "%u;%u;%u;%u;%.3f\n",
            xs[rowIndex],
            ys[rowIndex],
            xs[columnIndex],
            ys[columnIndex],
            distances[k]);
  }
  printf("main: Wrote %lu distances to file %s\n", distanceCount, outputfile);

  fclose(outputfh);
  printf("main: Closed file %s\n", outputfile);

  gettimeofday(&endtime, NULL);
  free(distances);
  free(distanceIndices);
  elapsedMs = (endtime.tv_sec - starttime.tv_sec) * 1000.0
    + (endtime.tv_usec - starttime.tv_usec) / 1000.0;
  printf("main: Sorted %lu distances in %.3f ms\n", distanceCount, elapsedMs);

  return 0;
}

#endif

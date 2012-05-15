#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

//#define INTEGER
#define PKGDATADIR "D:\\rules\\"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#include <string.h>
#include <limits>

#include "util.h"
#include "energy.h"

#ifndef isinf
# define isinf(x) (!isFinite(x) && x == x)
//#define isinf(x) (x==999999)
#endif

/* functions to load energy rules
 * programs should include energy.h and link with energy.o
 */

//const double R = .0019872;
const char BASES[5] = {'A', 'C', 'G', 'U', 'N'};
const char BASE_PAIRS[6][4] = {"A-U", "C-G", "G-C", "U-A", "G-U", "U-G"};


#ifdef INTEGER
#define scale(d) (isinf(d) ? INFINITY : floor((d) * PRECISION + 0.5))
const ENERGY INFINITY = 999999;
#else
#define scale(d) ((d) * PRECISION)
//const ENERGY INFINITY = 1.0 / 0.0;
const ENERGY INFINITY = std::numeric_limits<double>::infinity();
#endif

double ion(int NA, int polymer, double naConc, double mgConc)
{
  if (NA == 0)
    return 0;
  else
    if (polymer)
      {
	if (mgConc != 0.0)
	  fputs("Warning: [Mg++] correction ignored for polymer mode\n", stderr);
	return -(0.2 + 0.175 * log(naConc));
      }
    else
      return -0.114 * log(naConc + 3.3 * sqrt(mgConc));
}

int min(int a, int b)
{
  return (a < b) ? a : b;
}

FILE* openFile(char* name)
{
  FILE* file=NULL;
  char* buffer;

  file = fopen(name, "rt");

  if (!file && getenv("UNAFOLDDAT"))
    {
      buffer = new char[strlen(getenv("UNAFOLDDAT")) + strlen(name) + 2];
	  if (buffer==NULL)//Sorry you ran out of memory.
		  exit(EXIT_FAILURE);
      strcpy(buffer, getenv("UNAFOLDDAT"));
      strcat(buffer, "/");
      strcat(buffer, name);
      file = fopen(buffer, "rt");
      free(buffer);
    }

  if (!file)
    {
      buffer = (char*)malloc(strlen(PKGDATADIR) + strlen(name) + 1);
      strcpy(buffer, PKGDATADIR);
      strcat(buffer, name);
      if (!(file = fopen(buffer, "rt")))
	{
	  perror(name);
	  exit(EXIT_FAILURE);
	}
      free(buffer);
    }

  return file;
}

void loadStack(double stackEnergies[4][4][4][4], double stackEnthalpies[5][5][5][5], int NA, double saltCorrection)
{
  int i, j, ii, jj;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("stack.DGD");
      hFile = openFile("stack.DHD");
    }
  else
    {
      gFile = openFile("stack.DG");
      hFile = openFile("stack.DH");
    }

  for (i = 0; i < 5; ++i)
    for (ii = 0; ii < 5; ++ii)
      for (j = 0; j < 5; ++j)
	for (jj = 0; jj < 5; ++jj)
	  if (i == 4 || j == 4 || ii == 4 || jj == 4)
	    stackEnthalpies[i][j][ii][jj] = INFINITY;
	  else
	    {
	      fscanf( gFile, "%lf", &stackEnergies[i][j][ii][jj]);
	      stackEnergies[i][j][ii][jj] += saltCorrection;
	      fscanf( hFile, "%lf", &stackEnthalpies[i][j][ii][jj]);
	    }

  fclose(gFile);
  fclose(hFile);
}

void combineStack(double stackEnergies[4][4][4][4], double stackEnthalpies[5][5][5][5], double tRatio, ENERGY stack[5][5][5][5])
{
  int i, j, ii, jj;

  for (i = 0; i < 5; ++i)
    for (ii = 0; ii < 5; ++ii)
      for (j = 0; j < 5; ++j)
	for (jj = 0; jj < 5; ++jj)
	  if (i == 4 || j == 4 || ii == 4 || jj == 4)
	    stack[i][j][ii][jj] = INFINITY;
	  else if (!isFinite(stackEnergies[i][j][ii][jj]) || !isFinite(stackEnthalpies[i][j][ii][jj]))
	    stack[i][j][ii][jj] = INFINITY;
	  else
	    stack[i][j][ii][jj] = scale(tRatio * stackEnergies[i][j][ii][jj] + (1.0 - tRatio) * stackEnthalpies[i][j][ii][jj]);
}

void calculateStack(ENERGY stack[5][5][5][5], double tRatio, double scaleFactor)
{
  int i, j, ii, jj;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 5; ++ii)
	for (jj = 0; jj < 5; ++jj)
	  stack[i][j][ii][jj] = exp(-stack[i][j][ii][jj] / RT) / scaleFactor / scaleFactor;
}

void calculateStack2(ENERGY stack[5][5][6][6], double tRatio, double scaleFactor)
{
  int i, j, ii, jj;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  stack[i][j][ii][jj] = exp(-stack[i][j][ii][jj] / RT) / scaleFactor / scaleFactor;
}

void calculateZipStack2(ENERGY stack[5][5][6][6], double tRatio, ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double scaleFactor)
{
  int i, j, ii, jj;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  if (ii == 5 || jj == 5)
	    stack[i][j][ii][jj] = exp(-stack[i][j][ii][jj] / RT) / scaleFactor / scaleFactor;
	  else if (!isFinite(stack[i][j][ii][jj]))
	    stack[i][j][ii][jj] = (-dangle3[i][j][ii] - dangle5[i][j][jj]) / scaleFactor / scaleFactor;
	  else
	    stack[i][j][ii][jj] = (exp(-stack[i][j][ii][jj] / RT) - dangle3[i][j][ii] - dangle5[i][j][jj] - 1) / scaleFactor / scaleFactor;
}

void calculateZeroStack2(ENERGY stack[5][5][6][6])
{
  int i, j, ii, jj;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  stack[i][j][ii][jj] = 0;
}

void calculateInfStack2(ENERGY stack[5][5][6][6])
{
  int i, j, ii, jj;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (ii = 0; ii < 6; ++ii)
	for (jj = 0; jj < 6; ++jj)
	  stack[i][j][ii][jj] = INFINITY;
}

void loadStackSuffix(ENERGY stack[5][5][5][5], char* suffix)
{
  int i, j, ii, jj;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "stack.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i = 0; i < 5; ++i)
    for (ii = 0; ii < 5; ++ii)
      for (j = 0; j < 5; ++j)
	for (jj = 0; jj < 5; ++jj)
	  if (i == 4 || j == 4 || ii == 4 || jj == 4)
	    stack[i][j][ii][jj] = INFINITY;
	  else
	    {
	      fscanf( file, "%lf", &d);
	      stack[i][j][ii][jj] = scale(d);
	    }

  fclose(file);
}

void symmetryCheckStack(double stack[4][4][4][4], char* which)
{
  int i, j, ii, jj;

  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      for (ii = 0; ii < 4; ++ii)
	for (jj = 0; jj < 4; ++jj)
	  if (stack[i][j][ii][jj] != stack[jj][ii][j][i])
	    fprintf(stderr, "Warning: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n", BASES[i], BASES[j], BASES[ii], BASES[jj], which, stack[i][j][ii][jj], BASES[jj], BASES[ii], BASES[j], BASES[i], which, stack[jj][ii][j][i]);
}

double estimateScale(ENERGY stack[5][5][5][5])
{
  int i1, i2;
  double avg = 0.0;

  for (i1 = 0; i1 < 4; ++i1)
    for (i2 = 0; i2 < 4; ++i2)
      avg += stack[i1][3 - i1][i2][3 - i2];

  return avg / 80;
}

void loadDangle(double dangleEnergies3[4][4][4], double dangleEnthalpies3[5][5][6], double dangleEnergies5[4][4][4], double dangleEnthalpies5[5][5][6], int NA, double saltCorrection)
{
  int i, j, k;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("dangle.DGD");
      hFile = openFile("dangle.DHD");
    }
  else
    {
      gFile = openFile("dangle.DG");
      hFile = openFile("dangle.DH");
    }

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  dangleEnthalpies3[i][j][k] = INFINITY;
	else if (k == 4)
	  dangleEnthalpies3[i][j][k] = INFINITY;
	else if (k == 5)
	  dangleEnthalpies3[i][j][k] = INFINITY;
	else
	  {
		fscanf(gFile,"%lf",&dangleEnergies3[i][j][k]);
	    //fscanf( gFile, "%lf", &dangleEnergies3[i][j][k]);
	    dangleEnergies3[i][j][k] += 0.5 * saltCorrection;
	    //fscanf( hFile, "%lf", &dangleEnthalpies3[i][j][k]);
		fscanf(hFile,"%lf",&dangleEnthalpies3[i][j][k]);
	  }

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  dangleEnthalpies5[i][j][k] = INFINITY;
	else if (k == 4)
	  dangleEnthalpies5[i][j][k] = INFINITY;
	else if (k == 5)
	  dangleEnthalpies5[i][j][k] = INFINITY;
	else
	{
		fscanf(gFile,"%lf",&dangleEnergies5[i][j][k]);
		//fscanf( gFile, "%lf", &dangleEnergies5[i][j][k]);
		dangleEnergies5[i][j][k] += 0.5 * saltCorrection;
		//fscanf( hFile, "%lf", &dangleEnthalpies5[i][j][k]);
		fscanf(hFile,"%lf",&dangleEnthalpies5[i][j][k]);
	}

  fclose(gFile);
  fclose(hFile);
}

void combineDangle(double dangleEnergies3[4][4][4], double dangleEnergies5[4][4][4], double dangleEnthalpies3[5][5][6], double dangleEnthalpies5[5][5][6], double tRatio, ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6])
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  if (i == 4 || j == 4)
	    {
	      dangle3[i][j][k] = INFINITY;
	      dangle5[i][j][k] = INFINITY;
	    }
	  else if (k == 4)
	    {
	      dangle3[i][j][k] = INFINITY;
	      dangle5[i][j][k] = INFINITY;
	    }
	  else if (k == 5)
	    {
	      dangle3[i][j][k] = INFINITY;
	      dangle5[i][j][k] = INFINITY;
	    }
	  else
	    {
	      if (!isFinite(dangleEnergies3[i][j][k]) || !isFinite(dangleEnthalpies3[i][j][k]))
		dangle3[i][j][k] = INFINITY;
	      else
		dangle3[i][j][k] = scale(tRatio * dangleEnergies3[i][j][k] + (1.0 - tRatio) * dangleEnthalpies3[i][j][k]);
	      if (!isFinite(dangleEnergies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k]))
		dangle5[i][j][k] = INFINITY;
	      else
		dangle5[i][j][k] = scale(tRatio * dangleEnergies5[i][j][k] + (1.0 - tRatio) * dangleEnthalpies5[i][j][k]);
	    }
	}
}

void combineDangleNew(double dangleEnergies3[4][4][4], double dangleEnergies5[4][4][4], double dangleEnthalpies3[5][5][6], double dangleEnthalpies5[5][5][6], double tRatio, ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6])
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  if (i == 4 || j == 4)
	    {
	      dangle3[i][j][k] = 0;
	      dangle5[i][j][k] = 0;
	    }
	  else if (k == 4)
	    {
	      dangle3[i][j][k] = 0;
	      dangle5[i][j][k] = 0;
	    }
	  else if (k == 5)
	    {
	      dangle3[i][j][k] = 0;
	      dangle5[i][j][k] = 0;
	    }
	  else
	    {
	      if (!isFinite(dangleEnergies3[i][j][k]) || !isFinite(dangleEnthalpies3[i][j][k]))
		dangle3[i][j][k] = 0;
	      else
		{
		  double ea, b;
		  b = -dangleEnthalpies3[i][j][k] / R * exp(-dangleEnergies3[i][j][k] / R / 310.15) / (exp(-dangleEnergies3[i][j][k] / R / 310.15) - 1);
		  ea = (exp(-dangleEnergies3[i][j][k] / R / 310.15) - 1) / exp(b / 310.15);
		  dangle3[i][j][k] = ea * exp(b / tRatio / 310.15);
		}
	      if (!isFinite(dangleEnergies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k]))
		dangle5[i][j][k] = 0;
	      else
		{
		  double ea, b;
		  b = -dangleEnthalpies5[i][j][k] / R * exp(-dangleEnergies5[i][j][k] / R / 310.15) / (exp(-dangleEnergies5[i][j][k] / R / 310.15) - 1);
		  ea = (exp(-dangleEnergies5[i][j][k] / R / 310.15) - 1) / exp(b / 310.15);
		  dangle5[i][j][k] = ea * exp(b / tRatio / 310.15);
		}
	    }
	}
}

void calculateDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double tRatio, double scaleFactor)
{
  int i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  dangle3[i][j][k] = exp(-dangle3[i][j][k] / RT) / scaleFactor;
	  dangle5[i][j][k] = exp(-dangle5[i][j][k] / RT) / scaleFactor;
	}
}

void calculateZipDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double tRatio, double scaleFactor)
{
  int i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (k == 5)
	  {
	    dangle3[i][j][k] = exp(-dangle3[i][j][k] / RT) / scaleFactor;
	    dangle5[i][j][k] = exp(-dangle5[i][j][k] / RT) / scaleFactor;
	  }
	else
	  {
	    dangle3[i][j][k] = (exp(-dangle3[i][j][k] / RT) - 1.0) / scaleFactor;
	    dangle5[i][j][k] = (exp(-dangle5[i][j][k] / RT) - 1.0) / scaleFactor;
	  }
}

void calculateZeroDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6])
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  dangle3[i][j][k] = 0;
	  dangle5[i][j][k] = 0;
	}
}

void calculateInfDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6])
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  dangle3[i][j][k] = INFINITY;
	  dangle5[i][j][k] = INFINITY;
	}
}

void loadDangleSuffix(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], char* suffix)
{
  int i, j, k;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(7 + strlen(suffix) + 1);
  strcpy(buffer, "dangle.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  dangle3[i][j][k] = INFINITY;
	else if (k == 4)
	  dangle3[i][j][k] = INFINITY;
	else if (k == 5)
	  dangle3[i][j][k] = INFINITY;
	else
	  {
	    fscanf( file, "%lf", &d);
	    dangle3[i][j][k] = scale(d);
	  }
    
  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	if (i == 4 || j == 4)
	  dangle5[i][j][k] = INFINITY;
	else if (k == 4)
	  dangle5[i][j][k] = INFINITY;
	else if (k == 5)
	  dangle5[i][j][k] = INFINITY;
	else
	  {
	    fscanf( file, "%lf", &d);
	    dangle5[i][j][k] = scale(d);
	  }
    
  fclose(file);
}

void zipDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6])
{
  int i, j;
 
  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      dangle5[i][j][5] = dangle3[i][j][5] = 0.0;
}

void addZeroDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double tRatio)
{
  int i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  dangle3[i][j][k] = -RT * log(1 + exp(-dangle3[i][j][k] / RT));
	  dangle5[i][j][k] = -RT * log(1 + exp(-dangle5[i][j][k] / RT));
	}
}

void minZeroDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6])
{
  int i, j, k;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      for (k = 0; k < 6; ++k)
	{
	  if (dangle3[i][j][k] > 0.0)
	    dangle3[i][j][k] = 0.0;
	  if (dangle5[i][j][k] > 0.0)
	    dangle5[i][j][k] = 0.0;
	}
}

void loadLoop(double hairpinLoopEnergies[30], double interiorLoopEnergies[30], double bulgeLoopEnergies[30], double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], int NA, double saltCorrection)
{
  int k;
  /*FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("loop.DGD");
      hFile = openFile("loop.DHD");
    }
  else
    {
      gFile = openFile("loop.DG");
      hFile = openFile("loop.DH");
    }*/

  for (k = 0; k < 30; ++k)
    {
	//fscanf( gFile, "%*f%lf%lf%lf", &interiorLoopEnergies[k], &bulgeLoopEnergies[k], &hairpinLoopEnergies[k]);
      bulgeLoopEnergies[k] += saltCorrection * (1.0 + 0.5 * min(k, 10));
      interiorLoopEnergies[k] += saltCorrection * (1.0 + 0.5 * min(k, 10));
    //fscanf( hFile, "%*f%lf%lf%lf", &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k]);
    }

  /*fclose(gFile);
  fclose(hFile);*/
}

void combineLoop(double hairpinLoopEnergies[30], double interiorLoopEnergies[30], double bulgeLoopEnergies[30], double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], double tRatio, ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30])
{
  int k;

  for (k = 0; k < 30; ++k)
    {
      hairpinLoop[k] = scale(tRatio * hairpinLoopEnergies[k] + (1.0 - tRatio) * hairpinLoopEnthalpies[k]);
      interiorLoop[k] = scale(tRatio * interiorLoopEnergies[k] + (1.0 - tRatio) * interiorLoopEnthalpies[k]);
      bulgeLoop[k] = scale(tRatio * bulgeLoopEnergies[k] + (1.0 - tRatio) * bulgeLoopEnthalpies[k]);
    }
}

void calculateLoop(ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30], double tRatio, double scaleFactor)
{
  int k;
  const double RT = tRatio * 310.15 * R;

  for (k = 0; k < 30; ++k)
    {
      hairpinLoop[k] = exp(-hairpinLoop[k] / RT) / pow(scaleFactor, k + 3);
      interiorLoop[k] = exp(-interiorLoop[k] / RT) / pow(scaleFactor, k + 3);
      bulgeLoop[k] = exp(-bulgeLoop[k] / RT) / pow(scaleFactor, k + 3);
    }
}

void loadLoopSuffix(ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30], char* suffix)
{
  int k;
  double d1, d2, d3;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(5 + strlen(suffix) + 1);
  strcpy(buffer, "loop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (k = 0; k < 30; ++k)
    {
      fscanf( file, "%*f%lf%lf%lf", &d1, &d2, &d3);
      interiorLoop[k] = scale(d1);
      bulgeLoop[k] = scale(d2);
      hairpinLoop[k] = scale(d3);
    }

  fclose(file);
}

void loadSint2(double sint2Energies[6][6][4][4], double sint2Enthalpies[7][7][5][5], int NA, double saltCorrection)
{
  int b, c, i, j;
  FILE* gFile, *hFile;

  if (NA)
    {
      gFile = openFile("sint2.DGD");
      hFile = openFile("sint2.DHD");
    }
  else
    {
      gFile = openFile("sint2.DG");
      hFile = openFile("sint2.DH");
    }

  for (b = 0; b < 7; ++b)
    for (i = 0; i < 5; ++i)
      for (c = 0; c < 7; ++c)
	for (j = 0; j < 5; ++j)
	  if (b == 6 || c == 6)
	    sint2Enthalpies[b][c][i][j] = INFINITY;
	  else if (i == 4 || j == 4)
	    sint2Enthalpies[b][c][i][j] = 0;
	  else
	    {
	      fscanf( gFile, "%lf", &sint2Energies[b][c][i][j]);
	      sint2Energies[b][c][i][j] += 2.0 * saltCorrection;
	      fscanf( hFile, "%lf", &sint2Enthalpies[b][c][i][j]);
	    }

  fclose(gFile);
  fclose(hFile);
}

void combineSint2(double sint2Energies[6][6][4][4], double sint2Enthalpies[7][7][5][5], double tRatio, ENERGY sint2[7][7][5][5])
{
  int b, c, i, j;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  {
	    if (b == 6 || c == 6)
	      sint2[b][c][i][j] = INFINITY;
	    else if (i == 4 || j == 4)
	      sint2[b][c][i][j] = 0;
	    else if (!isFinite(sint2Energies[b][c][i][j]) || !isFinite(sint2Enthalpies[b][c][i][j]))
	      sint2[b][c][i][j] = INFINITY;
	    else
	      sint2[b][c][i][j] = scale(tRatio * sint2Energies[b][c][i][j] + (1.0 - tRatio) * sint2Enthalpies[b][c][i][j]);
	  }
}

void calculateSint2(ENERGY sint2[7][7][5][5], double tRatio, double scaleFactor)
{
  int b, c, i, j;
  const double RT = tRatio * 310.15 * R;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  sint2[b][c][i][j] = exp(-sint2[b][c][i][j] / RT) / pow(scaleFactor, 4);
}

void loadSint2Suffix(ENERGY sint2[7][7][5][5], char* suffix)
{
  int b, c, i, j;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "sint2.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (b = 0; b < 6; ++b)
    for (i = 0; i < 5; ++i)
      for (c = 0; c < 6; ++c)
	for (j = 0; j < 5; ++j)
	  if (b == 6 || c == 6)
	    sint2[b][c][i][j] = INFINITY;
	  else if (i == 4 || j == 4)
	    sint2[b][c][i][j] = 0;
	  else
	    {
	      fscanf( file, "%lf", &d);
	      sint2[b][c][i][j] = scale(d);
	    }

  fclose(file);
}

void symmetryCheckSint2(double sint2[6][6][4][4], char* which)
{
  int b, c, i, j;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i = 0; i < 4; ++i)
	for (j = 0; j < 4; ++j)
	  {
	    int bb, cc;
	    bb = (b > 3) ? 9 - b : 3 - b;
	    cc = (c > 3) ? 9 - c : 3 - c;

	    if (sint2[b][c][i][j] != sint2[cc][bb][j][i])
	      fprintf(stderr, "Warning: %s/%s/%c/%c sint2 %s is %g; %s/%s/%c/%c sint2 %s is %g\n", BASE_PAIRS[b], BASE_PAIRS[c], BASES[i], BASES[j], which, sint2[b][c][i][j], BASE_PAIRS[cc], BASE_PAIRS[bb], BASES[j], BASES[i], which, sint2[cc][bb][j][i]);
	  }
}

void loadAsint1x2(double asint1x2Energies[6][6][4][4][4], double asint1x2Enthalpies[7][7][5][5][5], int NA, double saltCorrection)
{
  int b, c, i, j, k;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("asint1x2.DGD");
      hFile = openFile("asint1x2.DHD");
    }
  else
    {
      gFile = openFile("asint1x2.DG");
      hFile = openFile("asint1x2.DH");
    }

  for (b = 0; b < 7; ++b)
    for (k = 0; k < 5; ++k)
      for (i = 0; i < 5; ++i)
	for (c = 0; c < 7; ++c)
	  for (j = 0; j < 5; ++j)
	    if (b == 6 || c == 6)
	      asint1x2Enthalpies[b][c][i][j][k] = INFINITY;
	    else if (i == 4 || j == 4 || k == 4)
	      asint1x2Enthalpies[b][c][i][j][k] = 0;
	    else
	      {
		fscanf( gFile, "%lf", &asint1x2Energies[b][c][i][j][k]);
		asint1x2Energies[b][c][i][j][k] += 2.5 * saltCorrection;
		fscanf( hFile, "%lf", &asint1x2Enthalpies[b][c][i][j][k]);
	      }

  fclose(gFile);
  fclose(hFile);
}

void combineAsint1x2(double asint1x2Energies[6][6][4][4][4], double asint1x2Enthalpies[7][7][5][5][5], double tRatio, ENERGY asint1x2[7][7][5][5][5])
{
  int b, c, i, j, k;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  for (k = 0; k < 5; ++k)
	    {
	      if (b == 6 || c == 6)
		asint1x2[b][c][i][j][k] = INFINITY;
	      else if (i == 4 || j == 4 || k == 4)
		asint1x2[b][c][i][j][k] = 0;
	      else if (!isFinite(asint1x2Energies[b][c][i][j][k]) || !isFinite(asint1x2Enthalpies[b][c][i][j][k]))
		asint1x2[b][c][i][j][k] = INFINITY;
	      else
		asint1x2[b][c][i][j][k] = scale(tRatio * asint1x2Energies[b][c][i][j][k] + (1.0 - tRatio) * asint1x2Enthalpies[b][c][i][j][k]);
	    }
}

void calculateAsint1x2(ENERGY asint1x2[7][7][5][5][5], double tRatio, double scaleFactor)
{
  int b, c, i, j, k;
  const double RT = tRatio * 310.15 * R;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i = 0; i < 5; ++i)
	for (j = 0; j < 5; ++j)
	  for (k = 0; k < 5; ++k)
	    asint1x2[b][c][i][j][k] = exp(-asint1x2[b][c][i][j][k] / RT) / pow(scaleFactor, 5);
}

void loadAsint1x2Suffix(ENERGY asint1x2[7][7][5][5][5], char* suffix)
{
  int b, c, i, j, k;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "asint1x2.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (b = 0; b < 7; ++b)
    for (k = 0; k < 5; ++k)
      for (i = 0; i < 5; ++i)
	for (c = 0; c < 7; ++c)
	  for (j = 0; j < 5; ++j)
	    if (b == 6 || c == 6)
	      asint1x2[b][c][i][j][k] = INFINITY;
	    else if (i == 4 || j == 4 || k == 4)
	      asint1x2[b][c][i][j][k] = 0;
	    else
	      {
		fscanf( file, "%lf", &d);
		asint1x2[b][c][i][j][k] = scale(d);
	      }

  fclose(file);
}

void loadSint4(double sint4Energies[6][6][4][4][4][4], double sint4Enthalpies[7][7][5][5][5][5], int NA, double saltCorrection)
{
  int b, c, i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("sint4.DGD");
      hFile = openFile("sint4.DHD");
    }
  else
    {
      gFile = openFile("sint4.DG");
      hFile = openFile("sint4.DH");
    }

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      if (b == 6 || c == 6)
		sint4Enthalpies[b][c][i1][j1][i2][j2] = INFINITY;
	      else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
		sint4Enthalpies[b][c][i1][j1][i2][j2] = 0;
	      else
		{
		  fscanf( gFile, "%lf", &sint4Energies[b][c][i1][j1][i2][j2]);
		  sint4Energies[b][c][i1][j1][i2][j2] += 3.0 * saltCorrection;
		  fscanf( hFile, "%lf", &sint4Enthalpies[b][c][i1][j1][i2][j2]);
		}

  fclose(gFile);
  fclose(hFile);
}

void combineSint4(double sint4Energies[6][6][4][4][4][4], double sint4Enthalpies[7][7][5][5][5][5], double tRatio, ENERGY sint4[7][7][5][5][5][5])
{
  int b, c, i1, j1, i2, j2;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      {
		if (b == 6 || c == 6)
		  sint4[b][c][i1][j1][i2][j2] = INFINITY;
		else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
		  sint4[b][c][i1][j1][i2][j2] = 0;
		else if (!isFinite(sint4Energies[b][c][i1][j1][i2][j2]) || !isFinite(sint4Enthalpies[b][c][i1][j1][i2][j2]))
		  sint4[b][c][i1][j1][i2][j2] = INFINITY;
		else
		  sint4[b][c][i1][j1][i2][j2] = scale(tRatio * sint4Energies[b][c][i1][j1][i2][j2] + (1.0 - tRatio) * sint4Enthalpies[b][c][i1][j1][i2][j2]);
	      }
}

void calculateSint4(ENERGY sint4[7][7][5][5][5][5], double tRatio, double scaleFactor)
{
  int b, c, i1, j1, i2, j2;
  const double RT = tRatio * 310.15 * R;

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      sint4[b][c][i1][j1][i2][j2] = exp(-sint4[b][c][i1][j1][i2][j2] / RT) / pow(scaleFactor, 6);
}

void loadSint4Suffix(ENERGY sint4[7][7][5][5][5][5], char* suffix)
{
  int b, c, i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "sint4.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (b = 0; b < 7; ++b)
    for (c = 0; c < 7; ++c)
      for (i1 = 0; i1 < 5; ++i1)
	for (j1 = 0; j1 < 5; ++j1)
	  for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2)
	      if (b == 6 || c == 6)
		sint4[b][c][i1][j1][i2][j2] = INFINITY;
	      else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
		sint4[b][c][i1][j1][i2][j2] = 0;
	      else
		{
		  fscanf( file, "%lf", &d);
		  sint4[b][c][i1][j1][i2][j2] = scale(d);
		}

  fclose(file);
}

void symmetryCheckSint4(double sint4[6][6][4][4][4][4], char* which)
{
  int b, c, i1, j1,i2, j2;

  for (b = 0; b < 6; ++b)
    for (c = 0; c < 6; ++c)
      for (i1 = 0; i1 < 4; ++i1)
	for (j1 = 0; j1 < 4; ++j1)
	  for (i2 = 0; i2 < 4; ++i2)
	    for (j2 = 0; j2 < 4; ++j2)
	      {
		int bb, cc;
		bb = (b > 3) ? 9 - b : 3 - b;
		cc = (c > 3) ? 9 - c : 3 - c;

		if (sint4[b][c][i1][j1][i2][j2] != sint4[cc][bb][j2][i2][j1][i1])
		  fprintf(stderr, "Warning: %s/%s/%c/%c/%c/%c sint4 %s is %g; %s/%s/%c/%c/%c/%c sint4 %s is %g\n", BASE_PAIRS[b], BASE_PAIRS[c], BASES[i1], BASES[j1], BASES[i2], BASES[j2], which, sint4[b][c][i1][j1][i2][j2], BASE_PAIRS[cc], BASE_PAIRS[bb], BASES[j2], BASES[i2], BASES[j1], BASES[i1], which, sint4[cc][bb][j2][i2][j1][i1]);
	      }
}

void loadTstacki(double tstackiEnergies[4][4][4][4], double tstackiEnthalpies[5][5][5][5], int NA)
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("tstacki.DGD");
      hFile = openFile("tstacki.DHD");
    }
  else
    {
      gFile = openFile("tstacki.DG");
      hFile = openFile("tstacki.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackiEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackiEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( gFile, "%lf", &tstackiEnergies[i1][j1][i2][j2]);
	      fscanf( hFile, "%lf", &tstackiEnthalpies[i1][j1][i2][j2]);
	    }

  fclose(gFile);
  fclose(hFile);
}

void loadTstacki23(double tstacki23Energies[4][4][4][4], double tstacki23Enthalpies[5][5][5][5], int NA)
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("tstacki23.DGD");
      hFile = openFile("tstacki23.DHD");
    }
  else
    {
      gFile = openFile("tstacki23.DG");
      hFile = openFile("tstacki23.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstacki23Enthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstacki23Enthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( gFile, "%lf", &tstacki23Energies[i1][j1][i2][j2]);
	      fscanf( hFile, "%lf", &tstacki23Enthalpies[i1][j1][i2][j2]);
	    }

  fclose(gFile);
  fclose(hFile);
}

void loadTstackh(double tstackhEnergies[4][4][4][4], double tstackhEnthalpies[5][5][5][5], int NA)
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("tstackh.DGD");
      hFile = openFile("tstackh.DHD");
    }
  else
    {
      gFile = openFile("tstackh.DG");
      hFile = openFile("tstackh.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackhEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackhEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( gFile, "%lf", &tstackhEnergies[i1][j1][i2][j2]);
	      fscanf( hFile, "%lf", &tstackhEnthalpies[i1][j1][i2][j2]);
	    }

  fclose(gFile);
  fclose(hFile);
}

void loadTstackm(double tstackmEnergies[4][4][4][4], double tstackmEnthalpies[5][5][6][6], int NA, double saltCorrection)
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("tstackm.DGD");
      hFile = openFile("tstackm.DHD");
    }
  else
    {
      gFile = openFile("tstackm.DG");
      hFile = openFile("tstackm.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackmEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    tstackmEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackmEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( gFile, "%lf", &tstackmEnergies[i1][j1][i2][j2]);
	      tstackmEnergies[i1][j1][i2][j2] += saltCorrection;
	      fscanf( hFile, "%lf", &tstackmEnthalpies[i1][j1][i2][j2]);
	    }

  fclose(gFile);
  fclose(hFile);
}

void loadTstacke(double tstackeEnergies[4][4][4][4], double tstackeEnthalpies[5][5][6][6], int NA, double saltCorrection)
{
  int i1, j1, i2, j2;
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("tstacke.DGD");
      hFile = openFile("tstacke.DHD");
    }
  else
    {
      gFile = openFile("tstacke.DG");
      hFile = openFile("tstacke.DH");
    }

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackeEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    tstackeEnthalpies[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackeEnthalpies[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( gFile, "%lf", &tstackeEnergies[i1][j1][i2][j2]);
	      tstackeEnergies[i1][j1][i2][j2] += saltCorrection;
	      fscanf( hFile, "%lf", &tstackeEnthalpies[i1][j1][i2][j2]);
	    }

  fclose(gFile);
  fclose(hFile);
}

void combineTstack(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][5][5], double tRatio, ENERGY tstack[5][5][5][5])
{
  int i1, j1, i2, j2;

  for (i1 = 0; i1 < 5; ++i1)
    for (j1 = 0; j1 < 5; ++j1)
      for (i2 = 0; i2 < 5; ++i2)
	for (j2 = 0; j2 < 5; ++j2)
	  {
	    if (i1 == 4 || j1 == 4)
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else if (i2 == 4 || j2 == 4)
	      tstack[i1][j1][i2][j2] = 0;
	    else if (!isFinite(tstackEnergies[i1][j1][i2][j2]) || !isFinite(tstackEnthalpies[i1][j1][i2][j2]))
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else
	      tstack[i1][j1][i2][j2] = scale(tRatio * tstackEnergies[i1][j1][i2][j2] + (1.0 - tRatio) * tstackEnthalpies[i1][j1][i2][j2]);
	  }
}

void combineTstack2(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][6][6], double tRatio, ENERGY tstack[5][5][6][6])
{
  int i1, j1, i2, j2;

  for (i1 = 0; i1 < 5; ++i1)
    for (j1 = 0; j1 < 5; ++j1)
      for (i2 = 0; i2 < 6; ++i2)
	for (j2 = 0; j2 < 6; ++j2)
	  {
	    if (i1 == 4 || j1 == 4)
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else if (i2 == 5 || j2 == 5)
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else if (i2 == 4 || j2 == 4)
	      tstack[i1][j1][i2][j2] = 0;
	    else if (!isFinite(tstackEnergies[i1][j1][i2][j2]) || !isFinite(tstackEnthalpies[i1][j1][i2][j2]))
	      tstack[i1][j1][i2][j2] = INFINITY;
	    else
	      tstack[i1][j1][i2][j2] = scale(tRatio * tstackEnergies[i1][j1][i2][j2] + (1.0 - tRatio) * tstackEnthalpies[i1][j1][i2][j2]);
	  }
}

void loadTstackiSuffix(ENERGY tstacki[5][5][5][5], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstacki.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstacki[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstacki[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      tstacki[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadTstacki23Suffix(ENERGY tstacki23[5][5][5][5], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(10 + strlen(suffix) + 1);
  strcpy(buffer, "tstacki23.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstacki23[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstacki23[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      tstacki23[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadTstackhSuffix(ENERGY tstackh[5][5][5][5], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstackh.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackh[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackh[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      tstackh[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadTstackmSuffix(ENERGY tstackm[5][5][6][6], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstackm.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackm[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    tstackm[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackm[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      tstackm[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadTstackeSuffix(ENERGY tstacke[5][5][6][6], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "tstacke.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 6; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 6; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstacke[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 5 || j2 == 5)
	    tstacke[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstacke[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      tstacke[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadCoaxialSuffix(ENERGY coaxial[5][5][5][5], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "coaxial.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    coaxial[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    coaxial[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      coaxial[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadTstackcoaxSuffix(ENERGY tstackcoax[5][5][5][5], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(11 + strlen(suffix) + 1);
  strcpy(buffer, "tstackcoax.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    tstackcoax[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    tstackcoax[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      tstackcoax[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadCoaxstackSuffix(ENERGY coaxstack[5][5][5][5], char* suffix)
{
  int i1, j1, i2, j2;
  double d;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(10 + strlen(suffix) + 1);
  strcpy(buffer, "coaxstack.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  for (i1 = 0; i1 < 5; ++i1)
    for (i2 = 0; i2 < 5; ++i2)
      for (j1 = 0; j1 < 5; ++j1)
	for (j2 = 0; j2 < 5; ++j2)
	  if (i1 == 4 || j1 == 4)
	    coaxstack[i1][j1][i2][j2] = INFINITY;
	  else if (i2 == 4 || j2 == 4)
	    coaxstack[i1][j1][i2][j2] = 0;
	  else
	    {
	      fscanf( file, "%lg", &d);
	      coaxstack[i1][j1][i2][j2] = scale(d);
	    }

  fclose(file);
}

void loadMulti(double multiEnergies[3], double multiEnthalpies[3], int NA)
{
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("miscloop.DGD");
      hFile = openFile("miscloop.DHD");
    }
  else
    {
      gFile = openFile("miscloop.DG");
      hFile = openFile("miscloop.DH");
    }

  fscanf( gFile, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnergies[0], &multiEnergies[1], &multiEnergies[2]);
  fscanf( hFile, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnthalpies[0], &multiEnthalpies[1], &multiEnthalpies[2]);

  fclose(gFile);
  fclose(hFile);
}

void loadMulti2(double multiEnergies[3], double multiEnthalpies[3], int NA)
{
  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("miscloop.DGD");
      hFile = openFile("miscloop.DHD");
    }
  else
    {
      gFile = openFile("miscloop.DG");
      hFile = openFile("miscloop.DH");
    }

  fscanf( gFile, "%*g%*g%*g%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnergies[0], &multiEnergies[1], &multiEnergies[2]);
  fscanf( hFile, "%*g%*g%*g%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnthalpies[0], &multiEnthalpies[1], &multiEnthalpies[2]);

  fclose(gFile);
  fclose(hFile);
}

void combineMulti(double multiEnergies[3], double multiEnthalpies[3], double tRatio, ENERGY multi[3])
{
  int i;

  for (i = 0; i < 3; ++i)
    multi[i] = scale(tRatio * multiEnergies[i] + (1.0 - tRatio) * multiEnthalpies[i]);
}

void calculateMulti(ENERGY multi[3], double tRatio, double scaleFactor)
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 3; ++i)
    multi[i] = exp(-multi[i] / RT);
  multi[0] /= scaleFactor * scaleFactor;
}

void loadMultiSuffix(ENERGY multi[3], char* suffix)
{
  double d1, d2, d3;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  fscanf( file, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &d1, &d2, &d3);
  multi[0] = scale(d1);
  multi[1] = scale(d2);
  multi[2] = scale(d3);

  fclose(file);
}

void loadMulti2Suffix(ENERGY multi[3], char* suffix)
{
  double d1, d2, d3;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  fscanf( file, "%*g%*g%*g%*g%*g%*g%*g%*g%*g%lg%lg%lg", &d1, &d2, &d3);
  multi[0] = scale(d1);
  multi[1] = scale(d2);
  multi[2] = scale(d3);

  fclose(file);
}

void loadMisc(double miscEnergies[13], double miscEnthalpies[13], int NA)
{
/*  FILE *gFile, *hFile;

  if (NA)
    {
      gFile = openFile("miscloop.DGD");
      hFile = openFile("miscloop.DHD");
    }
  else
    {
      gFile = openFile("miscloop.DG");
      hFile = openFile("miscloop.DH");
    }

  readOrDie(13, "miscloop", gFile, "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg", &miscEnergies[12], &miscEnergies[4], &miscEnergies[0], &miscEnergies[1], &miscEnergies[2], &miscEnergies[3], &miscEnergies[6], &miscEnergies[8], &miscEnergies[9], &miscEnergies[10], &miscEnergies[11], &miscEnergies[5], &miscEnergies[7]);
  readOrDie(13, "miscloop", hFile, "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg", &miscEnthalpies[12], &miscEnthalpies[4], &miscEnthalpies[0], &miscEnthalpies[1], &miscEnthalpies[2], &miscEnthalpies[3], &miscEnthalpies[6], &miscEnthalpies[8], &miscEnthalpies[9], &miscEnthalpies[10], &miscEnthalpies[11], &miscEnthalpies[5], &miscEnthalpies[7]);

  fclose(gFile);
  fclose(hFile);*/
}

void combineMisc(double miscEnergies[13], double miscEnthalpies[13], double tRatio, ENERGY misc[13])
{
  int i;

  for (i = 0; i < 7; ++i)
    misc[i] = scale(tRatio * miscEnergies[i] + (1.0 - tRatio) * miscEnthalpies[i]);
  misc[7] = miscEnergies[7] == 1 || miscEnthalpies[7] == 1;
  for (i = 8; i < 13; ++i)
    misc[i] = scale(tRatio * miscEnergies[i] + (1.0 - tRatio) * miscEnthalpies[i]);
}

void calculateMisc(ENERGY misc[13], double tRatio)
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < 7; ++i)
    misc[i] = exp(-misc[i] / RT);
  for (i = 8; i < 13; ++i)
    misc[i] = exp(-misc[i] / RT);
}

void loadMiscSuffix(ENERGY misc[13], char* suffix)
{
  int i;
  double d[13];
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  fscanf(file, "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg", &d[12], &d[4], &d[0], &d[1], &d[2], &d[3], &d[6], &d[8], &d[9], &d[10], &d[11], &d[5], &d[7]);
  for (i = 0; i < 13; ++i)
    misc[i] = scale(d[i]);

  fclose(file);
}

void makeAUPenalty(ENERGY misc[13], ENERGY aup[5][5], int isPF)
{
  int i, j;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      aup[i][j] = isPF ? 1.0 : 0.0;

  aup[0][3] = aup[3][0] = aup[2][3] = aup[3][2] = misc[6];
}

void makeAUPenaltyH(double misc[13], double aup[5][5], int isPF)
{
  int i, j;

  for (i = 0; i < 5; ++i)
    for (j = 0; j < 5; ++j)
      aup[i][j] = isPF ? 1.0 : 0.0;

  aup[0][3] = aup[3][0] = aup[2][3] = aup[3][2] = misc[6];
}

void loadTriloop(struct triloopE** triloopEnergies, struct triloopE** triloopEnthalpies, int* num, int NA)
{
  FILE *gFile, *hFile;
  int i, size;
  double energy;

  if (NA)
    gFile = openFile("triloop.DGD");
  else
    gFile = openFile("triloop.DG");

  *num = 0;
  size = 16;
  *triloopEnergies = (triloopE*)calloc(16, sizeof(struct triloopE));

  while (fscanf(gFile, "%5s %lg", (*triloopEnergies)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 5; ++i)
	(*triloopEnergies)[*num].loop[i] = util::toNum((*triloopEnergies)[*num].loop[i]);
      (*triloopEnergies)[*num].energy = energy;
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *triloopEnergies = (triloopE*)realloc(*triloopEnergies, size * sizeof(struct triloopE));
	}
    }

  *triloopEnergies = (triloopE*)realloc(*triloopEnergies, *num * sizeof(struct triloopE));

  fclose(gFile);

  if (NA)
    hFile = openFile("triloop.DHD");
  else
    hFile = openFile("triloop.DH");

  *num = 0;
  size = 16;
  *triloopEnthalpies = (triloopE*)calloc(16, sizeof(struct triloopE));

  while (fscanf(hFile, "%5s %lg", (*triloopEnthalpies)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 5; ++i)
	(*triloopEnthalpies)[*num].loop[i] = util::toNum((*triloopEnthalpies)[*num].loop[i]);
      (*triloopEnthalpies)[*num].energy = energy;
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *triloopEnthalpies = (triloopE*)realloc(*triloopEnthalpies, size * sizeof(struct triloopE));
	}
    }

  *triloopEnthalpies = (triloopE*)realloc(*triloopEnthalpies, *num * sizeof(struct triloopE));

  fclose(hFile);
}

void combineTriloop(const struct triloopE* triloopEnergies, const struct triloopE* triloopEnthalpies, double tRatio, struct triloop* triloop, int num)
{
  int i;

  for (i = 0; i < num; ++i)
    {
      memcpy(triloop[i].loop, triloopEnergies[i].loop, 5);
      triloop[i].energy = scale(tRatio * triloopEnergies[i].energy + (1.0 - tRatio) * triloopEnthalpies[i].energy);
    }
}

void calculateTriloop(struct triloop* triloop, int num, double tRatio)
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < num; ++i)
    triloop[i].energy = exp(-triloop[i].energy / RT);
}

void loadTriloopSuffix(struct triloop** triloop_v, int* num, char* suffix)
{
  int i, size;
  double energy;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(8 + strlen(suffix) + 1);
  strcpy(buffer, "triloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  *num = 0;
  size = 16;
  *triloop_v = (triloop*)calloc(16, sizeof(struct triloop));

  while (fscanf(file, "%5s %lg", (*triloop_v)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 5; ++i)
	(*triloop_v)[*num].loop[i] = util::toNum((*triloop_v)[*num].loop[i]);
      (*triloop_v)[*num].energy = scale(energy);
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *triloop_v = (triloop*)realloc(*triloop_v, size * sizeof(struct triloop));
	}
    }

  *triloop_v = (triloop*)realloc(*triloop_v, *num * sizeof(struct triloop));

  fclose(file);
}

int triloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = (const unsigned char*)loop1;
  const struct triloop *h2 = (const triloop*)loop2;

  for (i = 0; i < 5; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}

void loadTloop(struct tloopE** tloopEnergies, struct tloopE** tloopEnthalpies, int* num, int NA)
{
  FILE *gFile, *hFile;
  int i, size;
  double energy;

  if (NA)
    gFile = openFile("tloop.DGD");
  else
    gFile = openFile("tloop.DG");

  *num = 0;
  size = 16;
  *tloopEnergies = (tloopE*)calloc(16, sizeof(struct tloopE));

  while (fscanf(gFile, "%6s %lg", (*tloopEnergies)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 6; ++i)
	(*tloopEnergies)[*num].loop[i] = util::toNum((*tloopEnergies)[*num].loop[i]);
      (*tloopEnergies)[*num].energy = energy;
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *tloopEnergies = (tloopE*)realloc(*tloopEnergies, size * sizeof(struct tloopE));
	}
    }

  *tloopEnergies = (tloopE*)realloc(*tloopEnergies, *num * sizeof(struct tloopE));

  fclose(gFile);

  if (NA)
    hFile = openFile("tloop.DHD");
  else
    hFile = openFile("tloop.DH");

  *num = 0;
  size = 16;
  *tloopEnthalpies = (tloopE*)calloc(16, sizeof(struct tloopE));

  while (fscanf(hFile, "%6s %lg", (*tloopEnthalpies)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 6; ++i)
	(*tloopEnthalpies)[*num].loop[i] = util::toNum((*tloopEnthalpies)[*num].loop[i]);
      (*tloopEnthalpies)[*num].energy = energy;
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *tloopEnthalpies = (tloopE*)realloc(*tloopEnthalpies, size * sizeof(struct tloopE));
	}
    }

  *tloopEnthalpies = (tloopE*)realloc(*tloopEnthalpies, *num * sizeof(struct tloopE));

  fclose(hFile);
}

void combineTloop(const struct tloopE* tloopEnergies, const struct tloopE* tloopEnthalpies, double tRatio, struct tloop* tloop, int num)
{
  int i;

  for (i = 0; i < num; ++i)
    {
      memcpy(tloop[i].loop, tloopEnergies[i].loop, 6);
      tloop[i].energy = scale(tRatio * tloopEnergies[i].energy + (1.0 - tRatio) * tloopEnthalpies[i].energy);
    }
}

void calculateTloop(struct tloop* tloop, int num, double tRatio)
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < num; ++i)
    tloop[i].energy = exp(-tloop[i].energy / RT);
}

void loadTloopSuffix(struct tloop** tloop_v, int* num, char* suffix)
{
  int i, size;
  double energy;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(6 + strlen(suffix) + 1);
  strcpy(buffer, "tloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  *num = 0;
  size = 16;
  *tloop_v = (tloop*)calloc(16, sizeof(struct tloop));

  while (fscanf(file, "%6s %lg", (*tloop_v)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 6; ++i)
	(*tloop_v)[*num].loop[i] = util::toNum((*tloop_v)[*num].loop[i]);
      (*tloop_v)[*num].energy = scale(energy);
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *tloop_v = (tloop*)realloc(*tloop_v, size * sizeof(struct tloop));
	}
    }

  *tloop_v = (tloop*)realloc(*tloop_v, *num * sizeof(struct tloop));

  fclose(file);
}

int tloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = (const unsigned char*)loop1;
  const struct tloop *h2 = (const tloop*)loop2;

  for (i = 0; i < 6; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}

void loadHexaloop(struct hexaloopE** hexaloopEnergies, struct hexaloopE** hexaloopEnthalpies, int* num, int NA)
{
  FILE *gFile, *hFile;
  int i, size;
  double energy;

  if (NA)
    gFile = openFile("hexaloop.DGD");
  else
    gFile = openFile("hexaloop.DG");

  *num = 0;
  size = 16;
  *hexaloopEnergies = (hexaloopE*)calloc(16, sizeof(struct hexaloopE));

  while (fscanf(gFile, "%8s %lg", (*hexaloopEnergies)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 8; ++i)
	(*hexaloopEnergies)[*num].loop[i] = util::toNum((*hexaloopEnergies)[*num].loop[i]);
      (*hexaloopEnergies)[*num].energy = energy;
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *hexaloopEnergies = (hexaloopE*)realloc(*hexaloopEnergies, size * sizeof(struct hexaloopE));
	}
    }

  *hexaloopEnergies = (hexaloopE*)realloc(*hexaloopEnergies, *num * sizeof(struct hexaloopE));

  fclose(gFile);

  if (NA)
    hFile = openFile("hexaloop.DHD");
  else
    hFile = openFile("hexaloop.DH");

  *num = 0;
  size = 16;
  *hexaloopEnthalpies = (hexaloopE*)calloc(16, sizeof(struct hexaloopE));

  while (fscanf(hFile, "%8s %lg", (*hexaloopEnthalpies)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 8; ++i)
	(*hexaloopEnthalpies)[*num].loop[i] = util::toNum((*hexaloopEnthalpies)[*num].loop[i]);
      (*hexaloopEnthalpies)[*num].energy = energy;
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *hexaloopEnthalpies = (hexaloopE*)realloc(*hexaloopEnthalpies, size * sizeof(struct hexaloopE));
	}
    }

  *hexaloopEnthalpies = (hexaloopE*)realloc(*hexaloopEnthalpies, *num * sizeof(struct hexaloopE));

  fclose(hFile);
}

void combineHexaloop(const struct hexaloopE* hexaloopEnergies, const struct hexaloopE* hexaloopEnthalpies, double tRatio, struct hexaloop* hexaloop, int num)
{
  int i;

  for (i = 0; i < num; ++i)
    {
      memcpy(hexaloop[i].loop, hexaloopEnergies[i].loop, 8);
      hexaloop[i].energy = scale(tRatio * hexaloopEnergies[i].energy + (1.0 - tRatio) * hexaloopEnthalpies[i].energy);
    }
}

void calculateHexaloop(struct hexaloop* hexaloop, int num, double tRatio)
{
  int i;
  const double RT = tRatio * 310.15 * R;

  for (i = 0; i < num; ++i)
    hexaloop[i].energy = exp(-hexaloop[i].energy / RT);
}

void loadHexaloopSuffix(struct hexaloop** hexaloop_v, int* num, char* suffix)
{
  int i, size;
  double energy;
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "hexaloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  *num = 0;
  size = 16;
  *hexaloop_v = (hexaloop*)calloc(16, sizeof(struct hexaloop));

  while (fscanf(file, "%8s %lg", (*hexaloop_v)[*num].loop, &energy) == 2)
    {
      for (i = 0; i < 8; ++i)
	(*hexaloop_v)[*num].loop[i] = util::toNum((*hexaloop_v)[*num].loop[i]);
      (*hexaloop_v)[*num].energy = scale(energy);
      ++*num;
      if (*num == size)
	{
	  size *= 2;
	  *hexaloop_v = (hexaloop*)realloc(*hexaloop_v, size * sizeof(struct hexaloop));
	}
    }

  *hexaloop_v = (hexaloop*)realloc(*hexaloop_v, *num * sizeof(struct hexaloop));

  fclose(file);
}

int hexaloopcmp(const void* loop1, const void* loop2)
{
  int i;
  const unsigned char* h1 = (const unsigned char*)loop1;
  const struct hexaloop *h2 = (const hexaloop*)loop2;

  for (i = 0; i < 8; ++i)
    if (h1[i] < h2->loop[i])
      return -1;
    else if (h1[i] > h2->loop[i])
      return 1;

  return 0;
}

void loadRTSuffix(double* RT, char* suffix)
{
  char* buffer;
  FILE* file;

  buffer = (char*)malloc(9 + strlen(suffix) + 1);
  strcpy(buffer, "miscloop.");
  strcat(buffer, suffix);
  file = openFile(buffer);
  free(buffer);

  if (fscanf(file, "%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%*f%lf", RT) < 1)
    *RT = 37;

  *RT = R * (*RT + 273.15);

  fclose(file);
}

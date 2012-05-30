#pragma once
#include "energy.h"

class HybridMin
{
public:
	HybridMin(void);
	~HybridMin(void);

	double compute(double& dG, double& dH, const char* sequence1, const char* sequence2, double temperature);
private:
#define Lprime(i, j) lprime[g_len2 * (i) - (j)]
#define Rprime(i, j) rprime[g_len1 * (j) - (i)]

#define auPenaltyH(a, b) g_aupH[a][b]
#define hybAuPenalty(a, b) g_aup[a][b]

int min3(int a, int b, int c)
{
  if (a <= b && a <= c)
    return a;
  if (b <= c)
    return b;
  return c;
};

struct stackNode
{
  int i;
  int j;
  struct stackNode* next;
};

struct constraintListNode
{
  int i, j, k, l;
  struct constraintListNode* next;
} *forceList;
#if ENABLE_FORCE
char *g_ssok1, *g_ssok2;
#define ssOK1(i, j) g_ssok1[(i) * (g_len1 + 2) + j]
#define ssOK2(i, j) g_ssok2[(i) * (g_len2 + 2) + j]
#endif

struct pairListNode
{
  int i, j, length;
  ENERGY E;
  struct pairListNode* next;
} *pairList;

void initializeMatrices();
void force();
void prefilter();
void fillMatrixL(double);
void fillMatrixR(double);
void fillMatrixL_noI(double);
void fillMatrixR_noI(double);
void traceback(int, int, double, int*, int*, int*, int*, int*, int*, double*);
void traceback_noI(int, int, double, int*, int*, int*, int*, int*, int*, double*);
void traceback_rev(int, int, double, int*, int*, int*, int*, int*, int*);
void traceback_rev_noI(int, int, double, int*, int*, int*, int*, int*, int*);
void setStackBI(int, int, int, int, int*, int*, int*, int*);
int unique(int*, char**);
void makePairList(ENERGY, int*, int*);
void makePairList_noI(ENERGY, int*, int*);
ENERGY Es(int, int);
double Hs(int, int);
ENERGY Ebi(int, int, int, int, double);
double Hbi(int, int, int, int);
ENERGY R0(int, int);
ENERGY L0(int, int);

void push(struct stackNode**, int, int);
int equal(ENERGY, ENERGY);
ENERGY* recalloc2(ENERGY*, int, int);
#define min2(a, b) ((a) < (b) ? (a) : (b))
ENERGY min4(ENERGY, ENERGY, ENERGY, ENERGY);
void pushPairList(int, int, int, ENERGY);
void sortPairList();

ENERGY *lprime, *rprime;

int g_allPairs, g_maxLoop, g_nodangle, g_prefilter1, g_prefilter2, g_mfoldMax, g_mfoldP, g_mfoldW, g_quiet, g_zip, g_stream;
int g_numSeqs;
char *g_name1, *g_name2, *g_string1, *g_string2;
unsigned char *g_seq1, *g_seq2; /* [0-4] for [A,C,G,TU,N] */
char *g_file1, *g_file2, *g_prefix, *g_bpFile;
int g_len1, g_len2;
int g_oneTemp, g_firstSeq;
ENERGY g_homodimer;

double dangleEnergies3[4][4][4];
double dangleEnthalpies3[5][5][6];
double dangleEnergies5[4][4][4];
double dangleEnthalpies5[5][5][6];
double stackEnergies[4][4][4][4];
double stackEnthalpies[5][5][5][5];
double interiorLoopEnergies[30];
double bulgeLoopEnergies[30];
double hairpinLoopEnergies[30];
double interiorLoopEnthalpies[30];
double bulgeLoopEnthalpies[30];
double hairpinLoopEnthalpies[30];
double sint2Energies[6][6][4][4];
double sint2Enthalpies[7][7][5][5];
double asint1x2Energies[6][6][4][4][4];
double asint1x2Enthalpies[7][7][5][5][5];
double sint4Energies[6][6][4][4][4][4];
double sint4Enthalpies[7][7][5][5][5][5];
double tstackiEnergies[4][4][4][4];
double tstackiEnthalpies[5][5][5][5];
double tstacki23Energies[4][4][4][4];
double tstacki23Enthalpies[5][5][5][5];
double tstackeEnergies[4][4][4][4];
double tstackeEnthalpies[5][5][6][6];
double miscEnergies[13];
double miscEnthalpies[13];

ENERGY g_dangle3[5][5][6];
ENERGY g_dangle5[5][5][6];
ENERGY g_stack[5][5][5][5];
ENERGY g_interiorLoop[30];
ENERGY g_bulgeLoop[30];
ENERGY g_hairpinLoop[30];
ENERGY g_sint2[7][7][5][5];
ENERGY g_asint1x2[7][7][5][5][5];
ENERGY g_sint4[7][7][5][5][5][5];
ENERGY g_tstacki[5][5][5][5];
ENERGY g_tstacki23[5][5][5][5];
ENERGY g_tstacke[5][5][6][6];
ENERGY g_misc[13];
ENERGY g_aup[5][5];
double g_aupH[5][5];

int NA, polymer, skipTraceback, noIsolate, constraints;
  char* constraintsFile;
  double tMin, tInc, tMax;
  double naConc, mgConc;
  double saltCorrection;


  ENERGY Eleft;/*, Eright; */

  char gotSeqs;
  int count, i, j;
  int bestI, bestJ;
  double t, tRatio, RT;
  char *buffer, *suffix;
  struct constraintListNode* newTop;

};


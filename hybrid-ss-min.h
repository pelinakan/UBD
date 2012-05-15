#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <io.h>
#include <cmath>
#include <limits>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#include "energy.h"

class CtEnergy;

class HybridSSMin {

/* hybrid-ss-min
 * compute minimum energy folding of NA sequence and output as .ct files
 */
public:

#define Q(i, j) q[(g_len - 1) * (i - 1) + j - 1]
#define Qprime(i, j) qprime[(g_len - 1) * (i - 1) + j - 1]
#define QM(i, j) qm[(g_len - 1) * (i - 1) + j - 1]
#define Q5(i) q5[i]
#define Q3(j) q3[j - 1]

struct stackNode
{
  int i;
  int j;
  int matrix; /* [0, 1, 2, 3, 4] ~ [Q', QM, Q, Q5, Q3] */
  struct stackNode* next;
};

struct constraintListNode
{
  int i, j, k, l;
  struct constraintListNode* next;
} *prohibitList, *forceList;
#if ENABLE_FORCE
char* g_ssok;
#define ssOK(i, j) g_ssok[(i) * (g_len + 2) + j]
#else
#define ssOK(i, j) 1
#endif

struct pairListNode
{
  int i, j, length;
  ENERGY E;
  struct pairListNode* next;
} *pairList;

struct stack_node
{
	int i;
	int j;
	int open;
	struct stack_node* next;
};

HybridSSMin();

double computeGibsonFreeEnergy(double&, double&, const char*, double, double);
void initializeMatrices();
void fillMatrices1();
void fillMatrices2();
void computeQ53();
void traceback(int, int, int, int*, int*, int*);
void traceback_noI(int, int, int, int*, int*, int*);
void setStack(int, int, int*, int*);
void setDangle5(int, int*, int*);
void setDangle3(int, int*, int*);
void setBI(int, int, int, int, int*, int*);
void circularize(int*, int*, int*, int);
int unique(int*, char**);
void writeStructure(int*, int*, int*, double, ENERGY);
void writeStructureC(int*, int*, int*, double, ENERGY);
void writePlotExt(int*, double, ENERGY);
void makePairList(ENERGY, int*);
void makePairList_noI(ENERGY, int*);
void writePnum(int*, double);
ENERGY Q5_1(int);
ENERGY Q5_2(int);
ENERGY Q5_3(int);
ENERGY Q5_4(int);
ENERGY Q3_1(int);
ENERGY Q3_2(int);
ENERGY Q3_3(int);
ENERGY Q3_4(int);
ENERGY Ed5(int, int);
ENERGY Ed3(int, int);
ENERGY Etstackm(int, int);
ENERGY Etstacke(int, int);
ENERGY Eh(int, int);
ENERGY Es(int, int);
ENERGY Ebi(int, int, int, int);
#define auPenalty(i, j) g_aup[g_seq[i]][g_seq[j]]
ENERGY QBI(int, int);
ENERGY QBI2(int, int);

ENERGY* recalloc2(ENERGY*, int);
ENERGY* recalloc2_double(ENERGY*, int);
#define min2(a, b) ((a) < (b) ? (a) : (b))
ENERGY min4(ENERGY, ENERGY, ENERGY, ENERGY);
ENERGY min5(ENERGY, ENERGY, ENERGY, ENERGY, ENERGY);

int equal(ENERGY, ENERGY);
int helixLength(int i, int j, ENERGY* qprime);
void push(struct stackNode**, int, int, int);
void pushPairList(int, int, int, ENERGY);
void sortPairList();
int isCircular();
int isHomodimer();
double Ed3(int i, int j, int k);
double Ed5(int i, int j, int k);
double Ee(int i, int j, int ds);
double tstackOrDangle(int i, int j, int external);
double chooseDangle(int a, int b);

void prefilter();

int g_len;
ENERGY *q, *qprime, *qm, *q5, *q3;
double RT;

int g_debug, g_nodangle, g_allPairs, g_maxLoop, g_simple, g_noisolate, g_prefilter1, g_prefilter2, /*g_mfoldMax,*/ g_mfoldP, g_mfoldW, g_maxBP, g_circular, g_quiet;
int g_numSeqs;
const char *g_name, *g_string;
unsigned char* g_seq; /* [0-4] for [A,C,G,TU,N] */
char *g_file, *g_prefix, *g_bpFile;
int g_oneTemp, g_firstSeq;
int g_verbose;
int *g_bp, *g_numbers;
int *g_next, *g_prev;
int *g_upst, *g_dnst;
int g_hasStackingInfo;
stack_node *stack;
/*//For ct-energy...

int g_hasStackingInfo;

stack_node* new_top;
int ds, ss1, ss2, open, is_exterior;
double eloop, etotal;

//All the CT function declarations
int ct_isCircular();
int ct_isHomodimer();
double ct_tstackOrDangle(int i, int j, int external);
double ct_chooseDangle(int a, int b);
double ct_auPenalty(int i, int j);
double ct_Ee(int i, int j, int ds);
double ct_Ebi(int i, int j, int ii, int jj);
double ct_Es(int i, int j);
double ct_Eh(int i, int j);
double ct_Etstacke(int i, int j);
double ct_Etstackm(int i, int j);
double ct_Ed3(int i, int j, int k);
double ct_Ed5(int i, int j, int k);*/

ENERGY g_dangle3[5][5][6];
ENERGY g_dangle5[5][5][6];
ENERGY g_stack[5][5][5][5];
ENERGY g_hairpinLoop[30];
ENERGY g_interiorLoop[30];
ENERGY g_bulgeLoop[30];
ENERGY g_sint2[7][7][5][5];
ENERGY g_asint1x2[7][7][5][5][5];
ENERGY g_sint4[7][7][5][5][5][5];
ENERGY g_tstackh[5][5][5][5];
ENERGY g_tstacki[5][5][5][5];
ENERGY g_tstacki23[5][5][5][5];
ENERGY g_tstackm[5][5][6][6];
ENERGY g_tstacke[5][5][6][6];
struct triloop* g_triloop; 
int numTriloops;
struct tloop* g_tloop; 
int numTloops;
struct hexaloop* g_hexaloop; 
int numHexaloops;
ENERGY g_multi[3];
ENERGY g_misc[13];
ENERGY g_aup[5][5];

private:
	int min3(int a, int b, int c)
{
  if (a <= b && a <= c)
    return a;
  if (b <= c)
    return b;
  return c;
};

	//All the declarations we need!
	double dangleEnergies3[4][4][4];
	double dangleEnthalpies3[5][5][6];
	double dangleEnergies5[4][4][4];
	double dangleEnthalpies5[5][5][6];
	double stackEnergies[4][4][4][4];
	double stackEnthalpies[5][5][5][5];
	//Just put these here. There's no point in reading them from a file
	double sint2Energies[6][6][4][4];
	double sint2Enthalpies[7][7][5][5];
	double asint1x2Energies[6][6][4][4][4];
	double asint1x2Enthalpies[7][7][5][5][5];
	double sint4Energies[6][6][4][4][4][4];
	double sint4Enthalpies[7][7][5][5][5][5];
	double tstackhEnergies[4][4][4][4];
	double tstackhEnthalpies[5][5][5][5];
	double tstackiEnergies[4][4][4][4];
	double tstackiEnthalpies[5][5][5][5];
	double tstacki23Energies[4][4][4][4];
	double tstacki23Enthalpies[5][5][5][5];
	double tstackmEnergies[4][4][4][4];
	double tstackmEnthalpies[5][5][6][6];
	double tstackeEnergies[4][4][4][4];
	double tstackeEnthalpies[5][5][6][6];
	struct triloopE* triloopEnergies;
	struct triloopE* triloopEnthalpies;
	struct tloopE* tloopEnergies;
	struct tloopE* tloopEnthalpies;
	struct hexaloopE* hexaloopEnergies;
	struct hexaloopE* hexaloopEnthalpies;

	CtEnergy* enComputer;
};
//#include "options.h"

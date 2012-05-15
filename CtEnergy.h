#pragma once

class HybridSSMin;
#define ENERGY double

class CtEnergy
{

struct stack_node
{
	int i;
	int j;
	int open;
	struct stack_node* next;
};

public:
	CtEnergy(void);
	~CtEnergy(void);
	double compute(HybridSSMin* hComp,int* bp, int* upst, int* dnst);

int g_nodangle,g_len;
int *g_bp, *g_numbers;
unsigned char* g_seq;
int *g_next, *g_prev;
int *g_upst, *g_dnst;
int g_hasStackingInfo;
stack_node *stack;
stack_node* new_top;
int ds, ss1, ss2, open, is_exterior;
double eloop, etotal;
struct triloop* g_triloop; 
int numTriloops;
struct tloop* g_tloop; 
int numTloops;
struct hexaloop* g_hexaloop; 
int numHexaloops;

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
ENERGY g_multi[3];
ENERGY g_misc[13];

private:
	int min3(int a, int b, int c)
{
  if (a <= b && a <= c)
    return a;
  if (b <= c)
    return b;
  return c;
};


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
double ct_Ed5(int i, int j, int k);
};


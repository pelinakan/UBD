#ifndef UTIL_H
#define UTIL_H

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>

#if HAVE_FCNTL_H
# include <fcntl.h>
#endif

#include <stdlib.h>
#include <string.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define TURN 3 /* minimum size of hairpin loop */

//Fake the stupid xcalloc stuff
#define xrealloc realloc
#define xmalloc malloc
#define xcalloc	calloc

const double R = 0.0019872;

#ifdef NO_GU_BASEPAIRS
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 6, 6, 6},
		       {3, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#else
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 4, 6, 6},
		       {3, 6, 5, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#endif
#define basePairIndex(a, b) BPI[a][b]

class util
{

/* #define NO_GU_BASEPAIRS */
public:
static int roundInt(double d);
static void strcatc(char* str, char c);
static void checkArray(char** array, unsigned int* available, unsigned int used, unsigned int increment);
static unsigned char toNum(char c);
static int seqcmp(unsigned char* seq1, unsigned char* seq2, int length);
static int same(unsigned char* a, unsigned char* b, int len);

};
#endif

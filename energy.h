#ifndef ENERGY_H
#define ENERGY_H

/* include this file to use functions in energy.c
 * then link with energy.o
 */

#include <float.h>
#include <math.h>

#ifndef ENERGY
# define ENERGY double
#endif

#ifndef PRECISION
# define PRECISION 1
#endif

#ifdef INTEGER
# define isFinite(x) (x < INFINITY / 2)
#else
#define finite(x) (_finite(x) && !_isnan(x))
//#define finite(x) (finite(x) && !isnan(x))
#define isFinite(x) (finite(x))
#endif

#ifdef INFINITY
# undef INFINITY
#endif

extern const ENERGY INFINITY;
extern const double R;
extern const char BASES[5];
extern const char BASE_PAIRS[6][4];

double ion(int NA, int polymer, double naConc, double mgConc);

void loadStack(double stackEnergies[4][4][4][4], double stackEnthalpies[5][5][5][5], int NA, double saltCorrection);

void combineStack(double stackEnergies[4][4][4][4], double stackEnthalpies[5][5][5][5], double tRatio, ENERGY stack[5][5][5][5]);

void calculateStack(ENERGY stack[5][5][5][5], double tRatio, double scaleFactor);

void calculateStack2(ENERGY stack[5][5][6][6], double tRatio, double scaleFactor);

void calculateZipStack2(ENERGY stack[5][5][6][6], double tRatio, ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double scaleFactor);

void calculateZeroStack2(ENERGY stack[5][5][6][6]);

void calculateInfStack2(ENERGY stack[5][5][6][6]);

void loadStackSuffix(ENERGY stack[5][5][5][5], char* suffix);

void symmetryCheckStack(double stack[4][4][4][4], char* which);

double estimateScale(ENERGY stack[5][5][5][5]);

void loadDangle(double dangleEnergies3[4][4][4], double dangleEnthalpies3[5][5][6], double dangleEnergies5[4][4][4], double dangleEnthalpies5[5][5][6], int NA, double saltCorrection);

void combineDangle(double dangleEnergies3[4][4][4], double dangleEnergies5[4][4][4], double dangleEnthalpies3[5][5][6], double dangleEnthalpies5[5][5][6], double tRatio, ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6]);

void calculateDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double tRatio, double scaleFactor);

void calculateZipDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double tRatio, double scaleFactor);

void calculateZeroDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6]);

void calculateInfDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6]);

void loadDangleSuffix(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], char* suffix);

void zipDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6]);

void addZeroDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6], double tRatio);

void minZeroDangle(ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6]);

void loadLoop(double hairpinLoopEnergies[30], double interiorLoopEnergies[30], double bulgeLoopEnergies[30], double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], int NA, double saltCorrection);

void combineLoop(double hairpinLoopEnergies[30], double interiorLoopEnergies[30], double bulgeLoopEnergies[30], double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], double tRatio, ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30]);

void calculateLoop(ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30], double tRatio, double scaleFactor);

void loadLoopSuffix(ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30], char* suffix);

void loadSint2(double sint2Energies[6][6][4][4], double sint2Enthalpies[7][7][5][5], int NA, double saltCorrection);

void combineSint2(double sint2Energies[6][6][4][4], double sint2Enthalpies[7][7][5][5], double tRatio, ENERGY sint2[7][7][5][5]);

void calculateSint2(ENERGY sint2[7][7][5][5], double tRatio, double scaleFactor);

void loadSint2Suffix(ENERGY sint2[7][7][5][5], char* suffix);

void symmetryCheckSint2(double sint2[6][6][4][4], char* which);

void loadAsint1x2(double asint1x2Energies[6][6][4][4][4], double asint1x2Enthalpies[7][7][5][5][5], int NA, double saltCorrection);

void combineAsint1x2(double asint1x2Energies[6][6][4][4][4], double asint1x2Enthalpies[7][7][5][5][5], double tRatio, ENERGY asint1x2[7][7][5][5][5]);

void calculateAsint1x2(ENERGY asint1x2[7][7][5][5][5], double tRatio, double scaleFactor);

void loadAsint1x2Suffix(ENERGY asint1x2[7][7][5][5][5], char* suffix);

void loadSint4(double sint4Energies[6][6][4][4][4][4], double sint4Enthalpies[7][7][5][5][5][5], int NA, double saltCorrection);

void combineSint4(double sint4Energies[6][6][4][4][4][4], double sint4Enthalpies[7][7][5][5][5][5], double tRatio, ENERGY sint4[7][7][5][5][5][5]);

void calculateSint4(ENERGY sint4[7][7][5][5][5][5], double tRatio, double scaleFactor);

void loadSint4Suffix(ENERGY sint4[7][7][5][5][5][5], char* suffix);

void symmetryCheckSint4(double sint4[6][6][4][4][4][4], char* which);

void loadTstacki(double tstackiEnergies[4][4][4][4], double tstackiEnthalpies[5][5][5][5], int NA);

void loadTstacki23(double tstacki23Energies[4][4][4][4], double tstacki23Enthalpies[5][5][5][5], int NA);

void loadTstackh(double tstackhEnergies[4][4][4][4], double tstackhEnthalpies[5][5][5][5], int NA);

void loadTstackm(double tstackmEnergies[4][4][4][4], double tstackmEnthalpies[5][5][6][6], int NA, double saltCorrection);

void loadTstacke(double tstackeEnergies[4][4][4][4], double tstackeEnthalpies[5][5][6][6], int NA, double saltCorrection);

void combineTstack(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][5][5], double tRatio, ENERGY tstack[5][5][5][5]);

void combineTstack2(double tstackEnergies[4][4][4][4], double tstackEnthalpies[5][5][6][6], double tRatio, ENERGY tstack[5][5][6][6]);

void loadTstackiSuffix(ENERGY tstacki[5][5][5][5], char* suffix);

void loadTstacki23Suffix(ENERGY tstacki23[5][5][5][5], char* suffix);

void loadTstackhSuffix(ENERGY tstackh[5][5][5][5], char* suffix);

void loadTstackmSuffix(ENERGY tstackm[5][5][6][6], char* suffix);

void loadTstackeSuffix(ENERGY tstacke[5][5][6][6], char* suffix);

void loadCoaxialSuffix(ENERGY coaxial[5][5][5][5], char* suffix);

void loadTstackcoaxSuffix(ENERGY tstackcoax[5][5][5][5], char* suffix);

void loadCoaxstackSuffix(ENERGY coaxstack[5][5][5][5], char* suffix);

void loadMulti(double multiEnergies[3], double multiEnthalpies[3], int NA);

void loadMulti2(double multiEnergies[3], double multiEnthalpies[3], int NA);

void combineMulti(double multiEnergies[3], double multiEnthalpies[3], double tRatio, ENERGY multi[3]);

void calculateMulti(ENERGY multi[3], double tRatio, double scaleFactor);

void loadMultiSuffix(ENERGY multi[3], char* suffix);

void loadMulti2Suffix(ENERGY multi[3], char* suffix);

void loadMisc(double miscEnergies[13], double miscEnthalpies[13], int NA);

void combineMisc(double miscEnergies[13], double miscEnthalpies[13], double tRatio, ENERGY misc[13]);

void calculateMisc(ENERGY misc[13], double tRatio);

void loadMiscSuffix(ENERGY misc[13], char* suffix);

void makeAUPenalty(ENERGY misc[13], ENERGY aup[5][5], int isPF);

void makeAUPenaltyH(double misc[13], double aup[5][5], int isPF);

struct triloopE { char loop[5]; double energy; };
struct triloop { char loop[5]; ENERGY energy; };

void loadTriloop(struct triloopE**, struct triloopE**, int*, int);

void combineTriloop(const struct triloopE*, const struct triloopE*, double, struct triloop*, int);

void calculateTriloop(struct triloop*, int, double);

void loadTriloopSuffix(struct triloop**, int*, char*);

int triloopcmp(const void*, const void*);

struct tloopE { char loop[6]; double energy; };
struct tloop { char loop[6]; ENERGY energy; };

void loadTloop(struct tloopE**, struct tloopE**, int*, int);

void combineTloop(const struct tloopE*, const struct tloopE*, double, struct tloop*, int);

void calculateTloop(struct tloop*, int, double);

void loadTloopSuffix(struct tloop**, int*, char*);

int tloopcmp(const void*, const void*);

struct hexaloopE { char loop[8]; double energy; };
struct hexaloop { char loop[8]; ENERGY energy; };

void loadHexaloop(struct hexaloopE**, struct hexaloopE**, int*, int);

void combineHexaloop(const struct hexaloopE*, const struct hexaloopE*, double, struct hexaloop*, int);

void calculateHexaloop(struct hexaloop*, int, double);

void loadHexaloopSuffix(struct hexaloop**, int*, char*);

int hexaloopcmp(const void*, const void*);

void loadRTSuffix(double* RT, char* suffix);

#endif

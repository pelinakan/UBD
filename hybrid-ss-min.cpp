#include "hybrid-ss-min.h"
#include "util.h"
#include "CtEnergy.h"
#include "HybridMin.h"

double interiorLoopEnergies[30] = {999999,5,3.2,3.6,4,4.4,4.6,4.8,4.9,4.9,4.9,5,5.1,5.2,5.2,5.3,5.3,5.4,5.5,5.5,5.6,5.6,5.7,5.7,5.7,5.8,5.8,5.9,5.9,5.9};
double bulgeLoopEnergies[30] = {4,2.9,3.1,3.2,3.3,3.5,3.7,3.8,3.9,4,4.1,4.2,4.2,4.3,4.3,4.4,4.4,4.5,4.6,4.6,4.7,4.7,4.8,4.8,4.9,4.9,5,5.1,5.1,5.2};
double hairpinLoopEnergies[30] = {999999,999999,3.5,3.5,3.3,4,4.3,4.3,4.3,4.3,4.4,4.5,4.6,4.6,4.7,4.8,4.9,5,5,5.1,5.2,5.3,5.4,5.6,5.7,5.8,5.9,6,6.1,6.2};
double interiorLoopEnthalpies[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double bulgeLoopEnthalpies[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
double hairpinLoopEnthalpies[30] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double multiEnergies[3]={4.5,0.2,0.2};
double multiEnthalpies[3]={0,0,0};
double miscEnergies[13]={0.4,0.3,0.2,0.1,3.0,1.96,0.05,0.0,0.0,0.0,0.0,0.0,1.07857764};
double miscEnthalpies[13]={0.0,0.0,0.0,0.0,0.0,0.2 ,2.2 ,0.0,0.0,0.0,0.0,0.0,0.0};
/**
 * @brief Read all the stuff we need and get that done with!
 */

HybridSSMin::HybridSSMin()
{

	//Initialize ctEn computer
	enComputer = new CtEnergy();
	hybridComputer = new HybridMin();

	char gotSeq;
  int count, i, j;
  double t, tRatio;
  char *buffer;
  time_t now;
  int NA, polymer, constraints;
  //double tMin, tInc, tMax;
  double naConc, mgConc;
  double saltCorrection;

  q = qprime = qm = q5 = q3 = NULL;
#if ENABLE_FORCE
  g_ssok = NULL;
#endif
  //DNA!!!!!!!! NOT :) RNA
  NA = 1;
  gotSeq = 0;
  g_allPairs = 0;
  g_maxLoop = 30;
  g_nodangle = 0;
  g_simple = 0;
  g_prefix = NULL;
  naConc = 1.0;
  mgConc = 0.0;
  polymer = 0;
  g_prefilter1 = g_prefilter2 = 2;
  prohibitList = forceList = NULL;
  //g_mfoldMax = 0;
  g_mfoldP = 0;
  g_mfoldW = 0;
  g_maxBP = 0;
  g_circular = 0;
  g_quiet = 1;
  g_verbose = 0;
  /* initializations below are unnecessary but prevent compiler warnings */
  g_string = NULL;

	//	g_oneTemp = (tMin + tInc > tMax) ? 1 : 0;

	if (g_maxLoop < 0)
		g_maxLoop = 999999999;
	if (g_maxBP <= 0)
		g_maxBP = 999999999;

	saltCorrection = ion(NA, polymer, naConc, mgConc);

	/* read free energies and entropies */
  
    loadStack(stackEnergies, stackEnthalpies, NA, saltCorrection);
    symmetryCheckStack(stackEnergies, "energy");
    /* symmetryCheckStack(stackEnthalpies, "enthalpy"); */
    if (!g_nodangle)
		loadDangle(dangleEnergies3, dangleEnthalpies3, dangleEnergies5, dangleEnthalpies5, NA, saltCorrection);
    loadLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, NA, saltCorrection);
    loadSint2(sint2Energies, sint2Enthalpies, NA, saltCorrection);
    symmetryCheckSint2(sint2Energies, "energy");
    /* symmetryCheckSint2(sint2Enthalpies, "enthalpy"); */
    loadAsint1x2(asint1x2Energies, asint1x2Enthalpies, NA, saltCorrection);
    loadSint4(sint4Energies, sint4Enthalpies, NA, saltCorrection);
    symmetryCheckSint4(sint4Energies, "energy");
    /* symmetryCheckSint4(sint4Enthalpies, "enthalpy"); */
    loadTstackh(tstackhEnergies, tstackhEnthalpies, NA);
    loadTstacki(tstackiEnergies, tstackiEnthalpies, NA);
    loadTstacki23(tstacki23Energies, tstacki23Enthalpies, NA);
    if (!g_nodangle)
	{
		loadTstackm(tstackmEnergies, tstackmEnthalpies, NA, saltCorrection);
		loadTstacke(tstackeEnergies, tstackeEnthalpies, NA, saltCorrection);
	}
    loadTriloop(&triloopEnergies, &triloopEnthalpies, &numTriloops, NA);
    g_triloop = (struct triloop*) xcalloc(numTriloops, sizeof(struct triloop));
    loadTloop(&tloopEnergies, &tloopEnthalpies, &numTloops, NA);
    g_tloop = (struct tloop*) xcalloc(numTloops, sizeof(struct tloop));
    loadHexaloop(&hexaloopEnergies, &hexaloopEnthalpies, &numHexaloops, NA);
    g_hexaloop = (struct hexaloop*) xcalloc(numHexaloops, sizeof(struct hexaloop));
	//Never mind this one, we have them already!
	//TODO: Print out that this doesn't work for RNA!!!!!!!!!!!
    //loadMulti(multiEnergies, multiEnthalpies, NA);
    //loadMisc(miscEnergies, miscEnthalpies, NA);
	t=50.0;
	tRatio = (t + 273.15) / 310.15;
	RT = R * (t + 273.15);
	combineStack(stackEnergies, stackEnthalpies, tRatio, g_stack);
		if (!g_nodangle)
			combineDangle(dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, g_dangle3, g_dangle5);
		else
			calculateInfDangle(g_dangle3, g_dangle5);
		combineLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, tRatio, g_hairpinLoop, g_interiorLoop, g_bulgeLoop);
		combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
		combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
		combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
		combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
		combineTstack(tstacki23Energies, tstacki23Enthalpies, tRatio, g_tstacki23);
		combineTstack(tstackhEnergies, tstackhEnthalpies, tRatio, g_tstackh);
		if (g_nodangle)
		{
			calculateInfStack2(g_tstackm);
			calculateInfStack2(g_tstacke);
		}
		else
		{
			combineTstack2(tstackmEnergies, tstackmEnthalpies, tRatio, g_tstackm);
			combineTstack2(tstackeEnergies, tstackeEnthalpies, tRatio, g_tstacke);
		}
		combineMulti(multiEnergies, multiEnthalpies, tRatio, g_multi);
		combineMisc(miscEnergies, miscEnthalpies, tRatio, g_misc);
		combineTriloop(triloopEnergies, triloopEnthalpies, tRatio, g_triloop, numTriloops);
		combineTloop(tloopEnergies, tloopEnthalpies, tRatio, g_tloop, numTloops);
		combineHexaloop(hexaloopEnergies, hexaloopEnthalpies, tRatio, g_hexaloop, numHexaloops);

}

void HybridSSMin::computeTwoProbeHybridization(double& dG, double& dH, const char* seq1, const char* seq2, double temp)
{
	hybridComputer->compute(dG,dH,seq1,seq2,temp);
}

double HybridSSMin::computeGibsonFreeEnergy(double& gE, double& ctE, const char* sequence, double tMin=37.0, double tMax=37.0)
{

  int i,skipTraceback;
  skipTraceback=0;
  double t,tRatio;
  double tInc=1.0;
	  //Get sequence!!!!
	  g_string = sequence;
      g_len = strlen(g_string);
	  g_seq = NULL;
      /* convert sequence to numbers for speed */
      g_seq = (unsigned char*) xrealloc(g_seq, g_circular ? 2 * g_len + 2 : g_len + 2);
      for (i = 1; i <= g_len; ++i)
		g_seq[i] = util::toNum(g_string[i - 1]);
      if (g_circular)
	{
	  for (i = 1; i <= g_len; ++i)
	    g_seq[g_len + i] = util::toNum(g_string[i - 1]);
	  g_seq[0] = g_seq[2 * g_len + 1] = 5;
	}
      else
	g_seq[0] = g_seq[g_len + 1] = 5;

      if (g_circular)
	g_len *= 2;

      
	  q = recalloc2(q, g_len);
	  qprime = recalloc2(qprime, g_len);
	  qm = recalloc2(qm, g_len);

      q5 = (double*) xrealloc(q5, (g_len + 1) * sizeof(ENERGY));
      q3 = (double*) xrealloc(q3, (g_len + 1) * sizeof(ENERGY));
#if ENABLE_FORCE
      g_ssok = xrealloc(g_ssok, (g_len + 2) * (g_len + 2));
#endif

int *bp, *upst, *dnst;

for (t = tMin; t <= tMax; t += tInc)
{
	int bestI, bestJ;
	ENERGY mfe;
	
	tRatio = (t + 273.15) / 310.15;
	RT = R * (t + 273.15);
	
	makeAUPenalty(g_misc, g_aup, 0);
	if (g_simple)
		g_multi[1] = g_multi[2] = 0;
	
	bestI = bestJ = 0;
	
	initializeMatrices();
	fillMatrices1();
	
	computeQ53();

	if (fabs(Q5(g_len) - Q3(1)) > 1e-12)
		fprintf(stderr, "Warning: Q5(n) = %g but Q3(1) = %g.\n", (double) Q5(g_len) / PRECISION, (double) Q3(1) / PRECISION);
	
	mfe = Q5(g_len);
	
	int k;
	      
	bp = new int[g_len];
	upst = new int[g_len];
	dnst = new int[g_len];
	for (k = 0; k < g_len; ++k)
		bp[k] = upst[k] = dnst[k] = 0;

	
	if (g_noisolate && isFinite(Q5(g_len)))
		traceback_noI(0, 0, 0, bp, upst, dnst);
	else if (isFinite(Q5(g_len)))
		traceback(0, 0, 0, bp, upst, dnst);
		
}
  
	gE = (double) Q5(g_len) / PRECISION;
	ctE = enComputer->compute(this,bp,upst,dnst);

	delete[] bp;
	delete[] upst;
	delete[] dnst;

	return EXIT_SUCCESS;
}

int HybridSSMin::helixLength(int i, int j, ENERGY* qprime)
{
  int k, length;

  if (!isFinite(Qprime(i, j)))
    return 0;

  length = 1;
  for (k = 1; i + k < j - k && isFinite(Qprime(i + k, j - k)); ++k);
  length += k - 1;
  for (k = 1; i > k && j + k <= g_len && isFinite(Qprime(i - k, j + k)); ++k);
  length += k - 1;

  return length;
}

void HybridSSMin::prefilter()
{
  char** in;
  int i, j, k, count;

  in = (char**) xcalloc(g_len, sizeof(char*));
  for (i = 1; i <= g_len; ++i)
    in[i - 1] = (char*) xcalloc(g_len, 1);

  for (i = 1; i <= g_len - g_prefilter2 + 1; ++i)
    for (j = g_len; j >= g_prefilter2 && j >= i; --j)
      {
	count = 0;
	for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	  if (isFinite(Qprime(i + k, j - k)))
	    ++count;
	if (count >= g_prefilter1)
	  for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	    ++in[i + k - 1][j - k - 1];
      }

  for (i = 1; i <= g_len; ++i)
    {
      for (j = g_len; j >= i; --j)
	if (!in[i - 1][j - 1])
	  Qprime(i, j) = INFINITY;
      free(in[i - 1]);
    }
  free(in);
}

void HybridSSMin::initializeMatrices()
{
  int i, j, k;
  struct constraintListNode* top;

  /* Q' is initialized to +infinity iff base pair is illegal; 0 otherwise
     Q and QM are always initialized to +infinity */
  for (i = 1; i <= g_len; ++i)
    for (j = i; j <= g_len; ++j)
      if (j - i < TURN + 1 || (basePairIndex(g_seq[i], g_seq[j]) == 6 && !g_allPairs))
	Q(i, j) = Qprime(i, j) = QM(i, j) = INFINITY;
      else if (j - i > g_maxBP)
	Q(i, j) = Qprime(i, j) = QM(i, j) = INFINITY;
      else
	{
	  Q(i, j) = QM(i, j) = INFINITY;
	  Qprime(i, j) = 0;
	}

  top = prohibitList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len &&
	  top->k >= 1 && top->k <= g_len && top->l >= 1 && top->l <= g_len)
	for (i = top->i; i <= top->j; ++i)
	  for (j = top->k; j <= top->l; ++j)
	    {
	      if (i <= j)
		Qprime(i, j) = INFINITY;
	      else
		Qprime(j, i) = INFINITY;
	    }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len)
	for (k = 0; k < top->k; ++k)
	  {
	    if (top->i + k <= top->j - k)
	      Qprime(top->i + k, top->j - k) = INFINITY;
	    else
	      Qprime(top->j - k, top->i + k) = INFINITY;
	  }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j == 0)
	for (k = 0; k < top->k; ++k)
	  {
	    for (j = 1; j <= top->i + k; ++j)
	      Qprime(j, top->i + k) = INFINITY;
	    for (j = top->i + k; j <= g_len; ++j)
	      Qprime(top->i + k, j) = INFINITY;
	  }
      else if (top->l == 0 && top->j >= 1 && top->j <= g_len && top->i == 0)
	for (k = 0; k < top->k; ++k)
	  {
	    for (i = 1; i <= top->j + k; ++i)
	      Qprime(i, top->j + k) = INFINITY;
	    for (i = top->j + k; i <= g_len; ++i)
	      Qprime(top->j + k, i) = INFINITY;
	  }

      top = top->next;
    }

#if ENABLE_FORCE
  for (i = 0; i <= g_len + 1; ++i)
    for (j = 0; j <= g_len + 1; ++j)
      ssOK(i, j) = 1;
  top = forceList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len)
	for (i = 0; i <= g_len + 1; ++i)
	  for (j = i; j <= g_len + 1; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (i <= top->i + k && top->i + k <= j)
		ssOK(i, j) = 0;

      if (top->j >= 1 && top->j <= g_len)
	{
	  if (top->i == 0)
	    {
	      for (i = 0; i <= g_len + 1; ++i)
		for (j = i; j <= g_len + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j + k && top->j + k <= j)
		      ssOK(i, j) = 0;
	    }
	  else
	    {
	      for (i = 0; i <= g_len + 1; ++i)
		for (j = i; j <= g_len + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j - k && top->j - k <= j)
		      ssOK(i, j) = 0;
	    }
	}

      if (top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len)
	{
	  for (i = 1; i <= g_len; ++i)
	    for (k = 0; k < top->k; ++k)
	      if (i != top->i + k && i <= top->j - k)
		/* Qprime(i, top->j - k) = Qprime(top->j - k, i) = INFINITY; */
		Qprime(i, top->j - k) = INFINITY;
	  for (j = 1; j <= g_len; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (j != top->j - k && top->i + k <= j)
		/* Qprime(top->i + k, j) = Qprime(j, top->i + k) = INFINITY; */
		Qprime(top->i + k, j) = INFINITY;
	}
      top = top->next;
    }
#endif

  if (g_prefilter1 > 1 && !g_allPairs)
    prefilter();

    /* for (i = 1; i <= g_len; ++i) 
      for (j = i; j <= g_len; ++j)
	if (helixLength(i, j, qprime) <= g_prefilter)
	Qprime(i, j) = INFINITY; */

}

void HybridSSMin::fillMatrices1()
{
  int i, j, k;

  /* start at top left, fill each column bottom->top
     when Q' is +infinity, don't consider it */
  for (j = 2; j <= g_len; ++j)
    for (i = j - TURN - 1; i >= (g_circular && j > g_len / 2 ? j - g_len / 2: 1); --i)
      {
	ENERGY au;
	au = auPenalty(i, j);

	if (g_circular)
	  {
	    if (j - i > g_len / 2)
	      continue;
	    if (i > g_len / 2)
	      {
		Q(i, j) = Q(i - g_len / 2, j - g_len / 2);
		Qprime(i, j) = Qprime(i - g_len / 2, j - g_len / 2);
		QM(i, j) = QM(i - g_len / 2, j - g_len / 2);
		continue;
	      }
	  }

	if (isFinite(Qprime(i, j)))
	  {
	    Qprime(i, j) = min4(Eh(i, j),
				Es(i, j) + Qprime(i + 1, j - 1),
				QBI(i, j),
				g_multi[0] + g_multi[2] + au + QM(i + 1, j - 1));

	    if (!g_nodangle)
	      {
		if (j > 2)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed5(i, j) + QM(i + 1, j - 2));
		if (i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed3(i, j) + QM(i + 2, j - 1));
		if (j > 2 && i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + au + Etstackm(i, j) + QM(i + 2, j - 2));
	      }
	  }

	QM(i, j) = INFINITY;
	for (k = i + TURN + 1; k <= j - TURN - 2; ++k)
	  QM(i, j) = min2(QM(i, j), Q(i, k) + Q(k + 1, j));

	if (g_noisolate)
	  Q(i, j) = min4(ssOK(i, i) ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j, j) ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 i <= j - 2 ? g_multi[2] + au + Es(i, j) + Qprime(i + 1, j - 1) : INFINITY,
			 QM(i, j));
	else
	  Q(i, j) = min4(ssOK(i, i) ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j, j) ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 g_multi[2] + au + Qprime(i, j),
			 QM(i, j));
	  if (!g_nodangle)
	    {
	      if (g_noisolate)
		{
		  if (i < j - TURN - 3)
		    {
		      Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1));
		      Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Es(i, j - 1) + Qprime(i + 1, j - 2));
		    }
		  if (i < j - TURN - 4)
		    Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2));
		}
	      else
		{
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Qprime(i + 1, j));
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Qprime(i, j - 1));
		  if (i < j - TURN - 2)
		    Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Qprime(i + 1, j - 1));
		}
	    }
      }
}

void HybridSSMin::fillMatrices2()
{
  int i, j, k;
  FILE* file;

  for (j = g_len + 1; j < 2 * g_len; ++j)
    for (i = g_len; i > j - g_len; --i)
      {
	ENERGY au;
	au = auPenalty(i, j - g_len);

	if (isFinite(Qprime(i, j)))
	  {
	    Qprime(i, j) = min4(i < g_len ? Es(i, j) + Qprime(i + 1, j - 1) : INFINITY,
				QBI2(i, j),
				i < g_len ? g_multi[0] + g_multi[2] + au + QM(i + 1, j - 1) : INFINITY,
				au + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY));
	    if (!g_nodangle)
	      {
		if (i < g_len && j > g_len + 2)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed5(i, j - g_len) + QM(i + 1, j - 2));
		if (i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed3(i, j - g_len) + QM(i + 2, j - 1));
		if (j > g_len + 2 && i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + au + Etstackm(i, j - g_len) + QM(i + 2, j - 2));
		if (j > g_len + 1)
		  Qprime(i, j) = min2(Qprime(i, j), au + Ed5(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY));
		if (i < g_len)
		  Qprime(i, j) = min2(Qprime(i, j), au + Ed3(i, j - g_len) + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY));
		if (j > g_len + 1 && i < g_len)
		  Qprime(i, j) = min2(Qprime(i, j), au + Etstacke(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY));
	      }
	  }

	QM(i, j) = INFINITY;
	for (k = i + TURN + 1; k <= g_len - 1; ++k)
	  QM(i, j) = min2(QM(i, j), Q(i, k) + Q(k + 1, j));
	for (k = g_len + 1; k <= j - TURN - 2; ++k)
	  QM(i, j) = min2(QM(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len));

	if (g_noisolate)
	  Q(i, j) = min4(ssOK(i, i) && i < g_len ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j - g_len, j - g_len) && j > g_len + 1 ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 i < g_len && j > g_len + 1 ? g_multi[2] + au + Es(i, j) + Qprime(i + 1, j - 1) : INFINITY,
			 QM(i, j));
	else
	  Q(i, j) = min4(ssOK(i, i) && i < g_len ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j - g_len, j - g_len) && j > g_len + 1 ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 g_multi[2] + au + Qprime(i, j),
		       QM(i, j));
	if (!g_nodangle)
	  {
	    if (g_noisolate)
	      {
		if (i < g_len - 1 && j > g_len + 1)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1));
		if (i < g_len && j > g_len + 2)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Es(i, j - 1) + Qprime(i + 1, j - 2));
		if (i < g_len - 1 && j > g_len + 2)
		  Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2));
	      }
	    else
	      {
		if (i < g_len)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Qprime(i + 1, j));
		if (j > g_len + 1)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Qprime(i, j - 1));
		if (i < g_len && j > g_len + 1)
		  Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Qprime(i + 1, j - 1));
	      }
	  }
      }

  if (g_debug)
    {
      file = fopen("Qprime-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%f\t", (double) Qprime(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%f\t", (double) Q(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("QM-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%f\t", (double) QM(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void HybridSSMin::computeQ53()
{
  int i, j;

  Q5(0) = Q5(1) = INFINITY;
  Q3(g_len + 1) = Q3(g_len) = INFINITY;

  if (g_nodangle)
    {
      for (i = 2; i <= g_len; ++i)
	Q5(i) = min2(ssOK(i, i) ? Q5(i - 1) : INFINITY, Q5_1(i));

      for (j = g_len - 1; j >= 1; --j)
	Q3(j) = min2(ssOK(j, j) ? Q3(j + 1) : INFINITY, Q3_1(j));
    }
  else
    {
      for (i = 2; i <= g_len; ++i)
	Q5(i) = min5(ssOK(i, i) ? Q5(i - 1) : INFINITY, Q5_1(i), Q5_2(i), Q5_3(i), Q5_4(i));

      for (j = g_len - 1; j >= 1; --j)
	Q3(j) = min5(ssOK(j, j) ? Q3(j + 1) : INFINITY, Q3_1(j), Q3_2(j), Q3_3(j), Q3_4(j));
    }
}

void HybridSSMin::traceback(int i, int j, int force, int* bp, int* upst, int* dnst)
{
  int k, ii, jj;
  struct stackNode *top, *stack = NULL;

  /* find best folding from i to j; force i-j pair iff force
     if i=j=0, best folding including Q5 and Q3 is found */

  if (i == 0 && j == 0)
    push(&stack, g_len, 0, 3);
  else
    {
      if (force)
	push(&stack, i, j, j > g_len ? 5 : 0);
      else
	push(&stack, i, j, 2);
    }

  while (stack)
    {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;

      if (top->matrix == 1) /* QM */
	{
	  for (k = i + TURN + 1; k < j - TURN - 1; ++k)
	    if (equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 2);
		break;
	      }
	}
      else if (top->matrix == 2) /* Q */
	{
	  while ((ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j) + Qprime(i, j)))
	    push(&stack, i, j, 0);
	  else if (equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Qprime(i + 1, j)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 1, j, 0);
	    }
	  else if (equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Qprime(i, j - 1)))
	    {
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i, j - 1, 0);
	    }
	  else if (equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Qprime(i + 1, j - 1)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else
	    {
	      for (k = i + TURN + 1; k <= j - TURN - 2; ++k)
		if (equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 2);
		    break;
		  }
	    }
	}
      else if (top->matrix == 0) /* Q' */
	{
	  bp[i - 1] = j;
	  bp[j - 1] = i;
	  if (equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else if (equal(Qprime(i, j), Eh(i, j)))
	    ;
	  else if (equal(Qprime(i, j), QBI(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = i + 1; ii < j - d; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(i, j, ii, jj) + Qprime(ii, jj)))
		      {
			setBI(i, j, ii, jj, upst, dnst);
			push(&stack, ii, jj, 0);
			++done;
			break;
		      }		    
		  }
	    }
	  else if (i <= j - 2 && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 1);
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed5(i, j) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j, upst, dnst);
	      push(&stack, i + 1, j - 2, 1);
	    }
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed3(i, j) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 1);
	    }
	  else
#ifdef DEBUG
	  if (i <= j - 4 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j) + Etstackm(i, j) + QM(i + 2, j - 2)))
#endif
	    {
	      setDangle5(j, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 1);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
      else if (top->matrix == 3) /* Q5 */
	{
	  while (ssOK(i, i) && equal(Q5(i), Q5(i - 1)))
	    --i;

	  if (i == 0)
	    continue;

	  if (equal(Q5(i), Q5_1(i)))
	    {
	      for (k = 0; k <= i - TURN - 2; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i) + Qprime(k + 1, i)))
		  {
		    push(&stack, k + 1, i, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i) + Qprime(k + 1, i)))
		  {
		    push(&stack, k + 1, i, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_2(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i) + Ed5(i, k + 2) + Qprime(k + 2, i)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 2, i, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Qprime(k + 2, i)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 2, i, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_3(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Qprime(k + 1, i - 1)))
		  {
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 1, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Qprime(k + 1, i - 1)))
		  {
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 1, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q5(i), Q5_4(i)))
#endif
	    {
	      for (k = 0; k <= i - TURN - 4; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Qprime(k + 2, i - 1)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Qprime(k + 2, i - 1)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q5(%d)\n", i);
#endif
	}
      else if (top->matrix == 4) /* Q3 */
	{
	  while (ssOK(j, j) && equal(Q3(j), Q3(j + 1)))
	    ++j;

	  if (j == g_len + 1)
	    continue;

	  if (equal(Q3(j), Q3_1(j)))
	    {
	      for (k = j + TURN + 2; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 1) + Qprime(j, k - 1)))
		  {
		    push(&stack, j, k - 1, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 1) + Qprime(j, k - 1) + Q3(k)))
		  {
		    push(&stack, j, k - 1, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_2(j)))
	    {
	      for (k = j + TURN + 3; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Qprime(j, k - 2)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Qprime(j, k - 2) + Q3(k)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_3(j)))
	    {
	      for (k = j + TURN + 3; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Qprime(j + 1, k - 1)))
		  {
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 1, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Qprime(j + 1, k - 1) + Q3(k)))
		  {
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 1, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q3(j), Q3_4(j)))
#endif
	    {
	      for (k = j + TURN + 4; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Qprime(j + 1, k - 2)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Qprime(j + 1, k - 2) + Q3(k)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q3(%d)\n", j);
#endif
	}
      else if (top->matrix == 6) /* QM */
	{
	  for (k = i; k <= j - 1; ++k)
	    if (k < g_len && equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 7);
		break;
	      }
	    else if (k > g_len && equal(QM(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
	      {
		push(&stack, i, k, 7);
		push(&stack, k + 1 - g_len, j - g_len, 2);
		break;
	      }
	}
      else if (top->matrix == 7) /* Q */
	{
	  while ((i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j - g_len) + Qprime(i, j)))
	    push(&stack, i, j, 5);
	  else if (i < g_len && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Qprime(i + 1, j)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 1, j, 5);
	    }
	  else if (j > g_len + 1 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Qprime(i, j - 1)))
	    {
	      setDangle3(j - 1 - g_len, upst, dnst);
	      push(&stack, i, j - 1, 5);
	    }
	  else if (i < g_len && j > g_len + 1 && equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Qprime(i + 1, j - 1)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1 - g_len, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else
	    {
	      for (k = i; k <= j - 1; ++k)
		if (k < g_len && equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 7);
		    break;
		  }
		else if (k > g_len && equal(Q(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
		  {
		    push(&stack, i, k, 7);
		    push(&stack, k + 1 - g_len, j - g_len, 2);
		    break;
		  }
	    }
	}
      else /* Q' */
#ifdef DEBUG
      if (top->matrix == 5)
#endif
	{
	  bp[(i - 1) % g_len] = j - g_len;
	  bp[j - 1 - g_len] = i > g_len ? i - g_len : i;
	  if (i < g_len && equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else if (equal(Qprime(i, j), QBI2(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= 1 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii <= g_len; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(jj - g_len, ii, j - g_len, i) + Qprime(ii, jj)))
		      {
			setBI(i, j - g_len, ii, jj - g_len, upst, dnst);
			push(&stack, ii, jj, 5);
			++done;
			break;
		      }		    
		  }
	    }
	  else if (i < g_len && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j - g_len) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 6);
	  else if (i < g_len && j > g_len + 2 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed5(i, j - g_len) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 2, 6);
	    }
	  else if (i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed3(i, j - g_len) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 6);
	    }
	  else if (j > g_len + 2 && i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Etstackm(i, j - g_len) + QM(i + 2, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 6);
	    }
	  else if (equal(Qprime(i, j), auPenalty(i, j - g_len) +
			 min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) +
			 min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || Q3(i + 1) < 0.0)
		push(&stack, 0, i + 1, 4);
	    }
	  else if (j > g_len + 1 && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed5(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || (i < g_len && Q3(i + 1) < 0.0))
		push(&stack, 0, i + 1, 4);
	    }
	  else if (i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed3(i, j - g_len) + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
	  else
#ifdef DEBUG
	  if (j > g_len + 1 && i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Etstacke(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
#endif
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
#ifdef DEBUG
      else
	fputs("Error in traceback\n", stderr);
#endif
      free(top);
    }
}

void HybridSSMin::traceback_noI(int i, int j, int force, int* bp, int* upst, int* dnst)
{
  int k, ii, jj;
  struct stackNode *top, *stack = NULL;

  /* find best folding from i to j; force i-j pair iff force
     if i=j=0, best folding including Q5 and Q3 is found */

  if (i == 0 && j == 0)
    push(&stack, g_len, 0, 3);
  else
    {
      if (force)
	push(&stack, i, j, j > g_len ? 5 : 0);
      else
	push(&stack, i, j, 2);
    }

  while (stack)
    {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;

      if (top->matrix == 1) /* QM */
	{
	  for (k = i + TURN + 1; k < j - TURN - 1; ++k)
	    if (equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 2);
		break;
	      }
	}
      else if (top->matrix == 2) /* Q */
	{
	  while ((ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j) + Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      bp[i - 1] = j;
	      bp[j - 1] = i;
	      setStack(i, j, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else if (i < g_len - 1 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1)))
	    {
	      bp[i] = j;
	      bp[j - 1] = i + 1;
	      setStack(i + 1, j, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 2, j - 1, 0);
	    }
	  else if (j > 2 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Es(i, j - 1) + Qprime(i + 1, j - 2)))
	    {
	      bp[i - 1] = j - 1;
	      bp[j - 2] = i;
	      setStack(i, j - 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 1, j - 2, 0);
	    }
	  else if (i < g_len - 1 && j > 2 && equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2)))
	    {
	      bp[i] = j - 1;
	      bp[j - 2] = i + 1;
	      setStack(i + 1, j - 1, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 2, j - 2, 0);
	    }
	  else
	    {
	      for (k = i; k <= j - 1; ++k)
		if (equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 2);
		    break;
		  }
	    }
	}
      else if (top->matrix == 0) /* Q' */
	{
	  bp[i - 1] = j;
	  bp[j - 1] = i;
	  if (equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else if (equal(Qprime(i, j), Eh(i, j)))
	    ;
	  else if (equal(Qprime(i, j), QBI(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= TURN + 3 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = i + 1; ii < j - d; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(i, j, ii, jj) + Es(ii, jj) + Qprime(ii + 1, jj - 1)))
		      {
			bp[ii - 1] = jj;
			bp[jj - 1] = ii;
			setBI(i, j, ii, jj, upst, dnst);
			setStack(ii, jj, upst, dnst);
			push(&stack, ii + 1, jj - 1, 0);
			++done;
			break;
		      }
		  }
	    }
	  else if (i <= j - 2 && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 1);
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed5(i, j) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j, upst, dnst);
	      push(&stack, i + 1, j - 2, 1);
	    }
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed3(i, j) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 1);
	    }
	  else
#ifdef DEBUG
	  if (i <= j - 4 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j) + Etstackm(i, j) + QM(i + 2, j - 2)))
#endif
	    {
	      setDangle5(j, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 1);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
      else if (top->matrix == 3) /* Q5 */
	{
	  while (ssOK(i, i) && equal(Q5(i), Q5(i - 1)))
	    --i;

	  if (i == 0)
	    continue;

	  if (equal(Q5(i), Q5_1(i)))
	    {
	      for (k = 0; k <= i - TURN - 2; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i) + Es(k + 1, i) + Qprime(k + 2, i - 1)))
		  {
		    bp[k] = i;
		    bp[i - 1] = k + 1;
		    setStack(k + 1, i, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i) + Es(k + 1, i) + Qprime(k + 2, i - 1)))
		  {
		    bp[k] = i;
		    bp[i - 1] = k + 1;
		    setStack(k + 1, i, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_2(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i) + Ed5(i, k + 2) + Es(k + 2, i) + Qprime(k + 3, i - 1)))
		  {
		    bp[k + 1] = i;
		    bp[i - 1] = k + 2;
		    setStack(k + 2, i, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 3, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Es(k + 2, i) + Qprime(k + 3, i - 1)))
		  {
		    bp[k + 1] = i;
		    bp[i - 1] = k + 2;
		    setStack(k + 2, i, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 3, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_3(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Es(k + 1, i - 1) + Qprime(k + 2, i - 2)))
		  {
		    bp[k] = i - 1;
		    bp[i - 2] = k + 1;
		    setStack(k + 1, i - 1, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 2, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Es(k + 1, i - 1) + Qprime(k + 2, i - 2)))
		  {
		    bp[k] = i - 1;
		    bp[i - 2] = k + 1;
		    setStack(k + 1, i - 1, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 2, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q5(i), Q5_4(i)))
#endif
	    {
	      for (k = 0; k <= i - TURN - 4; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Es(k + 2, i - 1) + Qprime(k + 3, i - 2)))
		  {
		    bp[k + 1] = i - 1;
		    bp[i - 2] = k + 2;
		    setStack(k + 2, i - 1, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 3, i - 2, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Es(k + 2, i - 1) + Qprime(k + 3, i - 2)))
		  {
		    bp[k + 1] = i - 1;
		    bp[i - 2] = k + 2;
		    setStack(k + 2, i - 1, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 3, i - 2, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q5(%d)\n", i);
#endif
	}
      else if (top->matrix == 4) /* Q3 */
	{
	  while (ssOK(j, j) && equal(Q3(j), Q3(j + 1)))
	    ++j;

	  if (j == g_len + 1)
	    continue;

	  if (equal(Q3(j), Q3_1(j)))
	    {
	      for (k = j + TURN + 4; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 1) + Es(j, k - 1) + Qprime(j + 1, k - 2)))
		  {
		    bp[k - 2] = j;
		    bp[j - 1] = k - 1;
		    setStack(j, k - 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 1) + Es(j, k - 1) + Qprime(j + 1, k - 2) + Q3(k)))
		  {
		    bp[k - 2] = j;
		    bp[j - 1] = k - 1;
		    setStack(j, k - 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_2(j)))
	    {
	      for (k = j + TURN + 5; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Es(j, k - 2) + Qprime(j + 1, k - 3)))
		  {
		    bp[k - 3] = j;
		    bp[j - 1] = k - 2;
		    setStack(j, k - 2, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 1, k - 3, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Es(j, k - 2) + Qprime(j + 1, k - 3) + Q3(k)))
		  {
		    bp[k - 3] = j;
		    bp[j - 1] = k - 2;
		    setStack(j, k - 2, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 1, k - 3, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_3(j)))
	    {
	      for (k = j + TURN + 5; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Es(j + 1, k - 1) + Qprime(j + 2, k - 2)))
		  {
		    bp[k - 2] = j + 1;
		    bp[j] = k - 1;
		    setStack(j + 1, k - 1, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 2, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Es(j + 1, k - 1) + Qprime(j + 2, k - 2) + Q3(k)))
		  {
		    bp[k - 2] = j + 1;
		    bp[j] = k - 1;
		    setStack(j + 1, k - 1, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 2, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q3(j), Q3_4(j)))
#endif
	    {
	      for (k = j + TURN + 6; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Es(j + 1, k - 2) + Qprime(j + 2, k - 3)))
		  {
		    bp[k - 3] = j + 1;
		    bp[j] = k - 2;
		    setStack(j + 1, k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 2, k - 3, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Es(j + 1, k - 2) + Qprime(j + 2, k - 3) + Q3(k)))
		  {
		    bp[k - 3] = j + 1;
		    bp[j] = k - 2;
		    setStack(j + 1, k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 2, k - 3, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q3(%d)\n", j);
#endif
	}
      else if (top->matrix == 6) /* QM */
	{
	  for (k = i; k <= j - 1; ++k)
	    if (k < g_len && equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 7);
		break;
	      }
	    else if (k > g_len && equal(QM(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
	      {
		push(&stack, i, k, 7);
		push(&stack, k + 1 - g_len, j - g_len, 2);
		break;
	      }
	}
      else if (top->matrix == 7) /* Q */
	{
	  while ((i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j - g_len) + Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      bp[i - 1] = j - g_len;
	      bp[j - 1 - g_len] = i;
	      setStack(i, j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else if (i < g_len && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1)))
	    {
	      bp[i] = j - g_len;
	      bp[j - 1 - g_len] = i + 1;
	      setStack(i + 1, j - g_len, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 2, j - 1, 5);
	    }
	  else if (j > g_len + 1 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Es(i, j - 1) + Qprime(i + 1, j - 2)))
	    {
	      bp[i - 1] = j - 1 - g_len;
	      bp[j - 2 - g_len] = i;
	      setStack(i, j - 1 - g_len, upst, dnst);
	      setDangle3(j - 1 - g_len, upst, dnst);
	      push(&stack, i + 1, j - 2, 5);
	    }
	  else if (i < g_len && j > g_len + 1 && equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2)))
	    {
	      bp[i] = j - 1 - g_len;
	      bp[j - 2 - g_len] = i + 1;
	      setStack(i + 1, j - 1 - g_len, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1 - g_len, upst, dnst);
	      push(&stack, i + 2, j - 2, 5);
	    }
	  else
	    {
	      for (k = i; k <= j - 1; ++k)
		if (k < g_len && equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 7);
		    break;
		  }
		else if (k > g_len && equal(Q(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
		  {
		    push(&stack, i, k, 7);
		    push(&stack, k + 1 - g_len, j - g_len, 2);
		    break;
		  }
	    }
	}
      else /* Q' */
#ifdef DEBUG
      if (top->matrix == 5)
#endif
	{
	  bp[(i - 1) % g_len] = j - g_len;
	  bp[j - 1 - g_len] = i > g_len ? i - g_len : i;
	  if (i < g_len && equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else if (equal(Qprime(i, j), QBI2(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= 2 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii < g_len; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(jj - g_len, ii, j - g_len, i) + Es(ii, jj) + Qprime(ii + 1, jj - 1)))
		      {
			bp[ii - 1] = jj - g_len;
			bp[jj - 1 - g_len] = ii;
			setBI(i, j - g_len, ii, jj - g_len, upst, dnst);
			setStack(ii, jj - g_len, upst, dnst);
			push(&stack, ii + 1, jj - 1, 5);
			++done;
			break;
		      }		    
		  }
	    }
	  else if (i < g_len && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j - g_len) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 6);
	  else if (i < g_len && j > g_len + 2 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed5(i, j - g_len) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 2, 6);
	    }
	  else if (i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed3(i, j - g_len) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 6);
	    }
	  else if (j > g_len + 2 && i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Etstackm(i, j - g_len) + QM(i + 2, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 6);
	    }
	  else if (equal(Qprime(i, j), auPenalty(i, j - g_len) +
			 min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) +
			 min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || Q3(i + 1) < 0.0)
		push(&stack, 0, i + 1, 4);
	    }
	  else if (j > g_len + 1 && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed5(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || (i < g_len && Q3(i + 1) < 0.0))
		push(&stack, 0, i + 1, 4);
	    }
	  else if (i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed3(i, j - g_len) + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
	  else
#ifdef DEBUG
	  if (j > g_len + 1 && i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Etstacke(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
#endif
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d) = %g\n", i, j, (double) Qprime(i, j) / PRECISION);
#endif
	}
#ifdef DEBUG
      else
	fputs("Error in traceback\n", stderr);
#endif
      free(top);
    }
}

void HybridSSMin::setStack(int i, int j, int* upst, int* dnst)
{
  upst[i] = i;
  dnst[i - 1] = i + 1;
  upst[j - 1] = j - 1;
  dnst[j - 2] = j;
}

void HybridSSMin::setDangle5(int j, int* upst, int* dnst)
{
  upst[j - 1] = j - 1;
  dnst[j - 2] = j;
}

void HybridSSMin::setDangle3(int i, int* upst, int* dnst)
{
  upst[i] = i;
  dnst[i - 1] = i + 1;
}

void HybridSSMin::setBI(int i, int j, int ii, int jj, int* upst, int* dnst)
{
  int loopSize1, loopSize2;

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;

#ifdef DEBUG
  if (loopSize1 < 0 || loopSize2 < 0 || (loopSize1 == 0 && loopSize2 == 0))
    {
      fputs("Error: setBI() called with nonsense\n", stderr);
      return;
    }
  else
#endif

  if ((loopSize1 == 0 && loopSize2 == 1) || (loopSize2 == 0 && loopSize1 == 1))
    {
      upst[ii - 1] = i;
      dnst[i - 1] = ii;
      upst[j - 1] = jj;
      dnst[jj - 1] = j;
    }
  else if (loopSize1 && loopSize2 && (loopSize1 > 2 || loopSize2 > 2))
    {
      upst[i] = i;
      upst[ii - 1] = ii - 1;
      dnst[i - 1] = i + 1;
      dnst[ii - 2] = ii;
      upst[jj] = jj;
      upst[j - 1] = j - 1;
      dnst[jj - 1] = jj + 1;
      dnst[j - 2] = j;
    }
}

void HybridSSMin::circularize(int* bp, int* upst, int* dnst, int bestI)
{
  int i, n;

  n = g_len / 2;
  for (i = 1; i < bestI; ++i)
    {
      bp[i - 1] = bp[n + i - 1];
      if (bp[i - 1] > n)
	bp[i - 1] -= n;
      upst[i - 1] = upst[n + i - 1];
      if (upst[i - 1] > n)
	upst[i - 1] -= n;
      dnst[i - 1] = dnst[n + i - 1];
      if (dnst[i - 1] > n)
	dnst[i - 1] -= n;
    }
  upst[bestI - 1] = upst[n + bestI - 1];
  for (i = bestI; i <= n; ++i)
    {
      if (bp[i - 1] > n)
	bp[i - 1] -= n;
      if (upst[i - 1] > n)
	upst[i - 1] -= n;
      if (dnst[i - 1] > n)
	dnst[i - 1] -= n;
    }
}

int HybridSSMin::unique(int* bp, char** found)
{
  int i, unique_;

  unique_ = 0;
  for (i = 1; i <= (g_circular ? g_len / 2 : g_len); ++i)
    if (bp[i - 1] > i && !found[i - 1][bp[i - 1] - 1])
      ++unique_;

  return unique_;
}

void HybridSSMin::makePairList(ENERGY cutoff, int* pnum)
{
  int d, i, j, length, n;
  ENERGY E;

  n = g_circular ? g_len / 2 : g_len;

  pairList = NULL;
  for (d = 2; d < 2 * n - 1; ++d)
    {
      i = d / 2;
      j = (d + 1) / 2 + 1;
      length = 0;
      E = INFINITY;
      for (; i >= 1 && j <= n; --i, ++j)
	{
	  if (Qprime(i, j) + Qprime(j, i + n) <= cutoff)
	    {
	      ++pnum[i - 1];
	      ++pnum[j - 1];
	    }

	  if (length && (!isFinite(Qprime(i, j)) || !isFinite(Qprime(j, i + n))) && !isFinite(E))
	    ++length;
	  else if (length && equal(Qprime(i, j) + Qprime(j, i + n), E))
	    ++length;
	  else if (isFinite(Qprime(i, j)))
	    {
	      if (length && E <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 1;
	      E = Qprime(i, j) + Qprime(j, i + n);
	    }
	  else if (length)
	    {
	      if (E <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 0;
	    }
	  if ((i == 1 || j == n) && length)
	    if (E <= cutoff)
	      pushPairList(i, j, length, E);
	}
    }

}

void HybridSSMin::makePairList_noI(ENERGY cutoff, int* pnum)
{
  int d, i, j, length, n;
  ENERGY E;

  n = g_circular ? g_len / 2 : g_len;

  pairList = NULL;
  for (d = 4; d < 2 * n - 3; ++d)
    {
      i = d / 2;
      j = (d + 1) / 2 + 1;
      length = 0;
      E = INFINITY;
      for (; i >= 2 && j <= n - 1; --i, ++j)
	{
	  if (Qprime(i, j) + Es(i - 1, j + 1) + Qprime(j + 1, i - 1 + n) <= cutoff)
	    {
	      ++pnum[i - 1];
	      ++pnum[j - 1];
	      ++pnum[i - 2];
	      ++pnum[j];
	    }

	  if (length && equal(Qprime(i, j) + Es(i - 1, j + 1) + Qprime(j + 1, i - 1 + n), E))
	    ++length;
	  else if (isFinite(Qprime(i, j)) && isFinite(Es(i - 1, j + 1)) && isFinite(Qprime(j + 1, i - 1 + n)))
	    {
	      if (length && E <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 1;
	      E = Qprime(i, j) + Es(i - 1, j + 1) + Qprime(j + 1, i - 1 + n);
	    }
	  else if (length)
	    {
	      if (E <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 0;
	    }
	  if ((i == 2 || j == n - 1) && length)
	    if (E <= cutoff)
	      pushPairList(i, j, length + 1, E);
	}
    }

}

ENERGY HybridSSMin::Q5_1(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 4; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i) + Es(k + 1, i) + Qprime(k + 2, i - 1));
  else
    for (k = 0; k <= i - TURN - 2; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i) + Qprime(k + 1, i));

  return min;
}

ENERGY HybridSSMin::Q5_2(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 5; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Es(k + 2, i) + Qprime(k + 3, i - 1));
  else
    for (k = 0; k <= i - TURN - 3; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Qprime(k + 2, i));

  return min;
}

ENERGY HybridSSMin::Q5_3(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 5; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Es(k + 1, i - 1) + Qprime(k + 2, i - 2));
  else
    for (k = 0; k <= i - TURN - 3; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Qprime(k + 1, i - 1));
  
  return min;
}

ENERGY HybridSSMin::Q5_4(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 6; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Es(k + 2, i - 1) + Qprime(k + 3, i - 2));
  else
    for (k = 0; k <= i - TURN - 4; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Qprime(k + 2, i - 1));

  return min;
}

ENERGY HybridSSMin::Q3_1(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 4; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 1) + Es(j, k - 1) + Qprime(j + 1, k - 2) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 2; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 1) + Qprime(j, k - 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY HybridSSMin::Q3_2(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 5; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 2) + Es(j, k - 2) + Qprime(j + 1, k - 3) + Ed3(k - 2, j) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 3; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 2) + Qprime(j, k - 2) + Ed3(k - 2, j) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY HybridSSMin::Q3_3(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 5; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 1) + Es(j + 1, k - 1) + Qprime(j + 2, k - 2) + Ed5(k - 1, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 3; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 1) + Qprime(j + 1, k - 1) + Ed5(k - 1, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY HybridSSMin::Q3_4(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 6; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 2) + Es(j + 1, k - 2) + Qprime(j + 2, k - 3) + Etstacke(k - 2, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 4; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 2) + Qprime(j + 1, k - 2) + Etstacke(k - 2, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY HybridSSMin::Ed3(int i, int j)
{
  return ssOK(i + 1, i + 1) ? g_dangle3[g_seq[i]][g_seq[j]][g_seq[i + 1]] : INFINITY;
}

ENERGY HybridSSMin::Ed5(int i, int j)
{
  return ssOK(j - 1, j - 1) ? g_dangle5[g_seq[i]][g_seq[j]][g_seq[j - 1]] : INFINITY;
}

ENERGY HybridSSMin::Etstackm(int i, int j)
{
  return (ssOK(i + 1, i + 1) && ssOK(j - 1, j - 1)) ? g_tstackm[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] : INFINITY;
}

ENERGY HybridSSMin::Etstacke(int i, int j)
{
  return (ssOK(i + 1, i + 1) && ssOK(j - 1, j - 1)) ? g_tstacke[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] : INFINITY;
}

ENERGY HybridSSMin::Eh(int i, int j)
{
  ENERGY energy;
  int loopSize = j - i - 1;
  int k;

  if (loopSize < TURN)
    return INFINITY;

  if (i <= g_len && g_len < j)
    return INFINITY;
  else if (i > g_len)
    {
      i -= g_len;
      j -= g_len;
    }

#if ENABLE_FORCE
  if (!ssOK(i + 1, j - 1))
    return INFINITY;
#endif

  if (loopSize <= 30)
    energy = g_hairpinLoop[loopSize - 1];
  else
    energy = g_hairpinLoop[29] + g_misc[12] * log((double) loopSize / 30);

  if (loopSize > 3)
    energy += g_tstackh[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
  else
    energy += auPenalty(i, j);

  if (loopSize == 3)
    {
      struct triloop* loop;
      if (numTriloops)
	if ((loop = (triloop*) bsearch(g_seq + i, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
	  energy += loop->energy;
    }
  else if (loopSize == 4)
    {
      struct tloop* loop;
      if (numTloops)
	if ((loop = (tloop*) bsearch(g_seq + i, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
	  energy += loop->energy;
    }
  else if (loopSize == 6)
    {
      struct hexaloop* loop;
      if (numHexaloops)
	if ((loop = (hexaloop*) bsearch(g_seq + i, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
	  energy += loop->energy;
    }

  /* GGG */
  if (i >= 3 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[i] == 2 && g_seq[j] == 3)
    energy += g_misc[8];

  /* poly-C */
  if (loopSize == 3 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1 && g_seq[i + 3] == 1)
    energy += g_misc[11];
  else
    {
      for (k = 1; k <= loopSize; ++k)
	if (g_seq[i + k] != 1)
	  return energy;
      energy += g_misc[9] * loopSize + g_misc[10];
    }

  return energy;
}

ENERGY HybridSSMin::Es(int i, int j)
{
  if (i >= j)
    return INFINITY;
  /* fputs("Error in Es(): i isn't less than j\n", stderr); */

  if (i == g_len || j == g_len + 1)
    return INFINITY;

  if (i > g_len)
    i -= g_len;
  if (j > g_len)
    j -= g_len;

  return g_stack[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
}

ENERGY HybridSSMin::Ebi(int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  ENERGY loopEnergy, asPenalty;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
  if (ii >= jj)
    fputs("Error in Ebi(): jj isn't greater than ii\n", stderr);

  if ((i <= g_len && g_len < ii) || (jj <= g_len && g_len < j))
    return INFINITY;
#endif

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > g_maxLoop)
    return INFINITY;

#ifdef DEBUG
  if (i > g_len)
    i -= g_len;
  if (ii > g_len)
    ii -= g_len;
  if (j > g_len)
    j -= g_len;
  if (jj > g_len)
    jj -= g_len;
#endif

#if ENABLE_FORCE
  if (loopSize1 && !ssOK(i + 1, ii - 1))
    return INFINITY;
  if (loopSize2 && !ssOK(jj + 1, j - 1))
    return INFINITY;
#endif

#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Ebi() called with nonsense\n", stderr);
      return INFINITY;
    }
  else
#endif
  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return g_bulgeLoop[0] + g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]];
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] + auPenalty(i, j) + auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + auPenalty(i, j) + auPenalty(ii, jj);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] + g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]];
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] + auPenalty(i, j) + auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + auPenalty(i, j) + auPenalty(ii, jj);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(g_seq[jj], g_seq[ii])][basePairIndex(g_seq[j], g_seq[i])][g_seq[jj + 1]][g_seq[ii - 1]][g_seq[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[i + 2]][g_seq[j - 2]];
  else if ((loopSize1 == 2 && loopSize2 == 3) ||
	   (loopSize1 == 3 && loopSize2 == 2))
    {
      return g_tstacki23[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] +
	g_tstacki23[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]];
    }
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[g_seq[i]][g_seq[j]][0][0];
	  loopEnergy += g_tstacki[g_seq[jj]][g_seq[ii]][0][0];
	}
      else
	{
	  loopEnergy += g_tstacki[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
	  loopEnergy += g_tstacki[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }
}

ENERGY HybridSSMin::QBI(int i, int j)
{
  int d, ii, jj;
  ENERGY energy = INFINITY;

  if (g_noisolate)
    for (d = j - i - 3; d >= TURN + 3 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = i + 1; ii < j - d && ii < g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(i, j, ii, jj) + Es(ii, jj) + Qprime(ii + 1, jj - 1));
	}
  else
    for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = i + 1; ii < j - d && ii <= g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(i, j, ii, jj) + Qprime(ii, jj));
	}

  return energy;
}

ENERGY HybridSSMin::QBI2(int i, int j)
{
  int d, ii, jj;
  ENERGY energy = INFINITY;

  if (g_noisolate)
    for (d = j - i - 3; d >= 2 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii < g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(jj - g_len, ii, j - g_len, i) + Es(ii, jj) + Qprime(ii + 1, jj - 1));
	}
  else
    for (d = j - i - 3; d >= 1 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii <= g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(jj - g_len, ii, j - g_len, i) + Qprime(ii, jj));
	}

  return energy;
}

ENERGY* HybridSSMin::recalloc2(ENERGY* ptr, int n)
{
  return (double*) xrealloc(ptr, n * n * sizeof(ENERGY));
}

ENERGY* HybridSSMin::recalloc2_double(ENERGY* ptr, int n)
{
  return (double*) xrealloc(ptr, n * n * sizeof(ENERGY));
}

ENERGY HybridSSMin::min4(ENERGY a, ENERGY b, ENERGY c, ENERGY d)
{
  if (a < b && a < c && a < d)
    return a;
  else if (b < c && b < d)
    return b;
  else if (c < d)
    return c;
  else
    return d;
}

ENERGY HybridSSMin::min5(ENERGY a, ENERGY b, ENERGY c, ENERGY d, ENERGY e)
{
  if (a < b && a < c && a < d && a < e)
    return a;
  else if (b < c && b < d && b < e)
    return b;
  else if (c < d && c < e)
    return c;
  else if (d < e)
    return d;
  else
    return e;
}

int HybridSSMin::equal(ENERGY a, ENERGY b)
{
#ifdef INTEGER
  return a == b;
#endif

  if (!isFinite(a) || !isFinite(b))
    return 0;

  /* 2004-06-25: replaced relative difference with line below
     so that very small numbers compare equal to 0 */
  return fabs(a - b) < 1e-5;
}

void HybridSSMin::push(struct stackNode** stack, int i, int j, int matrix)
{
  struct stackNode* new_top;

  new_top = (stackNode*) xmalloc(sizeof(struct stackNode));
  new_top->i = i;
  new_top->j = j;
  new_top->matrix = matrix;
  new_top->next = *stack;
  *stack = new_top;
}

void HybridSSMin::pushPairList(int i, int j, int length, ENERGY E)
{
  struct pairListNode* node;

  node = (pairListNode*) xmalloc(sizeof(struct pairListNode));
  node->i = i;
  node->j = j;
  node->length = length;
  node->E = E;
  node->next = pairList;
  pairList = node;
}

int HybridSSMin::isCircular()
{
  return g_prev[0] == g_len && g_next[g_len - 1] % g_len == 1;
}

double HybridSSMin::Ed3(int i, int j, int k)
{
  if (g_nodangle)
    return HUGE_VAL;

  if (k == i + 1)
    return g_dangle3[g_seq[i - 1]][g_seq[j - 1]][g_seq[k - 1]];
  else if (k == j + 1)
    return g_dangle3[g_seq[j - 1]][g_seq[i - 1]][g_seq[k - 1]];
  else
    {
      fprintf(stderr, "Error: Ed3(%d, %d, %d)\n", i, j, k);
      return 0.0;
    }
}

double HybridSSMin::chooseDangle(int a, int b)
{
  double energy;

  if (b == a + 1 || g_nodangle)
    return 0;
  else if (b == a + 2)
    {
      double d5, d3;
      d5 = Ed5(b, g_bp[b - 1], b - 1);
      d3 = Ed3(g_bp[a - 1], a, a + 1);

      if (d3 <= d5 && (d3 <= 0.0))
	{
	  return d3;
	}
      else if (d5 <= d3 && (d5 <= 0.0))
	{
	  return d5;
	}
      else
	return 0.0;
    }
  else
    {
      energy = 0;
      if (Ed3(g_bp[a - 1], a, a + 1) <= 0.0)
	{
	  energy += Ed3(g_bp[a - 1], a, a + 1);
	}
      if (Ed5(b, g_bp[b - 1], b - 1) <= 0.0)
	{
	  energy += Ed5(b, g_bp[b - 1], b - 1);
	}
      return energy;
    }
}

double HybridSSMin::tstackOrDangle(int i, int j, int external)
{
  if (j > 1 &&
	   g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	   g_next[i - 1] == i + 1 && g_prev[i] == i &&
	   g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	   g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	   g_bp[j - 2] == 0 && g_bp[i] == 0)
    {
      return external ? Etstacke(i, j) : Etstackm(i, j);
    }
  else if (j > 1 &&
	   g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	   g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	   g_bp[j - 2] == 0)
    {
      return Ed5(i, j, j - 1);
    }
  else if (g_next[i - 1] == i + 1 && g_prev[i] == i &&
	   g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	   g_bp[i] == 0)
    {
      return Ed3(i, j, i + 1);
    }
  return 0.0;
}

double HybridSSMin::Ed5(int i, int j, int k)
{
  if (g_nodangle)
    return HUGE_VAL;

  if (k == j - 1)
    return g_dangle5[g_seq[i - 1]][g_seq[j - 1]][g_seq[k - 1]];
  else if (k == i - 1)
    return g_dangle5[g_seq[j - 1]][g_seq[i - 1]][g_seq[k - 1]];
  else
    {
      fprintf(stderr, "Error: Ed5(%d, %d, %d)\n", i, j, k);
      return 0.0;
    }
}

int HybridSSMin::isHomodimer()
{
  int i;

  if (g_len % 2 == 1)
    return 0;

  if (g_next[g_len / 2 - 1] || g_prev[g_len / 2])
    return 0;

  for (i = 1; i <= g_len / 2; ++i)
    if (g_seq[i - 1] != g_seq[g_len / 2 + i - 1])
      return 0;
    else if (i != g_len / 2 && g_next[i - 1] != i + 1)
      return 0;
    else if (g_prev[i - 1] != i - 1)
      return 0;
    else if (i != g_len / 2 && g_next[g_len / 2 + i - 1] != g_len / 2 + i + 1)
      return 0;
    else if (i != 1 && g_prev[g_len / 2 + i - 1] != g_len / 2 + i - 1)
     return 0;

  return 1;
}

double HybridSSMin::Ee(int i, int j, int ds)
{
  int count;
  double energy;
  struct stack_node* new_top;

  energy = g_misc[5] + auPenalty(i, j) + (isHomodimer() ? g_misc[12] / 1.75 * log(2.0) : 0.0);

  if (ds)
    {
    if (g_hasStackingInfo)
		energy += tstackOrDangle(i, j, 1);
    else
		energy += chooseDangle(stack->j, j);
    for (new_top = stack, count = 0; count < ds - 1; new_top = new_top->next, ++count)
	{
	  if (g_hasStackingInfo)
	    energy += tstackOrDangle(new_top->j, new_top->i, 1);
	  else
	    energy += chooseDangle(new_top->next->j, new_top->i);
	  energy += auPenalty(new_top->i, new_top->j);
	}
      energy += auPenalty(new_top->i, new_top->j);
      if (g_hasStackingInfo)
	energy += tstackOrDangle(new_top->j, new_top->i, 1);
      else
	energy += chooseDangle(i, new_top->i);
    }
  else if (!g_nodangle)
    {
		if (g_hasStackingInfo)
			energy += tstackOrDangle(i, j, 1);
		else
		{
			if (g_next[i - 1] && g_prev[i] && (Ed3(i, j, i + 1) <= 0.0))
			{
				energy += Ed3(i, j, i + 1);
			}
			if (g_prev[j - 1] && g_next[j - 2] && (Ed5(i, j, j - 1) <= 0.0))
			{
				energy += Ed5(i, j, j - 1);
			}
		}
	}

  return energy;
}

void HybridSSMin::sortPairList()
{
  struct pairListNode *a, *b;

  for (a = pairList; a; a = a->next)
    for (b = a->next; b; b = b->next)
      if (a->E > b->E || (equal(a->E, b->E) && a->length < b->length))
	{
	  int iTemp, jTemp, lengthTemp;
	  ENERGY eTemp;
	  iTemp = a->i;
	  jTemp = a->j;
	  lengthTemp = a->length;
	  eTemp = a->E;
	  a->i = b->i;
	  a->j = b->j;
	  a->length = b->length;
	  a->E = b->E;
	  b->i = iTemp;
	  b->j = jTemp;
	  b->length = lengthTemp;
	  b->E = eTemp;
	}
}

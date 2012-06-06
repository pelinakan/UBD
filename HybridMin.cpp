#include "HybridMin.h"
#include "util.h"

#include <stdlib.h>
#include <stdio.h>

HybridMin::HybridMin(void)
{
	lprime = rprime = NULL;
	#if ENABLE_FORCE
	g_ssok1 = g_ssok2 = NULL;
	#endif
	//This is DNA!!!
	NA = 1;

	gotSeqs = 0;
	g_allPairs = 0;
	g_nodangle = 0;
	g_maxLoop = 30;
	tMin = 37;
	tInc = 1;
	tMax = 37;
	suffix = NULL;
	g_prefix = NULL;
	naConc = 1;
	mgConc = 0;
	polymer = 0;
	g_prefilter1 = g_prefilter2 = 2;
	noIsolate = 0;
	forceList = NULL;
	skipTraceback = 0;
	g_mfoldMax = 0;
	g_mfoldP = 0;
	g_mfoldW = 0;
	g_quiet = 0;
	g_zip = 0;
	constraints = 0;
	constraintsFile = g_bpFile = NULL;
	bestI = bestJ = 0;

	g_oneTemp = (tMin + tInc > tMax) ? 1 : 0;
  
	if (g_maxLoop < 0)
		g_maxLoop = 999999999;

	saltCorrection = ion(NA, polymer, naConc, mgConc);

	loadStack(stackEnergies, stackEnthalpies, NA, saltCorrection);
	symmetryCheckStack(stackEnergies, "energy");
	if (!g_nodangle)
		loadDangle(dangleEnergies3, dangleEnthalpies3, dangleEnergies5, dangleEnthalpies5, NA, saltCorrection);
	loadLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, NA, saltCorrection);
	loadSint2(sint2Energies, sint2Enthalpies, NA, saltCorrection);
	symmetryCheckSint2(sint2Energies, "energy");
	loadAsint1x2(asint1x2Energies, asint1x2Enthalpies, NA, saltCorrection);
	loadSint4(sint4Energies, sint4Enthalpies, NA, saltCorrection);
	symmetryCheckSint4(sint4Energies, "energy");
	loadTstacki(tstackiEnergies, tstackiEnthalpies, NA);
	loadTstacki23(tstacki23Energies, tstacki23Enthalpies, NA);
	if (!g_nodangle)
		loadTstacke(tstackeEnergies, tstackeEnthalpies, NA, saltCorrection);
	loadMisc(miscEnergies, miscEnthalpies, NA);

	//Do other initializations
	tRatio = (t + 273.15) / 310.15;
	RT = R * (t + 273.15);

	if (!suffix && (!g_oneTemp || g_firstSeq))
	{
	    combineStack(stackEnergies, stackEnthalpies, tRatio, g_stack);
	    if (!g_nodangle)
			combineDangle(dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, g_dangle3, g_dangle5);
	    combineLoop(interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnergies, interiorLoopEnthalpies, bulgeLoopEnthalpies, hairpinLoopEnthalpies, tRatio, g_interiorLoop, g_bulgeLoop, g_hairpinLoop);
	    combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
	    combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
	    combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
	    combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
	    combineTstack(tstacki23Energies, tstacki23Enthalpies, tRatio, g_tstacki23);
	    if (!g_nodangle)
			combineTstack2(tstackeEnergies, tstackeEnthalpies, tRatio, g_tstacke);
	    combineMisc(miscEnergies, miscEnthalpies, tRatio, g_misc);
	}
	makeAUPenalty(g_misc, g_aup, 0);
	makeAUPenaltyH(miscEnthalpies, g_aupH, 0);

}
 


HybridMin::~HybridMin(void)
{
}

double HybridMin::compute(double& dG, double& dH, const char* sequence1, const char* sequence2, double temperature)
{
	tMin = temperature;
	tMax = temperature;

	g_len1 = strlen(sequence1);
	g_len2 = strlen(sequence2);
	g_string1 = new char[g_len1+1];
	g_string2 = new char[g_len2+1];
	strcpy(g_string1,sequence1);
	strcpy(g_string2,sequence2);
	/* convert sequence to numbers for speed */
      g_seq1 = new unsigned char[g_len1 + 2];
      g_seq2 = new unsigned char[g_len2 + 2];
      for (i = 1; i <= g_len1; ++i)
	g_seq1[i] = util::toNum(g_string1[i - 1]);
      for (i = 1; i <= g_len2; ++i)
	g_seq2[i] = util::toNum(g_string2[i - 1]);
      g_seq1[0] = g_seq1[g_len1 + 1] = g_seq2[0] = g_seq2[g_len2 + 1] = 5;

      lprime = new double[g_len1 * g_len2];
      if (g_mfoldMax)
		rprime = new double[g_len1 * g_len2];

    for (t = tMin; t <= tMax; t += tInc)
	{

	  g_homodimer = (g_len1 != g_len2 || util::seqcmp(g_seq1, g_seq2, g_len1 + 2)) ? 0.0 : floor(RT * log(2.0) * PRECISION + 0.5);

	  initializeMatrices();
#if ENABLE_FORCE
	  force();
#endif
	  if (g_prefilter1 > 1 && !g_allPairs)
	    prefilter();
	  if (noIsolate)
	    {
	      fillMatrixL_noI(RT);
	      if (g_mfoldMax)
		fillMatrixR_noI(RT);
	    }
	  else
	    {
	      fillMatrixL(RT);
	      if (g_mfoldMax)
		fillMatrixR(RT);
	    }

	  Eleft = INFINITY;
	  if (noIsolate)
	    {
	      for (i = 1; i < g_len1; ++i)
		for (j = g_len2; j > 1; --j)
		  if (Lprime(i, j) + Es(i, j) + R0(i + 1, j - 1) < Eleft)
		    {
		      Eleft = Lprime(i, j) + Es(i, j) + R0(i + 1, j - 1);
		      bestI = i;
		      bestJ = j;
		    }
	    }
	  else
	    {
	      for (i = 1; i <= g_len1; ++i)
			for (j = g_len2; j >= 1; --j)
				if (Lprime(i, j) + R0(i, j) < Eleft)
				{
				  Eleft = Lprime(i, j) + R0(i, j);
				  bestI = i;
				  bestJ = j;
				}
	    }

	  if (!isFinite(Eleft))
	    bestI = bestJ = 1;

	  if (!skipTraceback)
	    {
	      int *bp1, *bp2, *upst1, *upst2, *dnst1, *dnst2;

	      bp1 = new int[g_len1];
	      bp2 = new int[g_len2];
	      upst1 = new int[g_len1];
	      upst2 = new int[g_len2];
	      dnst1 = new int[g_len1];
	      dnst2 = new int[g_len2];
	      for (i = 0; i < g_len1; ++i)
		bp1[i] = upst1[i] = dnst1[i] = 0;
	      for (j = 0; j < g_len2; ++j)
		bp2[j] = upst2[j] = dnst2[j] = 0;

	    if (g_mfoldMax)
		{
		  char** found;
		  int structures, *pnum1, *pnum2;
		  ENERGY cutoff;

		  pnum1 = new int[g_len1];
		  pnum2 = new int[g_len2];
		  for (i = 1; i <= g_len1; ++i)
		    pnum1[i - 1] = 0;
		  for (j = 1; j <= g_len2; ++j)
		    pnum2[j - 1] = 0;

		  cutoff = (Eleft + g_misc[5] + g_homodimer) * (1.0 - g_mfoldP / 100.0);
		  if (cutoff > Eleft + g_misc[5] + g_homodimer + 12 * PRECISION)
		    cutoff = Eleft + g_misc[5] + g_homodimer + 12 * PRECISION;
		  else if (cutoff < Eleft + g_misc[5] + g_homodimer + 1 * PRECISION)
		    cutoff = Eleft + g_misc[5] + g_homodimer + 1 * PRECISION;

		  if (noIsolate)
		    makePairList_noI(cutoff, pnum1, pnum2);
		  else
		    makePairList(cutoff, pnum1, pnum2);

		  free(pnum1);
		  free(pnum2);

		  found = new char*[g_len1];
		  for (i = 1; i <= g_len1; ++i)
		    {
		      found[i - 1] = new char[g_len2];
		      for (j = 1; j <= g_len2; ++j)
			found[i - 1][j - 1] = 0;
		    }

		  sortPairList();

		  structures = 0;
		  do {
		    if (noIsolate)
		      {
			traceback_noI(pairList->i, pairList->j, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, NULL);
			traceback_rev_noI(pairList->i + 1, pairList->j - 1, RT, bp1, bp2, upst1, upst2, dnst1, dnst2);
		      }
		    else
		      {
			traceback(pairList->i, pairList->j, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, NULL);
			traceback_rev(pairList->i, pairList->j, RT, bp1, bp2, upst1, upst2, dnst1, dnst2);
		      }
		    if (unique(bp1, found) >= g_mfoldW)
		      {
			//writeStructure(bp1, bp2, upst1, upst2, dnst1, dnst2, t, pairList->E, NULL);
			g_firstSeq = 0;
		      }
		    /* if (structures == 0)
		       writePlotExt(bp1, bp2, t, pairList->E);*/
		    for (i = 1; i <= g_len1; ++i)
		      if (bp1[i - 1])
			{
			  int ii, jj;
			  for (ii = i > g_mfoldW ? i - g_mfoldW : 1; ii <= i + g_mfoldW && ii <= g_len1; ++ii)
			    for (jj = bp1[i - 1] > g_mfoldW ? bp1[i - 1] - g_mfoldW : 1; jj <= bp1[i - 1] + g_mfoldW && jj <= g_len2; ++jj)
			      ++found[ii - 1][jj - 1];
			}
		    ++structures;
		    for (i = 0; i < g_len1; ++i)
		      bp1[i] = upst1[i] = dnst1[i] = 0;
		    for (j = 0; j < g_len2; ++j)
		      bp2[j] = upst2[j] = dnst2[j] = 0;

		    while (pairList && found[pairList->i - 1][pairList->j - 1])
		      {
			struct pairListNode* newTop;
			newTop = pairList->next;
			free(pairList);
			pairList = newTop;
		      }
		  } while (pairList && structures < g_mfoldMax);

		  for (i = 1; i <= g_len1; ++i)
		    free(found[i - 1]);
		  free(found);
		}
	      else
		{
		  double enthalpy;
		  enthalpy = 0.0;
		  if (isFinite(Lprime(bestI, bestJ)))
		    {
		      if (noIsolate)
			traceback_noI(bestI, bestJ, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, suffix ? NULL : &enthalpy);
		      else
			traceback(bestI, bestJ, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, suffix ? NULL : &enthalpy);
		      if (!suffix)
			enthalpy += auPenaltyH(g_seq1[bestI], g_seq2[bestJ]) + miscEnthalpies[5];
		      if (!g_nodangle)
			{
			  if (bestI < g_len1 && bestJ > 1 &&
			      (g_zip || (g_tstacke[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]] <= g_dangle3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]] &&
					 g_tstacke[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]] <= g_dangle5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]] &&
					 g_tstacke[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]] <= 0.0)))
			    {
			      upst1[bestI] = bestI;
			      dnst1[bestI - 1] = bestI + 1;
			      upst2[bestJ - 1] = bestJ - 1;
			      dnst2[bestJ - 2] = bestJ;
			      if (!suffix)
				enthalpy += tstackeEnthalpies[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]];
			    }
			  else if (bestI < g_len1 &&
				   g_dangle3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]] <= g_dangle5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]] &&
				   g_dangle3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]] <= 0.0)
			    {
			      upst1[bestI] = bestI;
			      dnst1[bestI - 1] = bestI + 1;
			      if (!suffix)
				enthalpy += dangleEnthalpies3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]];
			    }
			  else if (bestJ > 1 && g_dangle5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]] <= 0.0)
			    {
			      upst2[bestJ - 1] = bestJ - 1;
			      dnst2[bestJ - 2] = bestJ;
			      if (!suffix)
				enthalpy += dangleEnthalpies5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]];
			    }
			}
		    }
			
			dG = (double) (g_misc[5] + Eleft + g_homodimer) / PRECISION;
			dH = enthalpy;

		   
		}

	      delete[] bp1;
	      delete[] bp2;
	      delete[] upst1;
	      delete[] upst2;
	      delete[] dnst1;
	      delete[] dnst2;
	    }
	  else if (g_quiet)
	    {
	      puts("");
	      fflush(stdout);
	    }
	}

	delete[] g_string1;
	delete[] g_string2;
	delete[] g_seq1;
	delete[] g_seq2;
	delete[] lprime;
	if (rprime)
		delete[] rprime;

    return 0;
}

void HybridMin::initializeMatrices()
{
  /* initialize L', R' to +infinity for illegal pairs, 0 otherwise */
  int i, j;

  for (i = 1; i <= g_len1; ++i)
    for (j = g_len2; j >= 1; --j)
      if (basePairIndex(g_seq1[i], g_seq2[j]) == 6 && !g_allPairs)
	{
	  Lprime(i, j) = INFINITY;
	  if (rprime)
	    Rprime(i, j) = INFINITY;
	}
      else
	{
	  Lprime(i, j) = 0.0;
	  if (rprime)
	    Rprime(i, j) = 0.0;
	}
}

#if ENABLE_FORCE
void force()
{
  int i, j, k;
  struct constraintListNode* top;

  top = forceList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len1)
	for (i = 0; i <= g_len1 + 1; ++i)
	  for (j = i; j <= g_len1 + 1; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (i <= top->i + k && top->i + k <= j)
		ssOK1(i, j) = 0;

      if (top->j >= 1 && top->j <= g_len2)
	{
	  if (top->i == 0)
	    {
	      for (i = 0; i <= g_len2 + 1; ++i)
		for (j = i; j <= g_len2 + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j + k && top->j + k <= j)
		      ssOK2(i, j) = 0;
	    }
	  else
	    {
	      for (i = 0; i <= g_len2 + 1; ++i)
		for (j = i; j <= g_len2 + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j - k && top->j - k <= j)
		      ssOK2(i, j) = 0;
	    }
	}

      if (top->i >= 1 && top->i <= g_len1 && top->j >= 1 && top->j <= g_len2)
	{
	  for (i = 1; i <= g_len1; ++i)
	    for (k = 0; k < top->k; ++k)
	      if (i != top->i + k)
		{
		  Lprime(i, top->j - k) = INFINITY;
		  if (rprime)
		    Rprime(i, top->j - k) = INFINITY;
		}
	  for (j = 1; j <= g_len2; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (j != top->j - k)
		{
		  Lprime(top->i + k, j) = INFINITY;
		  if (rprime)
		    Rprime(top->i + k, j) = INFINITY;
		}
	}

      top = top->next;
    }
}
#endif

void HybridMin::prefilter()
{
  char** in;
  int i, j, k, count;

  in = new char*[g_len1];
  for (i = 1; i <= g_len1; ++i)
    in[i - 1] = (char*)xcalloc(g_len2, 1);

  for (i = 1; i <= g_len1 - g_prefilter2 + 1; ++i)
    for (j = g_len2; j >= g_prefilter2; --j)
      {
	count = 0;
	for (k = 0; k < g_prefilter2; ++k)
	  if (isFinite(Lprime(i + k, j - k)))
	    ++count;
	if (count >= g_prefilter1)
	  for (k = 0; k < g_prefilter2; ++k)
	    ++in[i + k - 1][j - k - 1];
      }

  for (i = 1; i <= g_len1; ++i)
    {
      for (j = g_len2; j >= 1; --j)
	if (!in[i - 1][j - 1])
	  {
	    Lprime(i, j) = INFINITY;
	    if (rprime)
	      Rprime(i, j) = INFINITY;
	  }
      free(in[i - 1]);
    }
  delete[] in;
}

void HybridMin::fillMatrixL(double RT)
{
  int d, i, j, ii, jj;

  for (i = 1; i <= g_len1; ++i)
    for (j = g_len2; j >= 1; --j)
      if (isFinite(Lprime(i, j)))
	{
	  Lprime(i, j) = L0(i, j);
	  if (i > 1 && j < g_len2)
	    Lprime(i, j) = min2(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      ii = i - 1;
	      jj = ii + d + (j - i);
	      if (jj > g_len2)
		{
		  ii -= (jj - g_len2);
		  jj = g_len2;
		}
	      for (; ii > 0 && jj > j; --ii, --jj)
		if (isFinite(Lprime(ii, jj)))
		  Lprime(i, j) = min2(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Lprime(ii, jj));
	    }
	}
}

void HybridMin::fillMatrixR(double RT)
{
  int d, i, j, ii, jj;
  
  for (j = 1; j <= g_len2; ++j)
    for (i = g_len1; i >= 1; --i)
      if (isFinite(Rprime(i, j)))
	{
	  Rprime(i, j) = R0(i, j);
	  if (j > 1 && i < g_len1)
	    Rprime(i, j) = min2(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      jj = j - 1;
	      ii = jj + d + (i - j);
	      if (ii > g_len1)
		{
		  jj -= (ii - g_len1);
		  ii = g_len1;
		}
	      for (; jj > 0 && ii > i; --ii, --jj)
		if (isFinite(Rprime(ii, jj)))
		  Rprime(i, j) = min2(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Rprime(ii, jj));
	    }
	}
}

void HybridMin::fillMatrixL_noI(double RT)
{
  int d, i, j, ii, jj;

  for (i = 1; i <= g_len1; ++i)
    for (j = g_len2; j >= 1; --j)
      if (isFinite(Lprime(i, j)))
	{
	  Lprime(i, j) = L0(i, j);
	  if (i > 1 && j < g_len2)
	    Lprime(i, j) = min2(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      ii = i - 1;
	      jj = ii + d + (j - i);
	      if (jj >= g_len2)
		{
		  ii -= (jj - g_len2 + 1);
		  jj = g_len2 - 1;
		}
	      for (; ii > 1 && jj > j; --ii, --jj)
		if (isFinite(Lprime(ii, jj)))
		  Lprime(i, j) = min2(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Es(ii - 1, jj + 1) + Lprime(ii - 1, jj + 1));
	    }
	}
}

void HybridMin::fillMatrixR_noI(double RT)
{
  int d, i, j, ii, jj;

  for (j = 1; j <= g_len2; ++j)
    for (i = g_len1; i >= 1; --i)
      if (isFinite(Rprime(i, j)))
	{
	  Rprime(i, j) = R0(i, j);
	  if (j > 1 && i < g_len1)
	    Rprime(i, j) = min2(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      jj = j - 1;
	      ii = jj + d + (i - j);
	      if (ii >= g_len1)
		{
		  jj -= (ii - g_len1 + 1);
		  ii = g_len1 - 1;
		}
	      for (; jj > 1 && ii > i; --ii, --jj)
		if (isFinite(Rprime(ii, jj)))
		  Rprime(i, j) = min2(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Es(ii, jj) + Rprime(ii + 1, jj - 1));
	    }
	}
}

void HybridMin::traceback(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double* enthalpy)
{
  int d, ii, jj, done;

  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Lprime(i, j), L0(i, j)))
	{
	  if (enthalpy)
	    *enthalpy += auPenaltyH(g_seq1[i], g_seq2[j]);
	  if (!g_nodangle)
	    {
	      if (i > 1 && j < g_len2 && L0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += tstackeEnthalpies[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]];
		}
	      else if (i > 1 && L0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]];
		}
	      else if (j < g_len2 && L0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]])
		{
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]];
		}
	    }
	  break;
	}
      done = 0;
      if (i > 1 && j < g_len2 && equal(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1)))
	{
	  i = i - 1;
	  j = j + 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i] = i;
	  upst2[j - 1] = j - 1;
	  dnst1[i - 1] = i + 1;
	  dnst2[j - 2] = j;
	  if (enthalpy)
	    *enthalpy += Hs(i, j);
	  done = 1;
	}
      for (d = 3; !done && d <= g_maxLoop + 2; ++d)
	{
	  ii = i - 1;
	  jj = ii + d + (j - i);
	  if (jj > g_len2)
	    {
	      ii -= (jj - g_len2);
	      jj = g_len2;
	    }
	  for (; !done && ii > 0 && jj > j; --ii, --jj)
	    if (equal(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Lprime(ii, jj)))
	      {
		setStackBI(ii, jj, i, j, upst1, upst2, dnst1, dnst2);
		if (enthalpy)
		  *enthalpy += Hbi(ii, jj, i, j);
		i = ii;
		j = jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		done = 1;
		break;
	      }
	}
    }
}

void HybridMin::traceback_noI(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double* enthalpy)
{
  int d, ii, jj, done;

  bp1[i] = j - 1;
  bp2[j - 2] = i + 1;
  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Lprime(i, j), L0(i, j)))
	{
	  if (enthalpy)
	    *enthalpy += auPenaltyH(g_seq1[i], g_seq2[j]);
	  if (!g_nodangle)
	    {
	      if (i > 1 && j < g_len2 && L0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += tstackeEnthalpies[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]];
		}
	      else if (i > 1 && L0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]];
		}
	      else if (j < g_len2 && L0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]])
		{
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]];
		}
	    }
	  break;
	}
      done = 0;
      if (i > 1 && j < g_len2 && equal(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1)))
	{
	  i = i - 1;
	  j = j + 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i] = i;
	  upst2[j - 1] = j - 1;
	  dnst1[i - 1] = i + 1;
	  dnst2[j - 2] = j;
	  if (enthalpy)
	    *enthalpy += Hs(i, j);
	  done = 1;
	}
      for (d = 3; !done && d <= g_maxLoop + 2; ++d)
	{
	  ii = i - 1;
	  jj = ii + d + (j - i);
	  if (jj >= g_len2)
	    {
	      ii -= (jj - g_len2 + 1);
	      jj = g_len2 - 1;
	    }
	  for (; !done && ii > 1 && jj > j; --ii, --jj)
	    if (equal(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Es(ii - 1, jj + 1) + Lprime(ii - 1, jj + 1)))
	      {
		setStackBI(ii, jj, i, j, upst1, upst2, dnst1, dnst2);
		if (enthalpy)
		  *enthalpy += Hbi(ii, jj, i, j) + Hs(ii - 1, jj + 1);
		i = ii;
		j = jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		upst1[i - 1] = i - 1;
		dnst2[j - 1] = j + 1;
		--ii; ++jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		upst2[j - 1] = j - 1;
		dnst1[i - 1] = i + 1;
		done = 1;
		break;
	      }
	}
    }
}

void HybridMin::traceback_rev(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2)
{
  int ii, jj, done;

  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Rprime(i, j), R0(i, j)))
	{
	  if (!g_nodangle)
	    {
	      if (i < g_len1 && j > 1 && R0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	      else if (i < g_len1 && R0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		}
	      else if (j > 1 && R0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]])
		{
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	    }
	  break;
	}
      done = 0;
      if (j > 1 && i < g_len1 && equal(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1)))
	{
	  i = i + 1;
	  j = j - 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i - 1] = i - 1;
	  upst2[j] = j;
	  dnst1[i - 2] = i;
	  dnst2[j - 1] = j + 1;
	  done = 1;
	}

      for (ii = i + 1; !done && ii <= g_len1; ++ii)
	for (jj = j - 1; !done && jj >= 1; --jj)
	  if (ii - i + j - jj > 2 && ii - i + j - jj <= 2 + g_maxLoop)
	    if (equal(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Rprime(ii, jj)))
	      {
		setStackBI(i, j, ii, jj, upst1, upst2, dnst1, dnst2);
		i = ii;
		j = jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		done = 1;
		break;
	      }
    }
}

void HybridMin::traceback_rev_noI(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2)
{
  int ii, jj, done;

  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Rprime(i, j), R0(i, j)))
	{
	  if (!g_nodangle)
	    {
	      if (i < g_len1 && j > 1 && R0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	      else if (i < g_len1 && R0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		}
	      else if (j > 1 && R0(i, j) == hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]])
		{
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	    }
	  break;
	}
      done = 0;
      if (j > 1 && i < g_len1 && equal(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1)))
	{
	  i = i + 1;
	  j = j - 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i - 1] = i - 1;
	  upst2[j] = j;
	  dnst1[i - 2] = i;
	  dnst2[j - 1] = j + 1;	  
	  done = 1;
	}

      for (ii = i + 1; !done && ii < g_len1; ++ii)
	for (jj = j - 1; !done && jj > 1; --jj)
	  if (ii - i + j - jj > 2 && ii - i + j - jj <= 2 + g_maxLoop)
	    if (equal(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Es(ii, jj) + Rprime(ii + 1, jj - 1)))
	      {
		setStackBI(i, j, ii, jj, upst1, upst2, dnst1, dnst2);
		i = ii + 1;
		j = jj - 1;
		bp1[i - 2] = j + 1;
		bp2[j] = i - 1;
		upst2[j] = j;
		dnst1[i - 2] = i;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		upst1[i - 1] = i - 1;
		dnst2[j - 1] = j + 1;
		done = 1;
		break;
	      }
    }
}

void HybridMin::setStackBI(int i, int j, int ii, int jj, int* upst1, int* upst2, int* dnst1, int* dnst2)
{
  int loopSize1, loopSize2;

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;

#ifdef DEBUG
  if (loopSize1 < 0 || loopSize2 < 0 || (loopSize1 == 0 && loopSize2 == 0))
    {
      fputs("Error: setStackBI() called with nonsense\n", stderr);
      return;
    }
  else
#endif

  if ((loopSize1 == 0 && loopSize2 == 1) || (loopSize2 == 0 && loopSize1 == 1))
    {
      upst1[ii - 1] = i;
      dnst1[i - 1] = ii;
      upst2[j - 1] = jj;
      dnst2[jj - 1] = j;
    }
  else if (loopSize1 && loopSize2 && (loopSize1 > 2 || loopSize2 > 2))
    {
      upst1[i] = i;
      upst1[ii - 1] = ii - 1;
      dnst1[i - 1] = i + 1;
      dnst1[ii - 2] = ii;
      upst2[jj] = jj;
      upst2[j - 1] = j - 1;
      dnst2[jj - 1] = jj + 1;
      dnst2[j - 2] = j;
    }
}

int HybridMin::unique(int* bp1, char** found)
{
  int i, unique_;

  unique_ = 0;
  for (i = 1; i <= g_len1; ++i)
    if (bp1[i - 1] && !found[i - 1][bp1[i - 1] - 1])
      ++unique_;

  return unique_;
}

void HybridMin::makePairList(ENERGY cutoff, int* pnum1, int* pnum2)
{
  int d, i, j, length;
  ENERGY E;

  pairList = NULL;
  for (d = 1; d < g_len1 + g_len2; ++d)
    {
      if (d > g_len1)
	{
	  i = g_len1;
	  j = d + 1 - g_len1;
	}
      else
	{
	  i = d;
	  j = 1;
	}
      length = 0;
      E = INFINITY;
      for (; i >= 1 && j <= g_len2; --i, ++j)
	{
	  if (Lprime(i, j) + Rprime(i, j) + g_misc[5] + g_homodimer <= cutoff)
	    {
	      ++pnum1[i - 1];
	      ++pnum2[j - 1];
	    }

	  if (length && equal(Lprime(i, j) + Rprime(i, j), E))
	    ++length;
	  else if (isFinite(Lprime(i, j)))
	    {
	      if (length && E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 1;
	      E = Lprime(i, j) + Rprime(i, j);
	    }
	  else if (length)
	    {
	      if (E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 0;
	    }
	  if ((i == 1 || j == g_len2) && length)
	    if (E + g_misc[5] + g_homodimer <= cutoff)
	      pushPairList(i, j, length, E);
	}
    }
}

void HybridMin::makePairList_noI(ENERGY cutoff, int* pnum1, int* pnum2)
{
  int d, i, j, length;
  ENERGY E;

  pairList = NULL;
  for (d = 2; d < g_len1 + g_len2 - 1; ++d)
    {
      if (d > g_len1)
	{
	  i = g_len1 - 1;
	  j = d + 2 - g_len1;
	}
      else
	{
	  i = d - 1;
	  j = 2;
	}
      length = 0;
      E = INFINITY;
      for (; i >= 1 && j <= g_len2; --i, ++j)
	{
	  if (Lprime(i, j) + Es(i, j) + Rprime(i + 1, j - 1) + g_misc[5] + g_homodimer <= cutoff)
	    {
	      ++pnum1[i - 1];
	      ++pnum2[j - 1];
	      ++pnum1[i];
	      ++pnum2[j - 2];
	    }

	  if (length && equal(Lprime(i, j) + Es(i, j) + Rprime(i + 1, j - 1), E))
	    ++length;
	  else if (isFinite(Lprime(i, j) && isFinite(Es(i, j)) && isFinite(Rprime(i + 1, j - 1))))
	    {
	      if (length && E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 1;
	      E = Lprime(i, j) + Es(i, j) + Rprime(i + 1, j - 1);
	    }
	  else if (length)
	    {
	      if (E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 0;
	    }
	  if ((i == 1 || j == g_len2) && length)
	    if (E + g_misc[5] + g_homodimer <= cutoff)
	      pushPairList(i, j, length + 1, E);
	}
    }
}

ENERGY HybridMin::Es(int i, int j)
{
  return g_stack[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
}

double HybridMin::Hs(int i, int j)
{
  return stackEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
}

ENERGY HybridMin::Ebi(int i, int j, int ii, int jj, double RT)
{
  int loopSize1, loopSize2;
  ENERGY loopEnergy, asPenalty;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
#endif

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
#ifdef DEBUG
  if (loopSize1 + loopSize2 > g_maxLoop)
    return INFINITY;
#endif

#if ENABLE_FORCE
  if (loopSize1 && !ssOK1(i + 1, ii - 1))
    return INFINITY;
  if (loopSize2 && !ssOK2(jj + 1, j - 1))
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
	return g_bulgeLoop[0] + g_stack[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] + hybAuPenalty(g_seq1[i], g_seq2[j]) + hybAuPenalty(g_seq1[ii], g_seq2[jj]);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + hybAuPenalty(g_seq1[i], g_seq2[j]) + hybAuPenalty(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] + g_stack[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] + hybAuPenalty(g_seq1[i], g_seq2[j]) + hybAuPenalty(g_seq1[ii], g_seq2[jj]);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + hybAuPenalty(g_seq1[i], g_seq2[j]) + hybAuPenalty(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq2[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(g_seq2[jj], g_seq1[ii])][basePairIndex(g_seq2[j], g_seq1[i])][g_seq2[jj + 1]][g_seq1[ii - 1]][g_seq1[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq1[i + 2]][g_seq2[j - 2]];
  else if ((loopSize1 == 2 && loopSize2 == 3) ||
	   (loopSize1 == 3 && loopSize2 == 2))
    {
      return g_tstacki23[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]] +
	g_tstacki23[g_seq2[jj]][g_seq1[ii]][g_seq2[jj + 1]][g_seq1[ii - 1]];
    }
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[g_seq1[i]][g_seq2[j]][0][0];
	  loopEnergy += g_tstacki[g_seq2[jj]][g_seq1[ii]][0][0];
	}
      else
	{
	  loopEnergy += g_tstacki[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
	  loopEnergy += g_tstacki[g_seq2[jj]][g_seq1[ii]][g_seq2[jj + 1]][g_seq1[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }

}

double HybridMin::Hbi(int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Hbi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Hbi(): jj isn't less than j\n", stderr);
#endif

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > g_maxLoop)
    return INFINITY;

#if ENABLE_FORCE
  if (loopSize1 && !ssOK1(i + 1, ii - 1))
    return INFINITY;
  if (loopSize2 && !ssOK2(jj + 1, j - 1))
    return INFINITY;
#endif

#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Hbi() called with nonsense\n", stderr);
      return INFINITY;
    }
  else
#endif
  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return bulgeLoopEnthalpies[0] + stackEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize2 <= 30)
	return bulgeLoopEnthalpies[loopSize2 - 1] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
      else
	return bulgeLoopEnthalpies[29] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return bulgeLoopEnthalpies[0] + stackEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize1 <= 30)
	return bulgeLoopEnthalpies[loopSize1 - 1] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
      else
	return bulgeLoopEnthalpies[29] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return sint2Enthalpies[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return asint1x2Enthalpies[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq2[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return asint1x2Enthalpies[basePairIndex(g_seq2[jj], g_seq1[ii])][basePairIndex(g_seq2[j], g_seq1[i])][g_seq2[jj + 1]][g_seq1[ii - 1]][g_seq1[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return sint4Enthalpies[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq1[i + 2]][g_seq2[j - 2]];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = interiorLoopEnthalpies[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = interiorLoopEnthalpies[29];
      if (miscEnthalpies[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += tstackiEnthalpies[g_seq1[i]][g_seq2[j]][0][0];
	  loopEnergy += tstackiEnthalpies[g_seq2[jj]][g_seq1[ii]][0][0];
	}
      else
	{
	  loopEnergy += tstackiEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
	  loopEnergy += tstackiEnthalpies[g_seq2[jj]][g_seq1[ii]][g_seq2[jj + 1]][g_seq1[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * miscEnthalpies[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > miscEnthalpies[4])
	asPenalty = miscEnthalpies[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }

}

ENERGY HybridMin::R0(int i, int j)
{
  ENERGY energy;

  if (basePairIndex(g_seq1[i], g_seq2[j]) == 6)
    return INFINITY;

#if ENABLE_FORCE
  if (!ssOK1(i + 1, g_len1) || !ssOK2(1, j - 1))
    return INFINITY;
#endif

  if (g_nodangle)
    return hybAuPenalty(g_seq1[i], g_seq2[j]);

  energy = hybAuPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
  if (g_zip)
    return energy;

  energy = min2(energy, hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]]);
  energy = min2(energy, hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]]);
  energy = min2(energy, hybAuPenalty(g_seq1[i], g_seq2[j]));
  return energy;
}

ENERGY HybridMin::L0(int i, int j)
{
  ENERGY energy;

  if (basePairIndex(g_seq1[i], g_seq2[j]) == 6)
    return INFINITY;

#if ENABLE_FORCE
  if (!ssOK1(1, i - 1) || !ssOK2(j + 1, g_len2))
    return INFINITY;
#endif

  if (g_nodangle)
    return hybAuPenalty(g_seq1[i], g_seq2[j]);

  energy = hybAuPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]];
  if (g_zip)
    return energy;

  energy = min2(energy, hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]]);
  energy = min2(energy, hybAuPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]]);
  energy = min2(energy, hybAuPenalty(g_seq1[i], g_seq2[j]));
  return energy;
}

ENERGY* HybridMin::recalloc2(ENERGY* ptr, int m, int n)
{
  return (ENERGY*)xrealloc(ptr, m * n * sizeof(ENERGY));
}

ENERGY HybridMin::min4(ENERGY a, ENERGY b, ENERGY c, ENERGY d)
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

int HybridMin::equal(ENERGY a, ENERGY b)
{
#ifdef INTEGER
  return a == b;
#endif

  if (!finite(a) || !finite(b))
    return 0;

  /* 2004-06-25: replaced relative difference with line below
     so that very small numbers compare equal to 0 */
  return fabs(a - b) < 1e-5;
}

void HybridMin::push(struct stackNode** stack, int i, int j)
{
  struct stackNode* new_top;

  new_top = new stackNode();
  new_top->i = i;
  new_top->j = j;
  new_top->next = *stack;
  *stack = new_top;
}

void HybridMin::pushPairList(int i, int j, int length, ENERGY E)
{
  struct pairListNode* node;

  node = new pairListNode();
  node->i = i;
  node->j = j;
  node->length = length;
  node->E = E;
  node->next = pairList;
  pairList = node;
}

void HybridMin::sortPairList()
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

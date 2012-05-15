#include "CtEnergy.h"
#include <math.h>

#include "util.h"

#include <stdlib.h>
#include "hybrid-ss-min.h"
#include "energy.h"

CtEnergy::CtEnergy(void)
{
	//Set relevant variables 
	g_nodangle = 0;
	//Load all the crap here!
	char* suffix = "DHD";
	if (!g_nodangle)
	loadDangleSuffix(g_dangle3, g_dangle5, suffix);
    loadStackSuffix(g_stack, suffix);
    loadLoopSuffix(g_hairpinLoop, g_interiorLoop, g_bulgeLoop, suffix);
    loadSint2Suffix(g_sint2, suffix);
    loadAsint1x2Suffix(g_asint1x2, suffix);
    loadSint4Suffix(g_sint4, suffix);
    loadTstackhSuffix(g_tstackh, suffix);
    loadTstackiSuffix(g_tstacki, suffix);
    loadTstacki23Suffix(g_tstacki23, suffix);
    if (!g_nodangle)
	{
		loadTstackmSuffix(g_tstackm, suffix);
		loadTstackeSuffix(g_tstacke, suffix);
	}
    loadTriloopSuffix(&g_triloop, &numTriloops, suffix);
    loadTloopSuffix(&g_tloop, &numTloops, suffix);
    loadHexaloopSuffix(&g_hexaloop, &numHexaloops, suffix);
	loadMultiSuffix(g_multi, suffix);
    loadMiscSuffix(g_misc, suffix);
}

double CtEnergy::compute(HybridSSMin* hComp, int* bp, int* upst, int* dnst)
{
	g_len = hComp->g_len;
	g_bp = new int[g_len];
	memcpy(g_bp,bp,g_len*sizeof(int));
    g_upst = new int[g_len];
	memcpy(g_upst,upst,g_len*sizeof(int));
    g_dnst = new int[g_len];
	memcpy(g_dnst,dnst,g_len*sizeof(int));
    g_numbers = g_prev = g_next = NULL;

  //Do allocation of stuff!
  g_numbers = new int[g_len];
  g_prev = new int[g_len];
  g_next = new int[g_len];
  g_seq = new unsigned char[g_len];
  for (int i=0;i<g_len-1;++i) {
	  g_numbers[i]=i+1;
	  g_prev[i]=i;
	  g_next[i]=i+2;
	  g_seq[i]=hComp->g_seq[i+1];
  }
  g_seq[g_len-1]=util::toNum(hComp->g_string[g_len-1]);
  g_numbers[g_len-1]=g_len;
  g_prev[g_len-1]=g_len-1;
  g_next[g_len-1]=0;
  g_dnst[g_len-1]=0;

  //Now, reset some of the parameters
  double RT,t,tRatio;
  loadRTSuffix(&RT, "DHD");
  t = RT / R - 273.15;
  tRatio = (t + 273.15) / 310.15;
  //Recombine, as the tRatio has changed!!!!!!!!!!!!!!!!!!!


  /*for (i = 1; i < g_len; ++i)
	fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", i, g_string[i - 1], i - 1, i + 1, bp[i - 1], i, upst[i - 1], dnst[i - 1]);
      fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len, g_string[g_len - 1], g_len - 1, 0, bp[g_len - 1], g_len, upst[g_len - 1], 0);

for (i = 0; i < g_len; ++i)
    {
      if (!fgets(line, 80, file))
	return 0;
      g_upst[i] = g_dnst[i] = -1;
      count = sscanf(line, "%d %c %d %d %d %d %d %d", &num, &g_bases[i], &g_prev[i], &g_next[i], &g_bp[i], &g_numbers[i], &g_upst[i], &g_dnst[i]);
}*/
  int j,k,i;
  //We know this to be the case!
  g_hasStackingInfo = 1;
  double etotal = 0;
  stack = (stack_node*)xmalloc(sizeof(struct stack_node));

	  /* if molecule is circular, find a hairpin and start there */
	  if (ct_isCircular())
	    {
	      for (i = 1; i <= g_len; ++i)
		{
		  int flag = 0;
		  for (j = i + 1; j < g_bp[i - 1]; ++j)
		    if (g_bp[j - 1])
		      {
			++flag;
			break;
		      }
		  if (g_bp[i - 1] < i || flag)
		    continue;

		  /* circular permutation */
		  g_seq = (unsigned char*)realloc(g_seq, g_len + i);
		  g_bp = (int*)realloc(g_bp, (g_len + i) * sizeof(int));
		  g_numbers = (int*)realloc(g_numbers, (g_len + i) * sizeof(int));
		  g_prev = (int*)realloc(g_prev, (g_len + i) * sizeof(int));
		  g_next = (int*)realloc(g_next, (g_len + i) * sizeof(int));
		  g_upst = (int*)realloc(g_upst, (g_len + i) * sizeof(int));
		  g_dnst = (int*)realloc(g_dnst, (g_len + i) * sizeof(int));

		  for (j = 1; j <= i; ++j)
		    {
		      g_seq[g_len + j - 1] = g_seq[j - 1];
		      g_bp[g_len + j - 1] = 0;
		      g_numbers[g_len + j - 1] = g_numbers[j - 1];
		      g_prev[g_len + j - 1] = g_prev[j - 1];
		      g_next[g_len + j - 1] = g_next[j - 1];
		      g_upst[g_len + j - 1] = g_upst[j - 1];
		      g_dnst[g_len + j - 1] = g_dnst[j - 1];
		    }

		  for (j = 1; j <= i; ++j)
		    if (g_bp[j - 1] > 0)
		      {
			if (g_bp[j - 1] < i)
			  {
			    g_bp[g_len + j - 1] = g_len + g_bp[j - 1];
			    g_bp[g_len + g_bp[j - 1] - 1] = g_len + j;
			  }
			else
			  {
			    g_bp[g_bp[j - 1] - 1] = g_len + j;
			    g_bp[g_len + j - 1] = g_bp[j - 1];
			  }
		      }

		  for (j = i + 1; j <= g_len + i; ++j)
		    {
		      if (g_prev[j - 1] <= i)
			g_prev[j - 1] += g_len;
		      if (g_next[j - 1] <= i)
			g_next[j - 1] += g_len;
		      if (g_upst[j - 1] && g_upst[j - 1] <= i)
			g_upst[j - 1] += g_len;
		      if (g_dnst[j - 1] && g_dnst[j - 1] <= i)
			g_dnst[j - 1] += g_len;
		    }

		  etotal = ct_Eh(i, g_bp[i - 1]);
		  stack->i = g_bp[i - 1];
		  stack->j = g_len + i;
		  stack->open = 0;
		  stack->next = NULL;
		  break;
		}
	    }
	  else
	    {
	      stack->i = 1;
	      stack->j = g_len;
	      stack->open = 1;
	      stack->next = NULL;
	    }
	  int count;
	  /* start with (1,n), parse loop, push loops found onto stack */
	  while (stack)
	    {
	      i = stack->i;
	      j = stack->j;
	      open = stack->open;
	      new_top = stack->next;
	      free(stack);
	      stack = new_top;
	      ds = ss1 = ss2 = 0;
	      eloop = 0;
	      is_exterior = 0;
	      
	      for (k = open ? i : i + 1; k < (open ? j + 1 : j); ++k)
		if (g_bp[k - 1] > j)
		  {
		    fputs("Error: pseudoknot detected\n", stderr);
		    return EXIT_FAILURE;
		  }
		else if (g_bp[k - 1] > k)
		  {
		    ++ds;
		    new_top = (stack_node*)malloc(sizeof(struct stack_node));
		    new_top->i = k;
		    new_top->j = g_bp[k - 1];
		    new_top->open = 0;
		    new_top->next = stack;
		    stack = new_top;
		    k = g_bp[k - 1];
		  }
		else
		  {
		    if (!ct_isCircular() && (k == 1 || (g_next[k - 2] == 0 && g_prev[k - 1] == 0)))
		      ++is_exterior;
		    if (!ct_isCircular() && (k == g_len || (g_next[k - 1] == 0 && g_prev[k] == 0)))
		      ++is_exterior;
		    if (ds == 0)
		      ++ss1;
		    else
		      ++ss2;
		  }
	      
	      if (open)
		{
		  if (stack)
		    {
		    if (g_hasStackingInfo)
				for (new_top = stack, count = 0; count < ds; new_top = new_top->next, ++count)
				{
			    eloop += ct_tstackOrDangle(new_top->j, new_top->i, 1);
			    eloop += ct_auPenalty(new_top->i, new_top->j);
				}
		    else
			{
			  if (stack->j < g_len && !g_nodangle && (ct_Ed3(stack->i, stack->j, stack->j + 1) <= 0.0))
			    {
			      eloop += ct_Ed3(stack->i, stack->j, stack->j + 1);
			      
			    }
			  for (new_top = stack, count = 0; count < ds - 1; new_top = new_top->next, ++count)
			    {
			      eloop += ct_chooseDangle(new_top->next->j, new_top->i);
			      eloop += ct_auPenalty(new_top->i, new_top->j);
				}
			  eloop += ct_auPenalty(new_top->i, new_top->j);
			  if (new_top->i > 1 && !g_nodangle && (ct_Ed5(new_top->i, new_top->j, new_top->i - 1) <= 0.0))
			    {
			      eloop += ct_Ed5(new_top->i, new_top->j, new_top->i - 1);
			    }
			}
		    }
		}
	      else if (is_exterior)
		eloop = ct_Ee(i, j, ds);
	      else
		{
		  if (ds == 0)
		    eloop = ct_Eh(i, j);
		  else if (ds == 1)
		    {
		      if (ss1 == 0 && ss2 == 0)
			eloop = ct_Es(i, j);
		      else
			eloop = ct_Ebi(i, j, i + ss1 + 1, j - ss2 - 1);
		    }
		  else
		    {
			  eloop = g_multi[0] + g_multi[1] * (ss1 + ss2) + g_multi[2] * (ds + 1);

		      if (g_hasStackingInfo)
			eloop += ct_tstackOrDangle(i, j, 0);
		      else
			eloop += ct_chooseDangle(stack->j, j);
		      for (new_top = stack, count = 0; count < ds - 1; new_top = new_top->next, ++count)
			{
			  if (g_hasStackingInfo)
			    eloop += ct_tstackOrDangle(new_top->j, new_top->i, 0);
			  else
			    eloop += ct_chooseDangle(new_top->next->j, new_top->i);
			  eloop += ct_auPenalty(new_top->i, new_top->j);
			  
			}
		      eloop += ct_auPenalty(new_top->i, new_top->j);
		      
		      if (g_hasStackingInfo)
			eloop += ct_tstackOrDangle(new_top->j, new_top->i, 0);
		      else
			eloop += ct_chooseDangle(i, new_top->i);
		      eloop += ct_auPenalty(i, j);
		    }
		}
	      etotal += eloop;
	    }
	
	//Cleanup
	delete[] g_bp;
	delete[] g_upst;
	delete[] g_dnst;
	delete[] g_numbers;
	delete[] g_prev;
	delete[] g_next;
	delete[] g_seq;

	return etotal;
	
}

CtEnergy::~CtEnergy(void)
{
}

double CtEnergy::ct_Etstackm(int i, int j)
{
  if (g_nodangle)
    return HUGE_VAL;

  return g_tstackm[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
}

double CtEnergy::ct_Etstacke(int i, int j)
{
  if (g_nodangle)
    return HUGE_VAL;

  return g_tstacke[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
}

double CtEnergy::ct_Eh(int i, int j)
{
  double energy = 0.0;
  int loopSize = j - i - 1;
  int k;

  /* check for external loop in dimer */
  for (k = i; k < j; ++k)
    if (g_next[k - 1] == 0 && g_prev[k] == 0)
      return ct_Ee(i, j, 0);

  if (loopSize < 3)
    {
      return HUGE_VAL;
    }

  if (loopSize <= 30)
    {
      energy = g_hairpinLoop[loopSize - 1];
      
    }
  else
    {
      energy = g_hairpinLoop[29] + g_misc[12] * log((double) loopSize / 30);
      }

  if (loopSize > 3)
    {
      energy += g_tstackh[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
      
    }
  else
    {
      energy += ct_auPenalty(i, j);
      }

  if (loopSize == 3)
    {
      struct triloop* loop;
      if (numTriloops)
	if ((loop = (triloop*)bsearch(g_seq + i - 1, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
	  {
	    energy += loop->energy;
	    
	  }
    }
  else if (loopSize == 4)
    {
      struct tloop* loop;
      if (numTloops)
	if ((loop = (tloop*)bsearch(g_seq + i - 1, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
	  {
	    energy += loop->energy;
	    
	  }
    }
  else if (loopSize == 6)
    {
      struct hexaloop* loop;
      if (numHexaloops)
	if ((loop = (hexaloop*)bsearch(g_seq + i - 1, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
	  {
	    energy += loop->energy;
	    
	  }
    }

  /* GGG */
  if (i >= 3 && g_seq[i - 3] == 2 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[j - 1] == 3)
    {
      energy += g_misc[8];
      
    }

  /* poly-C */
  if (loopSize == 3 && g_seq[i] == 1 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1)
    {
      energy += g_misc[11];
      
    }
  else
    {
      for (k = 0; k < loopSize; ++k)
	if (g_seq[i + k] != 1)
	  {
	     return energy;
	  }
      energy += g_misc[9] * loopSize + g_misc[10];
      
    }

  
  return energy;
}

double CtEnergy::ct_Es(int i, int j)
{
  if (i >= j)
    fputs("Error in Es(): i isn't less than j\n", stderr);

  return g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
}

double CtEnergy::ct_Ebi(int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;

  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
  if (ii >= jj)
    fputs("Error in Ebi(): jj isn't greater than ii\n", stderr);

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;

  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Ebi() called with nonsense\n", stderr);
      return 1.0;
    }
  else if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	{
	  loopEnergy = g_bulgeLoop[0] + g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[ii - 1]][g_seq[jj - 1]];
	  
	}
      else if (loopSize2 <= 30)
	{
	  loopEnergy = g_bulgeLoop[loopSize2 - 1] + ct_auPenalty(i, j) + ct_auPenalty(ii, jj);
	  
	}
      else
	{
	  loopEnergy = g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + ct_auPenalty(i, j) + ct_auPenalty(ii, jj);
	  
	}
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	{
	  loopEnergy = g_bulgeLoop[0] + g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[ii - 1]][g_seq[jj - 1]];
	  
	}
      else if (loopSize1 <= 30)
	{
	  loopEnergy = g_bulgeLoop[loopSize1 - 1] + ct_auPenalty(i, j) + ct_auPenalty(ii, jj);
	  
	}
      else
	{
	  loopEnergy = g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + ct_auPenalty(i, j) + ct_auPenalty(ii, jj);
	  
	}
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    {
      loopEnergy = g_sint2[basePairIndex(g_seq[i - 1], g_seq[j - 1])][basePairIndex(g_seq[ii - 1], g_seq[jj - 1])][g_seq[i]][g_seq[j - 2]];
      
    }
  else if (loopSize1 == 1 && loopSize2 == 2)
    {
      loopEnergy = g_asint1x2[basePairIndex(g_seq[i - 1], g_seq[j - 1])][basePairIndex(g_seq[ii - 1], g_seq[jj - 1])][g_seq[i]][g_seq[j - 2]][g_seq[j - 3]];
      
    }
  else if (loopSize1 == 2 && loopSize2 == 1)
    {
      loopEnergy = g_asint1x2[basePairIndex(g_seq[jj - 1], g_seq[ii - 1])][basePairIndex(g_seq[j - 1], g_seq[i - 1])][g_seq[jj]][g_seq[ii - 2]][g_seq[ii - 3]];
      
    }
  else if (loopSize1 == 2 && loopSize2 == 2)
    {
      loopEnergy = g_sint4[basePairIndex(g_seq[i - 1], g_seq[j - 1])][basePairIndex(g_seq[ii - 1], g_seq[jj - 1])][g_seq[i]][g_seq[j - 2]][g_seq[i + 1]][g_seq[j - 3]];
      
    }
  else if ((loopSize1 == 2 && loopSize2 == 3) ||
	   (loopSize1 == 3 && loopSize2 == 2))
    {
      loopEnergy = g_tstacki23[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]] +
	g_tstacki23[g_seq[jj - 1]][g_seq[ii - 1]][g_seq[jj]][g_seq[ii - 2]];
      
    }
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	{
	  loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
	  
	}
      else
	{
	  loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
	 
	}
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[g_seq[i - 1]][g_seq[j - 1]][0][0];
	  loopEnergy += g_tstacki[g_seq[jj - 1]][g_seq[ii - 1]][0][0];
	  
	}
      else
	{
	  loopEnergy += g_tstacki[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
	  loopEnergy += g_tstacki[g_seq[jj - 1]][g_seq[ii - 1]][g_seq[jj]][g_seq[ii - 2]];
	  
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
		asPenalty = g_misc[4];
      loopEnergy += asPenalty;
      
    }
  
  
  return loopEnergy;
}

double CtEnergy::ct_Ee(int i, int j, int ds)
{
  int count;
  double energy;
  struct stack_node* new_top;

  energy = g_misc[5] + ct_auPenalty(i, j) + (ct_isHomodimer() ? g_misc[12] / 1.75 * log(2.0) : 0.0);

  if (ds)
    {
      if (g_hasStackingInfo)
	energy += ct_tstackOrDangle(i, j, 1);
      else
	energy += ct_chooseDangle(stack->j, j);
      for (new_top = stack, count = 0; count < ds - 1; new_top = new_top->next, ++count)
	{
	  if (g_hasStackingInfo)
	    energy += ct_tstackOrDangle(new_top->j, new_top->i, 1);
	  else
	    energy += ct_chooseDangle(new_top->next->j, new_top->i);
	  energy += ct_auPenalty(new_top->i, new_top->j);
	  
	}
      energy += ct_auPenalty(new_top->i, new_top->j);
      
      if (g_hasStackingInfo)
	energy += ct_tstackOrDangle(new_top->j, new_top->i, 1);
      else
	energy += ct_chooseDangle(i, new_top->i);
    }
  else if (!g_nodangle)
    {
      if (g_hasStackingInfo)
	energy += ct_tstackOrDangle(i, j, 1);
      else
	{
	  if (g_next[i - 1] && g_prev[i] && (ct_Ed3(i, j, i + 1) <= 0.0))
	    {
	     
	      energy += ct_Ed3(i, j, i + 1);
	    }
	  if (g_prev[j - 1] && g_next[j - 2] && (ct_Ed5(i, j, j - 1) <= 0.0))
	    {
	    
	      energy += ct_Ed5(i, j, j - 1);
	    }
	}
    }

 
  return energy;
}

double CtEnergy::ct_auPenalty(int i, int j)
{
  if (basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 0 ||
      basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 3 ||
      basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 4 ||
      basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 5)
    return g_misc[6];
  return 0;
}

double CtEnergy::ct_chooseDangle(int a, int b)
{
  double energy;

  if (b == a + 1 || g_nodangle)
    return 0;
  else if (b == a + 2)
    {
      double d5, d3;
      d5 = ct_Ed5(b, g_bp[b - 1], b - 1);
      d3 = ct_Ed3(g_bp[a - 1], a, a + 1);

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
      if (ct_Ed3(g_bp[a - 1], a, a + 1) <= 0.0)
	{
	  energy += ct_Ed3(g_bp[a - 1], a, a + 1);
	 }
      if (ct_Ed5(b, g_bp[b - 1], b - 1) <= 0.0)
	{
	  energy += ct_Ed5(b, g_bp[b - 1], b - 1);
	}
      return energy;
    }
}

double CtEnergy::ct_Ed5(int i, int j, int k)
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

double CtEnergy::ct_Ed3(int i, int j, int k)
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

double CtEnergy::ct_tstackOrDangle(int i, int j, int external)
{
  if (j > 1 &&
	   g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	   g_next[i - 1] == i + 1 && g_prev[i] == i &&
	   g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	   g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	   g_bp[j - 2] == 0 && g_bp[i] == 0)
    {
      
      return external ? ct_Etstacke(i, j) : ct_Etstackm(i, j);
    }
  else if (j > 1 &&
	   g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	   g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	   g_bp[j - 2] == 0)
    {
      
      return ct_Ed5(i, j, j - 1);
    }
  else if (g_next[i - 1] == i + 1 && g_prev[i] == i &&
	   g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	   g_bp[i] == 0)
    {
      
      return ct_Ed3(i, j, i + 1);
    }
  return 0.0;
}

int CtEnergy::ct_isHomodimer()
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

int CtEnergy::ct_isCircular()
{
  return g_prev[0] == g_len && g_next[g_len - 1] % g_len == 1;
}
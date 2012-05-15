#include "util.h"

int util::roundInt(double d)
{
  return (int) (d + .5);
}

void util::strcatc(char* str, char c)
{
  str[strlen(str) + 1] = 0;
  str[strlen(str)] = c;
}

void util::checkArray(char** array, unsigned int* available, unsigned int used, unsigned int increment)
{
  if (used == *available)
    {
      *array = (char*)  xrealloc(*array, *available + increment);
      *available += increment;
    }
}

unsigned char util::toNum(char c)
{
  c = toupper(c);
  switch (c)
    {
    case 'A': case '0':
      return 0;
    case 'C': case '1':
      return 1;
    case 'G': case '2':
      return 2;
    case 'T': case 'U': case '3':
      return 3;
    }
  return 4;
}

int util::seqcmp(unsigned char* seq1, unsigned char* seq2, int length)
{
  int i;

  for (i = 0; i < length; ++i)
    if (seq1[i] < seq2[i])
      return -1;
    else if (seq1[i] > seq2[i])
      return 1;
  return 0;
}

int util::same(unsigned char* a, unsigned char* b, int len)
{
  int i;

  for (i = 1; i <= len; ++i)
    if (a[i] != b[i])
      return 0;
  return 1;
}
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cstring>
#include <string>
#include <time.h>
#include <vector>
#include <omp.h>
#include <sstream>
using namespace std;

//****U S E R S  S H O U L D  E N T E R   T H E S E   P A R A M E T E R S**************************//
const int SeqLen=20; //LENGTH OF THE SEQUENCE IDENTIFIER (BARCODE)
const unsigned long int DesiredNofBarcodes=1000; //NUMBER OF BARCODES NEEDED
string LeftAdaptor="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"; //THE ADAPTOR SEQUENCE ADDED TO THE LEFT OF THE BARCODE (5')
string RightAdaptor=""; //THE ADAPTOR SEQUENCE ADDED TO THE RIGHT OF THE BARCODE (3')
const int homoplimit=4; //MAXIMUM NUMBER OF MONO, DI OR TRI-MERS ALLOWED WITHIN A BARCODE
const double MinGC=0.45; // MINIMUM GC PERCENTAGE OF THE BARCODE
const double MaxGC=0.65; //MAXIMUM GC PERCENTAGE OF THE BARCODE
int LenDiffThreshold; 
const double SelfHybT=50.0; //SELF HYBRIDIZATION TM THRESHOLD
double EditDistanceThreshold=4;
double EditDistanceThreshold_Self=5;
const double Hyb_Temperature=50;
//**************************************************************************************************
//***********C O N S T A N T S******************
const double R = 0.0019872; //For Tm calculations
//**********************************************
const int N_THREADS=7;
#include <pthread.h>
pthread_mutex_t poolMutex;
std::vector<string> pool;
bool die = false;
//**D E P E N D E N C I E S*********************
#include "LZW.h"
#include "hybrid-ss-min.h"
#include "GenerateRandomDNASequences_v6.h"
#include "FindIndepentSet_v2.h"
#include "DegreeDistribution_v1.h"
//**********************************************

//***************I/O FILE NAMES****************************************************
string PutSeqFNS="RandomSeqs_20mers";
//************************************************************************

void GenerateFileNames(string, string&, bool,unsigned long int);

int main(){

	string temp,seq, FN,FN1,FN2,s0;
	GenerateSequences Sequences;
	Network NW;
	DegreeDistribution DegreeDist;
	char PutSeqFNC[150],UniqueBarcodeFNC[150],EDDistFNC[150],DDFNC[150];


	unsigned long int i=0,j=0,SIZE, IndependentSetSize; // Keep number of barcodes
	bool keepbarcode;
	vector <string> randseqs;
	string randseq_revcomp;

	Sequences.initialisevars();
	Sequences.DetermineFilteringThresholds();
	//Startup generating threads!
	Sequences.GenerateRandomSequence_SetGC();
	//Create local buffer.
	vector<string> sequenceVector;
	do{
		//Wait for pool!
		pthread_mutex_lock(&poolMutex);
		//Check for size of pool
		if (pool.size() >= 100) {
			for (int i=0;i<pool.size();++i) {
			  sequenceVector.push_back(pool[i]);
			}
			pool.clear();
		}
		pthread_mutex_unlock(&poolMutex);
		if (sequenceVector.empty()) {
		  	//Do something silly for some time
			clock_t goal = 10 + clock();
			while (goal > clock());
			continue;
		}
		i=0;
		while (i<sequenceVector.size()) {
		  NW.RetrieveUniqueNodes(sequenceVector[i]);
		  if (NW.CommonSet.size()>DesiredNofBarcodes)
		    break;
		  ++i;
		}
		sequenceVector.clear();

		if(NW.CommonSet.size()%100==0)
			cout << NW.CommonSet.size() << "   Barcodes Selected" << endl;

	}while(NW.CommonSet.size()<=DesiredNofBarcodes);
	
	//Cleanup generating threads!
	pthread_mutex_lock(&poolMutex);
	die = true;
	cout << pool.size() << endl;
	pool.clear();
	pthread_mutex_unlock(&poolMutex);

	FN.append(PutSeqFNS);
	FN.append("_UniqueBarcodes.txt");
	strcpy(UniqueBarcodeFNC,FN.c_str());

	NW.PrintUniquePutBarcodes(UniqueBarcodeFNC);

	FN1.append(PutSeqFNS);
	FN1.append("_DegreeDistribution.txt");
	strcpy(DDFNC,FN1.c_str());
	

	FN2.append(PutSeqFNS);
	FN2.append("_EDDistribution.txt");
	strcpy(EDDistFNC,FN2.c_str());

	DegreeDist.NofBarcodes=NW.CommonSet.size();
	DegreeDist.initialisevars();
	DegreeDist.GenerateNetwork(NW.CommonSet);
	DegreeDist.DegreeDist(DDFNC,EDDistFNC,NW.CommonSet);


	return 0;
}

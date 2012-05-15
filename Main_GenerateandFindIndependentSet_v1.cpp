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
const unsigned long int DesiredNofBarcodes=135000; //NUMBER OF BARCODES NEEDED
string LeftAdaptor="ACACTCTTTCCCTACACGACGCTCTTCCGATCT"; //THE ADAPTOR SEQUENCE ADDED TO THE LEFT OF THE BARCODE (5')
string RightAdaptor=""; //THE ADAPTOR SEQUENCE ADDED TO THE RIGHT OF THE BARCODE (3')
const int homoplimit=4; //MAXIMUM NUMBER OF MONO, DI OR TRI-MERS ALLOWED WITHIN A BARCODE
const double MinGC=0.45; // MINIMUM GC PERCENTAGE OF THE BARCODE
const double MaxGC=0.65; //MAXIMUM GC PERCENTAGE OF THE BARCODE
int LenDiffThreshold; 
const double SelfHybT=50.0; //SELF HYBRIDIZATION TM THRESHOLD
double EditDistanceThreshold=5;
double EditDistanceThreshold_Self=5;
const double Hyb_Temperature=50;
//**************************************************************************************************
//***********C O N S T A N T S******************
const double R = 0.0019872; //For Tm calculations
//**********************************************
const int N_THREADS=7;
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
	
	vector <bool> passedfilter;
//	HybridSSMin* ceva=new HybridSSMin();

	Sequences.initialisevars();
	Sequences.DetermineFilteringThresholds();
	do{
		passedfilter.resize(N_THREADS); //Initialise randseq filter vector
		randseqs.resize(N_THREADS); // Initialise randseq vector
		randseq_revcomp.clear(); //Initialise randseq_revcomp
		Sequences.GenerateRandomSequence_SetGC(passedfilter,randseqs); // Generate N_THREADS x random sequences
		for(i=0;i<N_THREADS;i++){
			if(passedfilter[i])
				NW.RetrieveUniqueNodes(randseqs[i]);
		}
//		if(NW.CommonSet.size()%100==0)
//			cout << NW.CommonSet.size() << "   Barcodes Selected" << endl;
	}while(NW.CommonSet.size()<=DesiredNofBarcodes);

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

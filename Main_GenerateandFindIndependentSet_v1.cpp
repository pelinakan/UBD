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
#include <map>
#include <omp.h>
#include <sstream>

#include "getopt.h"
using namespace std;

//****U S E R S  S H O U L D  E N T E R   T H E S E   P A R A M E T E R S**************************//
//Value dictionary
typedef struct p{
	int type;
	void *var;
	p() {
		type = -1;
		var = NULL;
	};

	p (int type_p, void* var_p) {
		type = type_p;
		var = var_p;
	};
}parserEntry;

map<string,parserEntry*> parseMap;
//////////////////////////////////////

int SeqLen; //LENGTH OF THE SEQUENCE IDENTIFIER (BARCODE)
unsigned long int DesiredNofBarcodes; //NUMBER OF BARCODES NEEDED
string LeftAdaptor; //THE ADAPTOR SEQUENCE ADDED TO THE LEFT OF THE BARCODE (5')
string RightAdaptor; //THE ADAPTOR SEQUENCE ADDED TO THE RIGHT OF THE BARCODE (3')
int homoplimit; //MAXIMUM NUMBER OF MONO, DI OR TRI-MERS ALLOWED WITHIN A BARCODE
double MinGC; // MINIMUM GC PERCENTAGE OF THE BARCODE
double MaxGC; //MAXIMUM GC PERCENTAGE OF THE BARCODE
int LenDiffThreshold; 
double SelfHybT; //SELF HYBRIDIZATION TM THRESHOLD
double EditDistanceThreshold;
double EditDistanceThreshold_Self;
double Hyb_Temperature;
vector<string> AdaptorList;
vector<string> RestrictionList;
//**************************************************************************************************
//***********C O N S T A N T S******************
const double R = 0.0019872; //For Tm calculations
//**********************************************
const int N_THREADS=6;
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
string PutSeqFNS="";
//************************************************************************

/**
 * Print usage instructions
 */
static int print_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Version: 1.0\n");
  fprintf(stderr, "Usage:   designBarcodes [options] <output>\n");
  fprintf(stderr, "Options: \n");
  fprintf(stderr, "         -c [STR]            Configuration file\n");
  fprintf(stderr, "         -n [INT]            Number of barcodes to output.\n");
  return 1;
}

void buildParserMap() 
{
	//Integers
	parseMap["desiredBarcodeLength"] = new parserEntry(1,&SeqLen);
	parseMap["homopolymerLimit"] = new parserEntry(1,&homoplimit);

	//Strings
	parseMap["leftAdaptor"] = new parserEntry(2,&LeftAdaptor);
	parseMap["rightAdaptor"] = new parserEntry(2,&RightAdaptor);

	//Doubles
	parseMap["minGC"] = new parserEntry(3,&MinGC);
	parseMap["maxGC"] = new parserEntry(3,&MaxGC);
	parseMap["selfHybridizationTemperature"] = new parserEntry(3,&SelfHybT);
	parseMap["hybridizationTemperature"] = new parserEntry(3,&Hyb_Temperature);
	parseMap["minEditDistance"] = new parserEntry(3,&EditDistanceThreshold);
	parseMap["minEditDistanceSelf"] = new parserEntry(3,&EditDistanceThreshold_Self);
}

bool parseConfig(string path)
{
	FILE* cFile = fopen(path.c_str(),"r");
	if (cFile == NULL)
		return false;

	vector<string> * listRef = NULL;
	char line[1000];
	while (fgets(line,1000,cFile) != NULL) {
		if ((line[0] == '#') || (strlen(line) < 2))
			continue;
		string l = line;
		if (l.find("</") != string::npos) {//End of a definition
			listRef = NULL;
			continue;
		}
		if (listRef) {
			//Clean-up string
			l.erase(l.find_last_not_of(" \n\r\t")+1);
			listRef->push_back(l);
			continue;
		}
		if (l[0] == '<') {//Beginning of a definition list
			if (l.find("handles") != string::npos) {
				listRef = &AdaptorList;
			} else if (l.find("restrictionSites") != string::npos) {
				listRef = &RestrictionList;
			} else {
				//Now a known list!
				fprintf(stdout,"Wrong enumeration %s\n",l.c_str());
				continue;
			}
		}
		if (l.find("=") == string::npos)
			continue;
		string key = l.substr(0,l.find("="));
		key.erase(key.find_last_not_of(" \n\r\t")+1);
		string value = l.substr(l.find("=")+1,l.size()-l.find("="));
		value.erase(value.find_last_not_of(" \n\r\t")+1);
		//Check against known keys and populate the relevant values
		if (parseMap.find(key) != parseMap.end()) {
			switch (parseMap[key]->type) {
			case 1://Integer
				*static_cast<int*>(parseMap[key]->var) = atoi(value.c_str());
				break;
			case 2://String
				*static_cast<string*>(parseMap[key]->var) = value;
				break;
			case 3://Double
				*static_cast<double*>(parseMap[key]->var) = atof(value.c_str());
				break;
			default:
				fprintf(stderr,"read wrong variable type\n");
				break;
			}
		} else {
			fprintf(stderr,"Unable to parse %s variable\n",key.c_str());
		}
	}
}

int main(int argc, char *argv[]){
	//Build parsing map
	buildParserMap();
	//TODO: don't leak the parser map

	string configPath = "";
	//Parse arguments
	int arg;
	while ((arg = getopt(argc,argv,"c:n:")) >= 0) {
		switch (arg) {
		case 'c':
			configPath = optarg;
			break;
		case 'n':
			DesiredNofBarcodes = atoi(optarg);
			break;
		default:
			fprintf(stderr,"Unknown argument %c\n",arg);
			break;
		}
	}

	if (optind + 1 != argc) {//Something is wrong with the parameters
		print_usage();
		return -1;
	}
	//Output file name
	PutSeqFNS = argv[optind];

	//Parse config file
	parseConfig(configPath);

	string temp,seq, FN,FN1,FN2,s0;
	GenerateSequences Sequences;
	Network NW;
	DegreeDistribution DegreeDist;
	char /*PutSeqFNC[150],*/UniqueBarcodeFNC[150],EDDistFNC[150],DDFNC[150];


	unsigned long int i=0; // Keep number of barcodes
	//bool keepbarcode;
	vector <string> randseqs;
	string randseq_revcomp;
	
	omp_set_num_threads(N_THREADS);
	pthread_mutex_init(&poolMutex,NULL);

	unsigned long int old_count = 0;
	unsigned long int same_count = 0;
	//Sequences.initialisevars();
	Sequences.DetermineFilteringThresholds();
	//Startup generating threads!
	Sequences.GenerateRandomSequence_SetGC();
	//Create local buffer.
	vector<string> sequenceVector;
	do{
		//Wait for pool!
		pthread_mutex_lock(&poolMutex);
		//Check for size of pool
		if (pool.size() >= 50) {
			for (int i=0;i<pool.size();++i) {
			  sequenceVector.push_back(pool[i]);
			}
			pool.clear();
		}
		pthread_mutex_unlock(&poolMutex);
		if (sequenceVector.empty()) {
		  //Do something silly for some time
		  clock_t goal = 500 + clock();
		  while (goal > clock());
		  continue;
		}
		i=0;
		while (i < sequenceVector.size()) {
		  NW.RetrieveUniqueNodes(sequenceVector[i]);
		  if (NW.CommonSet.size() > DesiredNofBarcodes)
		    break;
		  ++i;
		}
		sequenceVector.clear();
		/*if (old_count != NW.CommonSet.size()) {
		  old_count = NW.CommonSet.size();
		  same_count = 0;
		} else {
		  ++same_count;
		  if (same_count > 100)
		    fprintf(stdout,"Found solution with %ld\n",old_count);
		    }*/
		
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
	if (!DegreeDist.initialisevars()) {
		fprintf(stderr,"Unable to allocate memory for degree distribution computation\n. Skipping...\n");
		return 0;
	}
	DegreeDist.GenerateNetwork(NW.CommonSet);
	DegreeDist.DegreeDist(DDFNC,EDDistFNC/*,NW.CommonSet*/);

	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <time.h>
#ifdef __APPLE__
#include <tr1/unordered_map>
using namespace tr1;
#else
#include <unordered_map>
#endif
#include <pthread.h>
#include "pol_fastq.h"

using namespace std;

static int print_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: findIndexes \n");
	fprintf(stderr, "Contact: Paul Costea <paul.igor.costea@embl.de>\n\n");
	fprintf(stderr, "Usage:   findIndexes [options] <ids.txt> <in.fastq> <out.fastq>\n\n");
	fprintf(stderr, "Options: \n");
	fprintf(stderr, "         -m INT     allowed mismatches [2]\n");
	fprintf(stderr, "         -k INT     kMer length [1/3*length]\n");
	fprintf(stderr, "         -s INT     start position of ID [0]\n");
	fprintf(stderr, "         -l INT     length of ID [0]\n");
	fprintf(stderr, "         -e INT     id positional error [0]\n");
	fprintf(stderr, "         -p STRING  name of pair file [NULL]\n");
	fprintf(stderr, "\n");
	return 1;
}

/*
 * Structure for ID definition. May be augmented with additional information.
 */
class IdStruct {
public:
	IdStruct() {
	};

	IdStruct(const IdStruct &other) {
		id = other.id;
	};

	string id;

};

long int count = 0;
long int perfectMatch = 0;
long int totalReads = 0;

long int ambiguous = 0;
long int editDistanceTooBig = 0;

//Parameter
int amismatch = 2;
int kLen = 0;
int probeStartPos = 0;
int probeLength = 0;
int positionalError = 0;

int qualMin = 1000;
int qualMax = 0;

std::string pairFile="";

pthread_mutex_t readMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t writeMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t countersMutex = PTHREAD_MUTEX_INITIALIZER;

#define GUARDED_INC(x) \
		pthread_mutex_lock(&countersMutex); \
		++x; \
		pthread_mutex_unlock(&countersMutex);

/////////////////Hash declaration//////////////////////
unordered_map<std::string,int> mainIds_map;
unordered_map<std::string,std::list<IdStruct>> kIds_map;
///////////////////////////////////////////////////////

bool loadIds(FILE* snpFile)
{
	char line[100000];
	char *tok = new char[10000];
	for (int i=0; i<10000; ++i)
		tok[i]='\0';
	while (fgets(line,100000,snpFile)) {
		//Replace new line with \0
		line[strlen(line)-1] = '\0';
		const char *t = line;
		int pos = 0;
		string id = "";
		while (*t) {
			t = pol_util::toksplit(t, '\t', tok, 10000);
			if (pos == 0) {//ID
				id = tok;
			}
			++pos;
		}
		mainIds_map[id] = 0;
	}
	delete[] tok;

	return true;
}

/*
 * Edit distance computation
 */
int editDistance(unsigned int **d, std::string s1, std::string s2, std::string qual="", int* alignScore=NULL, int maxDist=1000)
{
	int xLen = s1.size();
	int yLen = s2.size();

	if (d == NULL) {
		//The matrix given has not been initialized
		fprintf(stdout,"Null computation matrix!\n");
		return -1;
	}

	for(int x = 1; x <= xLen; ++x)
		for(int y = 1; y <= yLen; ++y) {
			d[x][y] = std::min( std::min(d[x - 1][y] + 1,d[x][y - 1] + 1),d[x - 1][y - 1] + (s1[x - 1] == s2[y - 1] ? 0 : 1) );
		}

	//Find min for sub-global alignment
	int min = 1000;
	int iPos = 0;
	for (int i=xLen;i>1/*xLen-extraBases-2*/;--i) {
		if (d[i][yLen] < min) {
			min = d[i][yLen];
			iPos = i;
		}
	}

	if ((alignScore != NULL) && (min <= maxDist)){//Do traceback
		//Compute score!
		*alignScore = 0;

		int i = iPos;
		int j = yLen;

		while ((i > 0) && (j > 0))
		{
			int Score = d[i][j];
			int ScoreDiag = d[i - 1][j - 1];
			int ScoreUp = d[i][j - 1];
			int ScoreLeft = d[i - 1][j];
			if (Score == ScoreDiag + (s1[i-1] == s2[j-1] ? 0 : 1))
			{
				if (s1[i-1] != s2[j-1])
					*alignScore -= qual[i-1];
				else
					*alignScore += qual[i-1];
				i = i - 1;
				j = j - 1;
			}
			else if (Score == ScoreLeft + 1)
			{
				//Add deletion score
				*alignScore -= 100;
				i = i - 1;
			}
			else if (Score == ScoreUp + 1)
			{
				//Add insertions score
				*alignScore -= qual[j-1];
				j = j - 1;
			}
		}
		//Now, add deletions in the beginning to the score!
		*alignScore -= 100 * j;
	}

	return min;
}

void createKmerMap(int kLen)
{
	unordered_map<std::string,int>::const_iterator it;
	for (it = mainIds_map.begin(); it != mainIds_map.end(); ++it) {//For each, introduce only part of it into the new map!
		for (unsigned int i=0; i<=(it->first.size()-kLen); ++i) {
			std::string newKey = it->first.substr(i,kLen);
			IdStruct idS;
			idS.id = it->first;
			kIds_map[newKey].push_back(idS);
		}
	}

}

typedef struct {
  FILE *inFile,*outFile,*pairFile,*pairFileOut;
}threadArgs;

/**
 * Thread entry
 */
void* idSearch(void* arguments)
{
	unordered_map<std::string,int>::iterator it;
	unordered_map<std::string,int> kMerHitMap;

	//Get longer sequence!!!!
	int leftShift = (probeStartPos-positionalError >= 0)? positionalError : 0;
	int rightShift = positionalError;
	int maxKmerHits = (probeLength+rightShift) + (leftShift);
	std::list<std::string>* orderedIds = new std::list<std::string>[maxKmerHits+1];

	//create matrix for alignment
	int xLen = probeLength+rightShift+leftShift;
	int yLen = probeLength;

	unsigned int **d = new unsigned int*[xLen+1];
	for (int i=0;i<=xLen;++i) {
		d[i] = new unsigned int[yLen+1];
	}
	d[0][0] = 0;
	for(int x = 1; x <= xLen; ++x) d[x][0] = 0;
	for(int x = 1; x <= yLen; ++x) d[0][x] = x;

	pol_util::FastqEntry* e,*e1;
	e1=NULL;
	threadArgs args = *(threadArgs*)arguments;
	while (true) {
		//Lock file mutex
		pthread_mutex_lock(&readMutex);
		e = pol_util::FastqEntry::readEntry(args.inFile);
		if (args.pairFile!=NULL)
		  e1 = pol_util::FastqEntry::readEntry(args.pairFile);
		pthread_mutex_unlock(&readMutex);
		if (e == NULL)
			break;
		int pos = probeStartPos+probeLength;
		bool has_seq = true;
		if (has_seq) {
			GUARDED_INC(totalReads)
			std::string id = e->getSequence(pos-probeLength,pos);
			std::string qual = e->getQuality(pos-probeLength,pos);
			//Search ID in map!
			it = mainIds_map.find(id);
			if (it != mainIds_map.end()) {//Found!
				GUARDED_INC(count)
				GUARDED_INC(perfectMatch)
				//Print this to output map!
				if (args.outFile != NULL) {
					e->setOptional("barcode="+id);
					pthread_mutex_lock(&writeMutex);
					e->write(args.outFile);
					if (e1!=NULL)
					  e1->write(args.pairFileOut);
					pthread_mutex_unlock(&writeMutex);
				}
			} else {//Probably with mismatch!
				//////////////////////////////////////////////////////////////////////////
				////////////////////////////Do kmer search////////////////////////////////
				id = e->getSequence(probeStartPos-leftShift,probeStartPos+probeLength+rightShift);
				qual = e->getQuality(probeStartPos-leftShift,probeStartPos+probeLength+rightShift);

				unordered_map<std::string,std::list<IdStruct>>::iterator splitIt;
				std::list<IdStruct>::const_iterator it;
				unordered_map<std::string,int>::iterator kIt;

				for (int displace=0; displace <= leftShift; ++displace) {
					int x = displace;
					bool forced = false;
					while ( (x <= id.size()-kLen) && (!forced)) {
						std::string id1 = id.substr(x,kLen);
						splitIt = kIds_map.find(id1);
						if (splitIt != kIds_map.end()) {
							//Add to counted map
							std::list<IdStruct>* l = &(splitIt->second);
							for (it = l->begin(); it != l->end(); ++it) {
								std::string locId = it->id;
								kIt = kMerHitMap.find(locId);
								if (kIt != kMerHitMap.end()) {
									kMerHitMap[locId]+=1;
								} else {//Add
									kMerHitMap[locId]=1;
								}
							}
						}
						//Make sure we get the last overlapping k-mer
						x+=kLen/2;

						if (x > id.size()-kLen) {
							x = id.size()-kLen;
							forced = true;
						}
					}
				}

				///////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////
				bool goodHit = false;
				int min_ed = 1000;
				int max_score = 0;
				std::string found_id = "";
				//Go through map and find maximum.
				//Order hits
				for (kIt = kMerHitMap.begin(); kIt != kMerHitMap.end(); ++kIt) {
					orderedIds[kIt->second].push_back(kIt->first);
				}
				int goOn = 0;
				int score = 0;
				int searching = 0;
				//Take only first x most kMer containing Id's
				for (int i=maxKmerHits-1;i>=0;--i) {
					if (orderedIds[i].size() != 0 && (goOn <= 4) ) {//This is the list of most hits
						for (std::list<std::string>::const_iterator it=orderedIds[i].begin(); it != orderedIds[i].end(); ++it) {
							double ed = editDistance(d,id,(*it),qual,&score,amismatch);
							++searching;
							if (ed < min_ed) {//New min
								min_ed = ed;
								max_score = score;
								goodHit = true;
								found_id = (*it);
							} else if (ed == min_ed) {//Two with same min :(
								//Discriminate by score!
								if (score == max_score) {
									goodHit = false;
									goOn = 0;
								} else if (score > max_score) {
									goodHit = true;
									max_score = score;
									found_id = (*it);
								}
							}
						}
						//++goOn;
					}
					orderedIds[i].clear();
				}

				if (goodHit && min_ed<= amismatch) {
					//We have a "best" hit!
					GUARDED_INC(count)
					if (args.outFile != NULL) {
						e->setOptional("barcode="+found_id);
						pthread_mutex_lock(&writeMutex);
						e->write(args.outFile);
						if (e1!=NULL)
						  e1->write(args.pairFileOut);
						pthread_mutex_unlock(&writeMutex);
					}
				} else {
					//Still print this!
					pthread_mutex_lock(&writeMutex);
					e->write(args.outFile);
					if (e1!=NULL)
					  e1->write(args.pairFileOut);
					pthread_mutex_unlock(&writeMutex);
					if (!goodHit) {
						GUARDED_INC(ambiguous)
					} else {
						GUARDED_INC(editDistanceTooBig)
					}
				}

				kMerHitMap.clear();

			}
		}
		delete e;
		if (e1 != NULL)
		  delete e1;
	}

	delete[] orderedIds;
	for (int i=0;i<=xLen;++i) {
		delete[] d[i];
	}
	delete[] d;

	return NULL;
}

/**
 * Main of app
 */
int main(int argc, char *argv[])
{

	FILE *out = NULL;
	int arg;
	//Get args
	while ((arg = getopt(argc, argv, "m:k:s:l:e:p:")) >= 0) {
		switch (arg) {
		case 'm': amismatch = atoi(optarg); break;
		case 'k': kLen = atoi(optarg);
		break;
		case 's': probeStartPos = atoi(optarg); break;
		case 'l':
			probeLength = atoi(optarg);
			if (kLen == 0) {//If not otherwise set
				kLen = probeLength / 3;
			}
			break;
		case 'e': positionalError = atoi(optarg); break;
		case 'p': pairFile = optarg; break;
		default:
			fprintf(stderr,"Read wrong arguments! \n");
			break;
		}
	}

	if (argc-optind != 3) {
		//Not enough paramters!
		print_usage();
		return 1;
	}
	if (probeLength <= 0) {
		fprintf(stderr,"You must input the length of the ID\n");
		return 1;
	}
	//Seed rand
	srand((unsigned)time(0));
	//Open files!
	FILE *ids,*inF,*pairF,*pairFOut;
	pairF = NULL;
	pairFOut = NULL;
	if ((ids = fopen(argv[optind], "r")) == 0) {
		fprintf(stderr, "Failed to open file %s\n", argv[optind]);
		return 1;
	}
	if ((inF = fopen(argv[optind+1], "r")) == 0) {
		fprintf(stderr, "Failed to open file %s\n", argv[optind+1]);
		return 1;
	}
	if (!pairFile.empty()) {
	  if ((pairF = fopen(pairFile.c_str(),"r")) == 0) {
	    fprintf(stderr, "Faile to open pair file %s\n", pairFile.c_str());
	    return 1;
	  }
	  std::string outFile = argv[optind+2];
	  outFile += "_pair";
	  if ((pairFOut = fopen(outFile.c_str(), "w"))==0) {
	    fprintf(stderr, "Failed to open pair output file: %s\n",outFile.c_str());
	    return 1;
	  }
	}
	if ((out = fopen(argv[optind+2], "w")) == 0) {
		fprintf(stderr, "Failed to open for writing: %s\n", argv[optind+2]);
		return 1;
	}
	fprintf(stdout,"Printing map to %s\n",argv[optind+2]);

	//Load IDS!
	loadIds(ids);
	//Create kMer map
	createKmerMap(kLen);
	pthread_t threadList[8];
	threadArgs args;
	args.inFile = inF;
	args.outFile = out;
	args.pairFile = pairF;
	args.pairFileOut = pairFOut;

	/*DO THREAD STUFF*/
	for (int i=0; i<8; ++i) {
		pthread_create( &threadList[i], NULL, idSearch, (void*) &args);
	}
	//Now join them...
	for (int i=0; i<8; ++i) {
		pthread_join(threadList[i],NULL);
	}

	fprintf(stdout,"Found: %ld\n",count);
	fprintf(stdout, "Perfect match: %ld\n",perfectMatch);
	fprintf(stdout, "Ambiguous ID's: %ld\n",ambiguous);
	fprintf(stdout, "Edit distance exceeding limit: %ld\n",editDistanceTooBig);

	fprintf(stdout, "Total reads: %ld\n",totalReads);


	fclose(ids);
	fclose(inF);
	if (out != NULL)
		fclose(out);
	if (pairF != NULL)
	  fclose(pairF);
	if (pairFOut != NULL)
	  fclose(pairFOut);
	return 0;
}



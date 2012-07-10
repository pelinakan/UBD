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

class GenerateSequences{

public:
	//Constructor
	GenerateSequences();
	~GenerateSequences();
	
	vector< vector <double> > Features;
	string randomseq_setGC(int);
	bool checkforruns(string);
	string RevComp(string);
	int CalculateEditDistance(string,string);
	static string AppendAdaptors(string&);
	void DetermineFilteringThresholds(void);
	void GenerateRandomSequence_SetGC();

private:
	unsigned int** d;
};

GenerateSequences::GenerateSequences()
{
	d = new unsigned int*[SeqLen+1];
	for (int i=0;i<=SeqLen;++i) {
		d[i] = new unsigned int[SeqLen+1];
	}
	d[0][0] = 0;
	for(int x = 1; x <= SeqLen; ++x) d[x][0] = x;
	for(int x = 1; x <= SeqLen; ++x) d[0][x] = x;
}

GenerateSequences::~GenerateSequences()
{
	//Don't leak d
	for (int i=0;i<=SeqLen;++i) {
		delete[] d[i];
	}
	delete[] d;
}

string GenerateSequences::randomseq_setGC(int GCcontent){

   int i,x,noffilledpos=0;
   double slen;
   int y;
   bool* filled = new bool[SeqLen];
   string qual,randseq;
   randseq.resize(SeqLen);
   for(i=0;i<SeqLen;i++) filled[i]=0; // Initialise the array keeping the state of each position
   slen=double(SeqLen);
   noffilledpos=0;
   while(noffilledpos<GCcontent){ // The GCcontent number of positions with C or G
	   x = (int) (slen*rand()/(RAND_MAX+1.0)); //which position will be C or G
	   if(!(filled[x])){
		   filled[x]=1;
		   //y = (double) (rand()/(RAND_MAX+1.0)); //is it C or G
		   //if(y<0.5)
		   y = rand();
		   if (y % 2 == 0)
			   randseq[x]='C';
		   else
			   randseq[x]='G';
		   noffilledpos++;
	   }
   }
   for(i=0;i<SeqLen;i++){
	   if(!filled[i]){
	     //y = (double) (rand()/(RAND_MAX+1.0)); //is it C or G
	     //if(y<0.5)
	     y = rand();
	     if (y % 2 == 0)
	       randseq[i]='A';
	     else
	       randseq[i]='T';
	   }
   }
   delete[] filled;
   return randseq;
}

bool GenerateSequences::checkforruns(string Tag){
	
int i,j,l=0,runs;
bool selected=1;
l=Tag.length();
	selected=1;
	
	for(i=0;i<=l-homoplimit;i++){
		runs=0;
		for(j=1;j<homoplimit;j++){
			if(Tag[i]==Tag[i+j])
				runs++;
		}
		if (runs==(homoplimit-1)){
			selected=0;
			break;
		}
	}
// Runs of nucleotides checked
	
	char r[4];
	for(j=2;j<l-3;j+=2){
		if(selected==0)
			break;
		runs=0;
		r[0]=Tag[j-2]; r[1]=Tag[j-1]; r[2]='\0';
		if(r[0]==Tag[j] && r[1]==Tag[j+1])
			runs++;
		if(r[0]==Tag[j+2] && r[1]==Tag[j+3])
			runs++;
		if (runs==2){
			selected=0;
			break;
		}
	}
// runs of dimers >2 checked

	for(j=3;j<l-5;j+=3){
		if(selected==0)
			break;
		runs=0;
		r[0]=Tag[j-3]; r[1]=Tag[j-2]; r[2]=Tag[j-1];
		if(r[0]==Tag[j] && r[1]==Tag[j+1] && r[2]==Tag[j+2])
			runs++;
		if(r[0]==Tag[j+3] && r[1]==Tag[j+4] && r[2]==Tag[j+5])
			runs++;
		if (runs==2){
			selected=0;
			break;
		}
	}
// runs of trimers >2 checked
	return selected;
}

string GenerateSequences::RevComp(string s){	
	int z;
	string s_revcom;
	for(z=(s.length()-1);z>=0;--z){
		switch (s[z]){
			case 'A' : { s_revcom.append("T"); break; }
			case 'C' : { s_revcom.append("G"); break; }
			case 'G' : { s_revcom.append("C"); break; }
			case 'T' : { s_revcom.append("A"); break; }
		}
	}

	return s_revcom;
}

int GenerateSequences::CalculateEditDistance(string s1, string s2){
	
	int x,y;
	string s2_revcom;
	//unsigned int **d = new unsigned int*[SeqLen+1];
	const int len=SeqLen;

	/*d[0][0] = 0;
	
	for(x = 1; x <= len; ++x) d[x][0] = x;
	for(x = 1; x <= len; ++x) d[0][x] = x;*/
	
	for(x = 1; x <= len; ++x)
			for(y = 1; y <= len; ++y)
				d[x][y] = std::min( std::min(d[x - 1][y] + 1,d[x][y - 1] + 1),d[x - 1][y - 1] + (s1[x - 1] == s2[y - 1] ? 0 : 1) );

	return d[len][len];
}

string GenerateSequences::AppendAdaptors(string &PutBarcode){

	string s;
	s.clear();

	if(LeftAdaptor.length()!=0)
		s.append(LeftAdaptor);
	
	s.append(PutBarcode);  //IlluminaAHandle+PutBarcode
	
	if(RightAdaptor.length()!=0)
		s.append(RightAdaptor);

	return s;
}



void GenerateSequences::DetermineFilteringThresholds(void){

int i,lwzscore;
string s1,s2,s3;

for(i=0;i<SeqLen;i++)	s1.append("A");

for(i=0;i<SeqLen/2;i++)		s2.append("A");
/*for(i=SeqLen/2;i<SeqLen;i++)	s2.append("T");
for(i=0;i<SeqLen;i+=4){
	s3.append("A"); s3.append("T"); s3.append("C"); s3.append("G");
}
*/
lwzscore=lzw(s1); //LWZ score of the least complex sequence
LenDiffThreshold=SeqLen-lwzscore;

}

typedef struct {
  int gcmin,gcmax;
	GenerateSequences* father;
}params;

bool addToPool(string seq, vector<string>& buffer)
{
	//Acquire lock
	int poolSize = 0;
	if (pthread_mutex_trylock(&poolMutex) != 0) {//Mutex held by someone else
		if (buffer.size() > 50) {//Don't buffer more than 100 sequences
			pthread_mutex_lock(&poolMutex);
			if (die)
				return false;
			//Drop entire buffer to pool
			for (int i=0; i< (int)buffer.size(); ++i) {
				pool.push_back(buffer[i]);
			}
			poolSize = pool.size();
			pthread_mutex_unlock(&poolMutex);
			buffer.clear();
		} else
			buffer.push_back(seq);
	} else {
		if (die)
			return false;
		//Do stuff
		pool.push_back(seq);
		poolSize = pool.size();
		pthread_mutex_unlock(&poolMutex);
	}
	//Is the pool completely full?
	if (poolSize > 5000) {//Sleep around for some time... (:P)
		do {
#ifdef __linux__
			pthread_yield();
			usleep(1000);
#endif
		}while (pool.size() != 0);
	}
	return true;
}

void* generateRandomChecked(void* args)
{
	//Seed random
	srand(time(NULL));
	string seq,seq_rc,probe;
	int lzwscore,ed;
	double dG, dS,dH,Tm;
	bool passedrepeatcheck;
	params *p = (params*)args;
	int gcmin = p->gcmin;
	int gcmax = p->gcmax;
	int GCcontent = 0;
	HybridSSMin* hybMin=new HybridSSMin(SelfHybT,Hyb_Temperature);
	GenerateSequences* mother = p->father;
	vector<string> buffer;
	while (true){
	  int x = (int) (rand()%(gcmax-gcmin+1)); //SET THE GC CONTENT
	  GCcontent = gcmin + x; //SET THE GC CONTENT
	  seq = mother->randomseq_setGC(GCcontent); //Generate Random Seq

	  probe= GenerateSequences::AppendAdaptors(seq);
	  //Check self-hybridization
	  hybMin->computeGibsonFreeEnergy(dG,dH,probe.c_str());
	  dS=(dH-dG)/(273.15+ SelfHybT);
	  Tm=(dH/dS)-273.15;
	  if(Tm<=(SelfHybT+(0.1*SelfHybT))) {
	    //CHECK ILLUMINA HANDLE VS BARCODE HYB
	    bool failedHandleHyb = false;
	    for (int i=0; i < (int)AdaptorList.size(); ++i) {
	      hybMin->computeTwoProbeHybridization(dG,dH,seq.c_str(),AdaptorList[i].c_str());
	      dS=(dH-dG)/(273.15+ Hyb_Temperature);
	      Tm=dH/(dS+R*log(0.00001/4));
	      Tm-=273.15;
	      if(Tm>(Hyb_Temperature+(0.1*Hyb_Temperature))) {//Forget this one!
		failedHandleHyb = true;
		break;
	      }
	    }
	    if (!failedHandleHyb) {
	      //Check for repeats
	      passedrepeatcheck=mother->checkforruns(seq); // Check for repeats
	      if (passedrepeatcheck) {
		//Check for complexity
		lzwscore=lzw(seq);
		if(lzwscore<=LenDiffThreshold) {
		  //Check for edit distance to reverse complement
		  seq_rc = mother->RevComp(seq);
		  ed = mother->CalculateEditDistance(seq,seq_rc);
		  if(ed>=EditDistanceThreshold_Self){
		    //Check if it contains restriction sites from given list.
		    bool failedRestrictionSites = false;
		    for (int i=0; i<(int)RestrictionList.size();++i) {
		      if (seq.find(RestrictionList[i]) != string::npos) {//Found overlap with restriction site
			failedRestrictionSites = true;
			break;
		      }
		    }
		    if (!failedRestrictionSites) {//Passed everything
		      if (!addToPool(seq,buffer))//Will return false if job is over.
			break;
		    }
		  }
		}
	      }
	    }
	  }
	}
	delete hybMin;
	delete mother;
	return NULL;
}

void GenerateSequences::GenerateRandomSequence_SetGC(){
	
	string Tstr,passedseq;

	int gcmin,gcmax;
	gcmin=MinGC*SeqLen;	gcmax=MaxGC*SeqLen;

	std::stringstream temperature;
	temperature << SelfHybT;

	//Start threads!
	pthread_attr_t attr;
	/* Initialize and set thread detached attribute */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
	pthread_t someThread;
    //Create the world!
	for (int i=0;i<1;++i) {
	  	params* p = new params();
		p->gcmin = gcmin;
		p->gcmax = gcmax;
		p->father = new GenerateSequences();
		int rc = pthread_create(&someThread,&attr,&generateRandomChecked,p);
		if (rc != 0) {//Thread creation failed!
			fprintf(stderr,"Unable to create thread. Error code %d returned\n",rc);
		}
	}

}


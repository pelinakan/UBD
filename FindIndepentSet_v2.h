// NO BARCODE FILTERING, NODES ARE SELECTED AT RANDOM TO ADD TO THE UNIQUE SET
class Network{

public:
	Network();
	~Network();
	vector<string> PutBarcodes; //Keep the putative barcodes
	unsigned long int NofPutBarcodes;
	vector <string> CommonSet;
	vector< vector <double> > Features;
	bool RetrieveUniqueNodes(string);
	void PrintUniquePutBarcodes(const char*);
	static string RevComp(string);
	void PrintEditDistanceDistribution(char *fname);
private:
	int CalculateEditDistance(string, string, unsigned int **d);
	unsigned int* editDistDistribution;
	unsigned int*** d;
};

Network::Network(){
	editDistDistribution = new unsigned int[SeqLen+1];
	for (int i=0; i<=SeqLen; ++i)
		editDistDistribution[i] = 0;
	//Create matrices for all edit distance computations
	int maxThreads = omp_get_max_threads();
	d = new unsigned int**[maxThreads];
	for (int k=0; k < maxThreads; ++k ) {
		d[k] = new unsigned int*[SeqLen+1];
		for (int i=0;i<=SeqLen;++i) {
			d[k][i] = new unsigned int[SeqLen+1];
		}
		d[k][0][0] = 0;
		for(int x = 1; x <= SeqLen; ++x) d[k][x][0] = x;
		for(int x = 1; x <= SeqLen; ++x) d[k][0][x] = x;
	}
}

Network::~Network() {
	delete[] editDistDistribution;
	int maxThreads = omp_get_max_threads();
	for (int k=0; k < maxThreads; ++k ) {
		for (int i=0;i<=SeqLen;++i) {
			delete[] d[k][i];
		}
		delete[] d[k];
	}
	delete[] d;
}

string Network::RevComp(string s){	
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

int Network::CalculateEditDistance(string s1, string s2, unsigned int **d){
	
	int x,y;
	string s2_revcom;
	const int len=SeqLen;

	if (d == NULL) {
		fprintf(stdout,"Something went very wrong\n");
		return 0;
	}

	for(x = 1; x <= len; ++x)
			for(y = 1; y <= len; ++y)
				d[x][y] = std::min( std::min(d[x - 1][y] + 1,d[x][y - 1] + 1),d[x - 1][y - 1] + (s1[x - 1] == s2[y - 1] ? 0 : 1) );
	
	if(d[len][len] > EditDistanceThreshold){
	  s2_revcom=Network::RevComp(s2);

		for(x = 1; x <= len; ++x)
			for(y = 1; y <= len; ++y)
				d[x][y] = std::min( std::min(d[x - 1][y] + 1,d[x][y - 1] + 1),d[x - 1][y - 1] + (s1[x - 1] == s2_revcom[y - 1] ? 0 : 1) );
	}

	int ret = d[len][len];

	return ret;
}

bool Network::RetrieveUniqueNodes(string putbarcode){

long int k;
bool addtoset=true;
bool done=false;

unsigned int* localDistribution = new unsigned int[SeqLen+1];
for (int i=0; i<=SeqLen; ++i) localDistribution[i] = 0;

#pragma omp parallel for default(shared) schedule(dynamic,200)
 for(k=0;k<CommonSet.size();++k){ // Make sure if the last selected node is not connected to already present nodes in the common set
 
	 int id = omp_get_thread_num();
   if(!done){
     unsigned int editdist = CalculateEditDistance(putbarcode,CommonSet[k],d[id]); // Check reverse complement also
	 ++localDistribution[editdist];
     if(editdist<=EditDistanceThreshold){
#pragma omp critical
       {
	 addtoset=false;
	 done=true;
       }
     }
   }
}

if(addtoset) {
  CommonSet.push_back(putbarcode); //Add it to the set if it not connected to the previously added members
  //Also add the distribution parameters
  for (int i=0; i<SeqLen; ++i) {
	  editDistDistribution[i] += localDistribution[i];
	}
 }
 delete[] localDistribution;

 if(CommonSet.size()%500==0) {
	cout << CommonSet.size() << "       Unique Barcodes Selected" << endl;
 }


 return addtoset;

}

void Network::PrintUniquePutBarcodes(const char *fname){

	int i=0;
	ofstream outf(fname);
	
	for(i=0;i<CommonSet.size();i++){
		outf << CommonSet[i] << endl;
	}
	outf.close();

}

void Network::PrintEditDistanceDistribution(char *fname)
{
	// Print Edit Distance Distribution
	ofstream outf2(fname);

	for(int ed=0;ed<SeqLen;ed++){
		outf2 << ed << '\t' << editDistDistribution[ed] << endl;
	}
	outf2.close();
}

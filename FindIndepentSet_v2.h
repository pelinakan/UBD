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
	int CalculateEditDistance(string, string);
	unsigned int* editDistDistribution;
};

Network::Network(){
	editDistDistribution = new unsigned int[SeqLen+1];
	for (int i=0; i<=SeqLen; ++i)
		editDistDistribution[i] = 0;
}

Network::~Network() {
	delete[] editDistDistribution;
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

int Network::CalculateEditDistance(string s1, string s2){
	
	int x,y;
	string s2_revcom;
	unsigned int **d = new unsigned int*[SeqLen+1];
	for (int i=0;i<=SeqLen;++i) {
		d[i] = new unsigned int[SeqLen+1];
	}
	d[0][0] = 0;
	for(int x = 1; x <= SeqLen; ++x) d[x][0] = x;
	for(int x = 1; x <= SeqLen; ++x) d[0][x] = x;
	const int len=SeqLen;

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
	for (int i=0;i<=SeqLen;++i) {
		delete[] d[i];
	}
	delete[] d;

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
   //#pragma omp critical
   if(!done){
     unsigned int editdist = CalculateEditDistance(putbarcode,CommonSet[k]); // Check reverse complement also
     if (editdist > SeqLen)
       fprintf(stdout,"oups\n");
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

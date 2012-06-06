class DegreeDistribution{
	friend class Network;
	unsigned long long int *EDDistribution; // EDDistribution: EDDistribution[E]: How many pairs is with an edit distance E
	unsigned long long int *EditDist; // EditDist[index] : The index of this array code for the pair and its value is the edit distance.
	unsigned int *dd; //Degree Distribution
public:
	unsigned long long int NofIndexes;
	unsigned long int NofBarcodes;
	DegreeDistribution();
	vector<bool> Connected; // Connected[x] : If pair x is connected or not.
	void initialisevars();
	int CalculateEditDistance(string, string);
	void GenerateNetwork(vector <string>&);
	unsigned long long FindMin(void);
	unsigned long long FindMax(void);
	void DegreeDist(char*, char * /*,vector <string>&*/);
	string RevComp(string);

};

DegreeDistribution::DegreeDistribution(){
}

void DegreeDistribution::initialisevars(){

	unsigned long long int x;
	x=(NofBarcodes*(NofBarcodes+1))/2;
	NofIndexes=(NofBarcodes*(NofBarcodes-1)/2);

	EditDist=(unsigned long long int*)  malloc((x)*sizeof(unsigned long long int));
	EDDistribution=(unsigned long long int*) malloc(SeqLen*sizeof(unsigned long long int));
	dd=(unsigned int*) malloc(NofBarcodes*sizeof(unsigned int));

}
string DegreeDistribution::RevComp(std::string s){	
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
int DegreeDistribution::CalculateEditDistance(string s1, string s2){
	
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
		s2_revcom=RevComp(s2);

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

void DegreeDistribution::GenerateNetwork(vector <string> &Barcodes){
//	ofstream tempf("temp_debug.txt");

NofBarcodes=Barcodes.size();
initialisevars();

	long int i,j;
	for(i=0;i<SeqLen;i++)
		EDDistribution[i]=0;
	
	cout << "Within GenerateNetwork Func" << endl;

	for(i=0;i<NofBarcodes-1;++i){
#pragma omp parallel for
		for(j=i+1;j<NofBarcodes;++j){
			unsigned long long int index=0;
			index=((((NofBarcodes*(NofBarcodes-1)/2)- ((NofBarcodes-i)*((NofBarcodes-i-1))/2))))+ (j-i-1); // index of the pair
			EditDist[index]=CalculateEditDistance(Barcodes[i],Barcodes[j]);
			++EDDistribution[EditDist[index]];
		}
		if(i%50==0)
			cout << i << endl;
	}
	NofIndexes=((NofBarcodes)*(NofBarcodes-1))/2;
	
	cout << "Network Connected" << endl;

}
unsigned long long DegreeDistribution::FindMin(void){

	unsigned long long i, m=EditDist[0],mineditdist=0;

	for(i=1;i<NofIndexes;i++){
		if(m>EditDist[i]){
			m=EditDist[i];
			mineditdist=EditDist[i];
		}
	}
	return mineditdist;
}
unsigned long long DegreeDistribution::FindMax(void){

	unsigned long long i;
	unsigned long long m,maxeditdist=0;

	m=EditDist[0];
	for(i=1;i<NofIndexes;i++){
		if(m<EditDist[i]){
			m=EditDist[i];
			maxeditdist=EditDist[i];
		}
		if(m<0)
			break;
	}
	return maxeditdist;
}
void DegreeDistribution::DegreeDist(char *fname, char *fname2 /*,vector <string> &IndependentSet*/){

unsigned long long i;
unsigned long int edi=0,j,nodeindex=0,nodedegree;
unsigned long long MinED, MaxED,ed,lastnodedd=0;

vector< vector < unsigned long int> > dd; //Degree Distribution

cout << "Entered Degree Distribution Function and created the vector" << endl;

MinED=FindMin();
cout << "MinED   " << MinED << endl;
MaxED=FindMax();
cout << "MaxED   " << MaxED << endl;

for(ed=MinED;ed<=MaxED;ed++){ // For each edit distance
	lastnodedd=0;
	dd.push_back(vector<unsigned long int>());  // Create one vector array
	dd[edi].push_back(0); // The first element
	nodeindex=0; // Index of the node need to be kept outside the for loop
	for(i=0;i<NofIndexes;){ // look through all pairs 
		nodedegree=0;       // Initialise the degree of the node
		for(j=0;j<NofBarcodes-nodeindex-1;j++){ // For each node, go through the other nodes (N-nodeindex-1)
			if(EditDist[i]<=ed) // If it is connected for that particular ed, i keeps the pair index
				++nodedegree; // Increase its degree by one
			i++; // 
		}
		//for the last node degree
		if(EditDist[i-1]<=ed)
			lastnodedd++;
		++nodeindex;
		if(nodedegree==0)
			++dd[edi][nodedegree];
		else{
			if(dd[edi].size()<=nodedegree){ // if the size of the array is smaller than required
				dd[edi].resize(nodedegree+1); // resize it
				++dd[edi][nodedegree]; // increase the number of nodes with degree "nodedegree" by 1
			}
			else
				++dd[edi][nodedegree]; // increase the number of nodes with degree "nodedegree" by 1
		}
	}
	if(dd[edi].size()>100)
		break;
	if(dd[edi].size()<=lastnodedd){
		dd[edi].resize(lastnodedd+1);
			++dd[edi][lastnodedd];
	}
	else
		++dd[edi][lastnodedd];

	++edi;
}

cout << "Degree Distribution Vector Array Filled" << endl;

ofstream outf(fname); // Degree Distribution

edi=0;
for(ed=MinED;ed<=MaxED;ed++){
	if(dd[edi].size()>100)
		break;
	outf << ed << '\t';
	for(i=0;i<dd[edi].size();i++)
		outf << dd[edi][i] << '\t';
	outf << endl;
	edi++;
}
outf.close();

// Print Edit Distance Distribution
ofstream outf2(fname2);

	for(ed=MinED;ed<=MaxED;ed++){
		outf2 << ed << '\t' << EDDistribution[ed] << endl;
	}
outf2.close();

}



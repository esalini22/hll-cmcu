#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <errno.h>

using namespace std;

void readtxtbin(string filename, int n) {
    ifstream *infile;
    uint64_t *buffer = (uint64_t *)malloc(n*sizeof(uint64_t));
    string out=filename+".bin";
    FILE *fout = fopen(out.c_str(), "wb");

    infile = new ifstream(filename.c_str());
    if (!infile) {
        printf("Cannot open file %s ", filename.c_str());
	    return;
    }
    string line; // current line
    uint64_t i=0;
    while (getline(*infile, line)) {
    	stringstream ss(line);
            string word;
    	ss>>word;
    	uint64_t ip = stoll(word);
    	buffer[i] = ip;
    	i++;
    }
    delete infile;
    fwrite(buffer,sizeof(uint64_t),n,fout);
    fclose(fout);
    free(buffer);
}

void readbin(string filename, int n) {
    FILE *fr = fopen(filename.c_str(), "rb");
    uint64_t *buffer = (uint64_t *)malloc(n*sizeof(uint64_t));
    fread(buffer,sizeof(uint64_t),n,fr);
    for(int i=0; i<n; i++){
	    cout<<buffer[i]<<endl;
    }
    free(buffer);
}

    
int main(int argc, char *argv[]) {

    if(argc !=3){
	    cout<<" faltan argumentos\n";
	return 1; 
    }

    int n = atoi(argv[2]); //cantidad de lineas
    string bname = string(argv[1]);

    //cuenta lineas
    /*FILE *fp;
    fp = fopen(argv[1], "rb");
    while(!feof(fp))
        if(fgetc(fp)=='\n') ++n;
    fclose(fp);*/

    readtxtbin(bname,n);
    //readbin(bname+".bin",n);
    
}

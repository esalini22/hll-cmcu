#include <stdio.h>
#include <string.h>
#include <string>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <thread>
#include <vector>
#include <algorithm>
#include <vector>

#include "pq_array.hpp"
#include "count_min_sketch.hpp"
#include "murmurhash.hpp"

#include <unordered_map>
#include <map>
#include <memory>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>

#include "safe_ptr.h"
#include "hll_fun.h"
#include <omp.h>

using namespace std;

vector<uint64_t> hll_sketch;
vector<uint64_t> bit_mask;
vector<CountMinSketch*> cu;
vector<pq_array*> pq_arr;
//unordered_map<uint64_t, uint32_t> map_hh;
vector<unordered_map<uint64_t, uint32_t>> map_hh;
int M = 0;
int fpint = -1;

vector<sf::contention_free_shared_mutex<>> s_m;

int reader(string fname, int chunk, int id, int p)
{
    uint64_t *bufferA = (uint64_t *)malloc(chunk*sizeof(uint64_t));
    int readSz = pread(fpint, bufferA, chunk*sizeof(uint64_t), id*chunk*sizeof(uint64_t));
	int M_i=0;

    for(int i=0; i<chunk; i++){
    	M_i++;
	    uint64_t ip = bufferA[i];
		//printf("ip: %lu i: %d chunk: %d\n",ip,i,chunk);
	    uint64_t hh = murmur64(ip);
		uint32_t hash = hh&((1UL<<32)-1);
		uint32_t v1 = hash&((1<<p)-1); //indice

		s_m[v1].lock();

		//insertamos en sketch hll
		//restamos 32 ya que el primer resultado es long long, es decir 64 bits
		uint8_t w = __builtin_clzll(hash >> p) - static_cast<uint8_t>(p) + 1 - 32;
		unsigned int indice=v1/12; //indice de celda
		unsigned int reg=v1%12; //12 registros por celda, encontramos el registros correspondiente
		uint64_t temp=(hll_sketch[indice]>>(5*reg))&0x1F; //obtenemos el registro guardado en el sketch
		if(w > temp) hll_sketch[indice]=(hll_sketch[indice]&(bit_mask[reg]))|((uint64_t)w<<(5*reg));

		//insertamos en sketch cm-cu y pq
		cu[v1]->updatecu(ip,1);
		int est = cu[v1]->estimate(ip);

		//heavy hitters
		pq_arr[v1]->add(ip,est);

		//insertamos en tabla hash para entropia real
		if (map_hh[v1].find(ip) != map_hh[v1].end())
			map_hh[v1][ip] ++;
		else map_hh[v1][ip] = 1;

		s_m[v1].unlock();
    }
    
    free(bufferA);

    return M_i;
}

//obtiene archivo de la linea de argumentos
string getTraces(char** argv, int argc){
	for(int i=1;i<argc;++i){
		if(!strcmp(argv[i],"-k") || !strcmp(argv[i],"-p") || !strcmp(argv[i],"-d") || !strcmp(argv[i],"-w") || !strcmp(argv[i],"-t") || !strcmp(argv[i],"-o") || !strcmp(argv[i],"-d") || !strcmp(argv[i],"-r")) ++i;
		else if(strcmp(argv[i],"-s")) return argv[i];
	}
	return "";
}

//formato: ./hll -opcion valor trazas.txt
//o bien ./hll trazas.txt -opcion valor
//no detecta caso en que se introduza opcion o valor invalido
int main(int argc, char *argv[]){
	if(argc<2) {
		printf("No hay suficientes argumentos\n");
		exit(1);
	}
	unsigned char p=7;
	int d_cmcu=6,w_cmcu=64,nth=1;
	int lines=0;
	
	char** option;
	char** end=argv+argc;
	option=std::find((char**)argv,end,(const std::string&)"-p");
	if(option!=end){
		char val=atoi(*(option+1));
		if(val<16 && val>5) p=val;
	}

	option=std::find((char**)argv,end,(const std::string&)"-d");
	if(option!=end){
		int val=atoi(*(option+1));
		d_cmcu=val;
	}
	option=std::find((char**)argv,end,(const std::string&)"-w");
	if(option!=end){
		int val=atoi(*(option+1));
		w_cmcu=val;
	}
	option=std::find((char**)argv,end,(const std::string&)"-t");
	if(option!=end){
		int val=atoi(*(option+1));
		nth=val;
	}
	option=std::find((char**)argv,end,(const std::string&)"-l");
	if(option!=end){
		int val=atoi(*(option+1));
		lines=val;
	}

	string traces = getTraces(argv,argc);

	if(traces=="") {
		printf("No se indico trazas\n");
		exit(1);
	}

	printf("p:%d d:%d w:%d\n",p,d_cmcu,w_cmcu);

	int N = 1<<p;
	cu.resize(N);
	pq_arr.resize(N);
	int elems_per_queue = log2(8192/(8*N));
	for(int i=0; i<N; i++){
		cu[i] = new CountMinSketch(w_cmcu,d_cmcu);
		pq_arr[i] = new pq_array(8, elems_per_queue, 25, 1);
	}

	//inicializamos el sketch hll en ceros
	for(unsigned int i=0;i<N/12;++i) //cada celda tendra 12 buckets del sketch array
		hll_sketch.emplace_back(0);
	//en caso de que sobre una celda
	if(N%12) hll_sketch.emplace_back(0);

	//inicializamos mascaras de bits
	for(unsigned char i=0;i<12;++i)
		bit_mask.emplace_back(~((uint64_t)0x1F<<(5*i)));

	//habra 1 mutex por cada celda del sketch
	vector<sf::contention_free_shared_mutex<>> list1(N);
	s_m.swap(list1);

	if(!lines){
		FILE *fp;
		fp = fopen(traces.c_str(), "rb");
		while(!feof(fp))
			if(fgetc(fp)=='\n') ++lines;
		fclose(fp);
	}

	int chunk = lines/nth;
    string fname = traces+".bin";
    fpint = open(fname.c_str(), O_RDWR | O_CREAT, S_IREAD | S_IWRITE | S_IRGRP | S_IROTH);

    int M_i[nth];
	map_hh.resize(N);

    omp_set_num_threads(nth);
	#pragma omp parallel
	{
		#pragma omp single
		{
			for(int i=0;i<nth;++i){
				#pragma omp task
				M_i[i]=reader(fname,chunk,i,p);
			}
		}
	}
	close(fpint);

	for(int i=0;i<nth;++i) M+=M_i[i];

	unordered_map<uint64_t, uint32_t> map_hh_real;
	for(int i=0;i<N;++i)
		for(unordered_map<uint64_t, uint32_t>::iterator it=map_hh[i].begin();it!=map_hh[i].end();++it){
			if (map_hh_real.find(it->first) != map_hh_real.end())
				map_hh_real[it->first]+= it->second;
			else map_hh_real[it->first] = it->second;
		}

	//entropia estimada
	vector<uint32_t> counters_est;
	for(uint32_t i=0;i<N;++i){
		vector<uint32_t> temp_vec = pq_arr[i]->get_data();
		counters_est.insert(counters_est.begin(),temp_vec.begin(),temp_vec.end());
	}
	sort(counters_est.begin(),counters_est.end(),std::greater<uint32_t>());
	unsigned K=8192,L=0;
	double est_entropy=0;
	for(uint32_t i=0;i<K;++i){ //solo top-k
		L+=counters_est[i];
		est_entropy+=(counters_est[i]/(double)M)*log2(counters_est[i]/(double)M);
	}

	size_t card=est_card(p,hll_sketch);

	double est_left=est_entropy;
	double est_right=((M-L)/(double)M)*log2((M-L)/(double)(M*(card-K)));
	est_entropy+=est_right;
	est_entropy=-est_entropy/(double)log2(card);

	vector<uint32_t> counters_real;
	for(unordered_map<uint64_t,uint32_t>::iterator it=map_hh_real.begin();it!=map_hh_real.end();++it)
		counters_real.emplace_back(it->second);
	sort (counters_real.begin(), counters_real.end(), greater<int>());
	double left_entropy=0,right_entropy=0;
	for (size_t i = 0; i < counters_real.size (); ++i){
		if(i<8192) left_entropy -= (counters_real[i]/(double)M)*(log2(counters_real[i]/(double)M));
		else right_entropy -= (counters_real[i]/(double)M)*(log2(counters_real[i]/(double)M));
	}
	double true_entropy=left_entropy+right_entropy;
	true_entropy=true_entropy/(double)log2(counters_real.size());
	printf("True entropy: %lf\n",true_entropy);
	printf("Estimated entropy: %lf\n",est_entropy);
	printf("ER: %lf\n",abs(est_entropy-true_entropy)/true_entropy);
	return 0;
}

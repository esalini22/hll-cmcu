#include "HyperCount.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <omp.h>
#include <zlib.h>
#include <algorithm>
#include "wyhash32.h"
//#include "murmurhash32.hpp"
//#include "wyhash.h"

#include <map>
#include <set>
#include <unordered_map>

using namespace std;
#define seed 5 //seed para hash
#define lim 4294967296 //2^32

inline float HyperCount::hsum_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 maxs = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, maxs); // high half -> low half
    maxs        = _mm_add_ss(maxs, shuf);
    return        _mm_cvtss_f32(maxs);
}

inline float HyperCount::hsum_avx(__m256 v) {
    __m128 lo = _mm256_castps256_ps128(v);   // low 128
    __m128 hi = _mm256_extractf128_ps(v, 1); // high 128
           lo = _mm_add_ps(lo, hi);          // max the low 128
    return hsum_sse3(lo);                    // and inline the sse3 version
}

HyperCount::HyperCount(unsigned char n1, unsigned char n2, unsigned char d, unsigned char w){
	//omp_set_num_threads(numThreads);
	p=n1;
	b=n2;
	bits_v2=(1<<b)-1; //(2^b)-1 ,causa stack smashing con 64 bits
	N=1<<p; //2^V1
	for(unsigned char i=0;i<16;++i)
		bit_mask.emplace_back(~((ullint)0xF<<(4*i)));
	sketch_size=to_string(p);
	min_val=1<<(b-15);

	a_m=(0.7213/(1+(1.079/N)))*N*N;
	ciclos_red=(b+2)/8+(((b+2)%8)>0);
	M=0;
	//samples.resize(N);

	cmcu.resize(N);
	for(unsigned int i=0;i<N;++i){
		CountMinCU* cu = new CountMinCU(d,w);
		cmcu[i] = cu;
	}
}
HyperCount::~HyperCount(){}

void HyperCount::addSketch(string genome){
	//inicializa sketch en 0s, y añade su nombre a un vector
	for(unsigned int j=0;j<N/16;++j) //cada celda tendra 12 buckets del sketch array
		sketch.emplace_back(0);
	//if((long)N%16)
		//sketch.emplace_back(0);
	name=genome;
	printf("%s\n",name.c_str());
}

void HyperCount::insert(ullint kmer){
	M++;
	//calculamos el hash usando wyhash de 32 bits (hashing muy rapido)
	const ullint *key = &kmer;
	//printf("%llu %llu\n",kmer,(ullint)key);
	//uint64_t hash;
	//wyhash((const uint64_t*)key,8,5,&hash);

	uint32_t hash=wyhash32(key,8,0x5);
	//uint32_t hash=murmurhash((const uint64_t*)key,5);

	//cout<<kmer<<" "<<key<<" "<<hash<<endl;


	unsigned int v1=(hash>>b); //se desplaza el hash b bits a la derecha
	unsigned int v2=hash&bits_v2; // AND 2^b -1 
	unsigned char w;
	if(v2<min_val) w=15;
	else w=__builtin_clz(v2)+1-p; /*cuenta cantidad de bits 0 mas significativo
	lo cual sirve para encontrar posicion de 1 mas a la izquierda
	se resta p, ya que al trabajar con int, se tiene 32 bits,
	y siempre habran al menos p ceros a la izquierda*/

	unsigned int indice=v1/16; //indice de celda
	unsigned char bucket=v1%16; //12 buckets por celda
	ullint temp=(sketch[indice]>>(4*bucket))&0xF;

	cmcu[v1]->insertCMinCU(v2);

	if(w > temp) sketch[indice]=(sketch[indice]&(bit_mask[bucket]))|(w<<(4*bucket));
}

void HyperCount::loadSketch(char *infilename){ //¿que pasa si tamaño del sketch es distinto al de argumento?
	string genome=infilename; //se debe determinar nombre del genoma a partir del archivo .hll
	printf("infilename: %s\n",infilename);
	genome.erase(genome.find(".w."),string::npos);
	printf("genome: %s\n",genome.c_str());

	//debemos descomprimir el archivo y ponerlo en un vector
	vector<ullint> v;
	gzFile infile = gzopen(infilename, "rb");
	//if(!infile) return -1;
	char buffer[128];
	int num_read = 0;
	char temp[19];
	char* end;
	int cont=0;
	while((num_read = gzread(infile,buffer, sizeof(buffer))) > 0){
		for(int i=0;i<num_read;++i){
			temp[cont%19]=buffer[i];
			++cont;
			if(cont%19==0){
				cont=0;
				v.emplace_back(strtoull(temp,&end,10));
			}
		}
	}
	if(cont>0) v.emplace_back(strtoull(temp,&end,10));
	v.emplace_back(0);
	gzclose(infile);

	sketch=v;
	name=genome;
	printf("done\n");
}

void HyperCount::saveSketch(){ //usar zlib para comprimir
	printf("save sketch\n");
	string temp_name=name;
	temp_name+=".spacing."+sketch_size+".hll";
	char* outfilename=(char*)temp_name.c_str();
	printf("%s\n",outfilename);
	gzFile outfile = gzopen(outfilename, "wb");
	//if(s.empty() || !outfile) return -1;
	char inbuffer[128];
	char cont=0;
	for(vector<ullint>::iterator it=sketch.begin();it!=sketch.end();++it){
		string temp=to_string(*it)+'\n';
		int len=temp.length();
		for(char i=0;i<len;++i){
			inbuffer[cont]=temp[i];
			++cont;
			if(cont==127){
				gzwrite(outfile, inbuffer, cont);
				cont=0;
			}
		}
	}
	if(cont>0) gzwrite(outfile, inbuffer, cont);
	gzclose(outfile);
}

float HyperCount::estCard(vector<ullint> v){ //cardinalidad individual
	__m256i vec3; //vector de mascaras de bits
	ullint bits_and[4]={0xF,0xF,0xF,0xF};
	vec3=_mm256_loadu_si256((const __m256i*)&bits_and[0]);

	vector<float> w(32,0.0);
	vector<ullint>::iterator it1=v.begin();
	vector<ullint>::iterator fin=v.end();
	//contamos las repeticiones de w en cada sketch
	while(it1!=fin){ 
		ullint i1=*it1;
		ullint it_array[16]={i1,(i1>>4),(i1>>8),(i1>>12),(i1>>16),(i1>>20),(i1>>24),(i1>>28),(i1>>32),(i1>>36),(i1>>40),(i1>>44),(i1>>48),(i1>>52),(i1>>56),(i1>>60)};
		__m256i vec4[4];
		vec4[0]=_mm256_loadu_si256((const __m256i*)&it_array[0]);
		vec4[1]=_mm256_loadu_si256((const __m256i*)&it_array[4]);
		vec4[2]=_mm256_loadu_si256((const __m256i*)&it_array[8]);
		vec4[3]=_mm256_loadu_si256((const __m256i*)&it_array[12]);
		vec4[0]=_mm256_and_si256(vec3,vec4[0]);
		vec4[1]=_mm256_and_si256(vec3,vec4[1]);
		vec4[2]=_mm256_and_si256(vec3,vec4[2]);
		vec4[3]=_mm256_and_si256(vec3,vec4[3]);
		__attribute__ ((aligned (32))) ullint out[4];
		for(int c=0;c<4;++c){
			_mm256_store_si256((__m256i *)&out[0],vec4[c]);
			w[out[0]]++;
			w[out[1]]++;
			w[out[2]]++;
			w[out[3]]++;
		}
		++it1;
	}
	/*while(it1!=fin){ 
		ullint i1=*it1;
		for(char i=0;i<16;++i){ //12 registros por celda
			ullint temp1=i1&0xF;
			w[temp1]++;
			i1=i1>>4;
		}
		++it1;
	}*/
	float card=0.0;
	float w2[32];
	for(int i=0;i<32;++i) w2[i]=1.0;
	int respow=1;
	for(int i=0;i<b+2;++i){
		w2[i]=(float)respow;
		respow=respow<<1;
	}
	__m256 vec,vec2;
	for(int i=0;i<ciclos_red;++i){
		vec=_mm256_loadu_ps((const float *)&w[i*8]);
		vec2=_mm256_loadu_ps((const float *)&w2[i*8]);
		vec=_mm256_div_ps(vec,vec2);
		card+=hsum_avx(vec);
	}

	//media armonica
	card=(float)a_m/card;
	int ceros = w[0];
	if(ceros && card<=5*N/2) //C_HLL, ln cuando hay muchos ceros;
		card=N*log(N/ceros);
	else if(card>lim/30)
		card=-lim*log(1-(card/lim));
	//printf("estimacion cardinalidad %s: %f\n",name.c_str(),card);
	return card;
}

vector<ullint> HyperCount::merge(HyperCount *hll){ //hace la union y devuelve un nuevo sketch
	vector<ullint> sketch2=hll->getSketch();
	vector<ullint> ret=max(sketch,sketch2);
	ret.resize(max(sketch.size(),sketch2.size()));
	vector<ullint>::iterator it2=ret==sketch?sketch2.begin():sketch.begin();
	
	/*for(unsigned int j=0;j<N;++j){
		if(){
		}
		else if(){
		}
		else{
		}
	}*/

	return ret;
}

const vector<ullint> &HyperCount::getSketch(){
	return sketch;
}
string HyperCount::getName(){
	return name;
}

void HyperCount::print(){
	
}

float HyperCount::entropy(){	
	vector<vector<uint32_t>> counters;
	for(uint32_t i=0;i<N;++i){
		//cmcu[i]->print();
		//printf("----------------\n");
		counters.push_back(cmcu[i]->topK());
	}
	unsigned L=0,K=0;
	float ent=0;
	for(uint32_t i=0;i<N;++i){
		for(vector<uint32_t>::iterator it=counters[i].begin();it!=counters[i].end();++it){
			L+=(float)*it;
			K++;
			ent+=((float)*it/(float)M)*log2((float)*it/(float)M);
		}
	}
	ent+=((float)(M-L)/(float)M)*log2((float)(M-L)/(float)(M*(estCard(sketch)-K)));
	return -(float)ent/(float)log2(estCard(sketch));
}
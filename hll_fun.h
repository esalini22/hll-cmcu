#ifndef HLL_FUN_H
#define HLL_FUN_H

#include <immintrin.h>
#include <vector>

#define lim 4294967296 //2^32

using namespace std;

inline float hsum_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 maxs = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, maxs); // high half -> low half
    maxs        = _mm_add_ss(maxs, shuf);
    return        _mm_cvtss_f32(maxs);
}

inline float hsum_avx(__m256 v) {
    __m128 lo = _mm256_castps256_ps128(v);   // low 128
    __m128 hi = _mm256_extractf128_ps(v, 1); // high 128
           lo = _mm_add_ps(lo, hi);          // max the low 128
    return hsum_sse3(lo);                    // and inline the sse3 version
}

/*float est_card(int p, vector<uint64_t> hll_sketch){
	int N=1<<p;
	int b=32-p;
	int a_m=(0.7213/(1+(1.079/N)))*N*N;
	if(N==64) a_m=0.709*N*N;
	vector<float> w(32,0.0);

	vector<uint64_t>::iterator it1=hll_sketch.begin();
	vector<uint64_t>::iterator fin=hll_sketch.end();
	while(it1!=fin-1){ 
		for(int i=0;i<12;++i){
			uint64_t i1=*it1;
			w[i1&0x1F]++;
			i1=i1>>5;
		}
		++it1;
	}
	unsigned char tam=(long)N%12;
	uint64_t i1=*it1;
	for(char i=0;i<tam;++i){ //12 registros por celda
		w[i1&0x1F]++;
		i1=i1>>5;
	}
	//eliminamos las repeticiones por celda extra vacía creada
	//w[0]-=12;

	float card=0;
	for(unsigned char i=0;i<b+2;++i)
		if(w[i]) card+=(float)w[i]/(float)(1<<i);

	int ceros=w[0];
	card=(float)a_m/card;
	if(ceros && card<=5*N/2) //C_HLL, ln cuando hay muchos ceros;
		card=N*log(N/ceros);
	else if(card>lim/30)
		card=-lim*log(1-(card/lim));
	printf("estimacion cardinalidad: %f\n",card);
	return card;
}*/

float est_card(int p, vector<uint64_t> hll_sketch){
	int N=1<<p;
	int b = 32-p;
	int ciclos_red=(b+2)/8+(((b+2)%8)>0);
	int a_m=(0.7213/(1+(1.079/N)))*N*N;
	if(N==64) a_m=0.709*N*N;

	__m256i vec3; //vector de mascaras de bits
	uint64_t bits_and[4]={0x1F,0x1F,0x1F,0x1F};
	vec3=_mm256_loadu_si256((const __m256i*)&bits_and[0]);

	vector<float> w(32,0.0);
	vector<uint64_t>::iterator it1=hll_sketch.begin();
	vector<uint64_t>::iterator fin=hll_sketch.end();
	//contamos las repeticiones de w en cada sketch
	while(it1!=fin-1){ 
		uint64_t i1=*it1;
		uint64_t it_array[12]={i1,(i1>>5),(i1>>10),(i1>>15),(i1>>20),(i1>>25),(i1>>30),(i1>>35),(i1>>40),(i1>>45),(i1>>50),(i1>>55)};
		__m256i vec4[3];
		vec4[0]=_mm256_loadu_si256((const __m256i*)&it_array[0]);
		vec4[1]=_mm256_loadu_si256((const __m256i*)&it_array[4]);
		vec4[2]=_mm256_loadu_si256((const __m256i*)&it_array[8]);
		vec4[0]=_mm256_and_si256(vec3,vec4[0]);
		vec4[1]=_mm256_and_si256(vec3,vec4[1]);
		vec4[2]=_mm256_and_si256(vec3,vec4[2]);
		__attribute__ ((aligned (32))) uint64_t out[4];
		for(int c=0;c<3;++c){
			_mm256_store_si256((__m256i *)&out[0],vec4[c]);
			w[out[0]]++;
			w[out[1]]++;
			w[out[2]]++;
			w[out[3]]++;
		}
		++it1;
	}
	unsigned char tam=(long)N%12;
	uint64_t i1=*it1;
	for(char i=0;i<tam;++i){ //12 registros por celda
		uint64_t temp1=i1&0x1F;
		w[temp1]++;
		i1=i1>>5;
	}
	//eliminamos las repeticiones por celda extra vacía creada
	//w[0]-=12;

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
	printf("estimacion cardinalidad: %f\n",card);
	return card;
}

#endif
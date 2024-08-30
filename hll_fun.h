#ifndef HLL_FUN_H
#define HLL_FUN_H

#include <vector>

#define lim 4294967296 //2^32

using namespace std;

float est_card(int p, vector<uint64_t> hll_sketch){
	float N=1<<p;
	int b=32-p;
	float a_m=(0.7213/(1+(1.079/N)));
	if(N==64) a_m=0.709;
	vector<float> w(b+2,0.0);

	vector<uint64_t>::iterator it1=hll_sketch.begin();
	vector<uint64_t>::iterator fin=hll_sketch.end();
	while(it1!=fin-1){ 
		uint64_t i1=*it1;
		for(int i=0;i<12;++i){
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

	float card=0.0;
	for(unsigned char i=0;i<b+2;++i)
		if(w[i]) card+=(float)w[i]/(float)(1<<i);

	int ceros=w[0];
	card=a_m*N*N/card;
	if(ceros && card<=5*N/2) //C_HLL, ln cuando hay muchos ceros;
		card=N*log(N/ceros);
	else if(card>lim/30)
		card=-lim*log(1-(card/lim));
	printf("card: %f\n",card);
	return card;
}

#endif

#include "CountMinCU.h"
#include "wyhash32.h"
#include <bits/stdc++.h>
#include <algorithm>
using namespace std;
CountMinCU::CountMinCU(int d, int w){
    filas=d;
    columnas=w;
    //this->C=new int* [filas]; //se crea la matriz de contadores
    C.resize(filas);
    for(int i=0;i<filas;i++) //se inicializan los contadores en 0
        C[i].resize(columnas);
    pq = new pq_array(8, 4, 25, 1);
}
CountMinCU::~CountMinCU(){}

void CountMinCU::insertCMinCU(ullint n){
    const ullint *key = &n;
    for(uint32_t i=0;i<filas;i++){ //se hace hashing a celdas para aumentar contadores
        uint32_t hash=wyhash32(key,8,i);
        if((C[i][hash%columnas])==estimarFreq(n)) //incrementa solo los contadores de los que poseen freq est
            C[i][hash%columnas]++;
    }
    pq->add(n,estimarFreq(n));
}
ullint CountMinCU::estimarFreq(ullint n){
    const ullint *key = &n;
    Freq_est=numeric_limits<int>::max();
    for(uint32_t i=0;i<filas;i++){
        uint32_t hash=wyhash32(key,8,i);
        Freq_est = min(Freq_est, C[i][hash%columnas]); //encuentra el contador mas pequeño
    }
    return Freq_est;
}

void CountMinCU::print(){
    cout<<"print"<<endl;
    for(uint32_t i=0;i<columnas;++i){
        for(uint32_t j=0;j<filas;++j)
            cout<<C[j][i]<<" ";
    }
    cout<<endl;
}

vector<uint32_t> CountMinCU::topK(){
    return pq->get_data();
}
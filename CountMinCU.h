#ifndef COUNTMINCU_H
#define COUNTMINCU_H
#include <iostream>
#include <vector>
#include <limits>
#include <bits/stdc++.h>
#include "pq_array.hpp"
using namespace std;
typedef unsigned long long int ullint;
class CountMinCU{
    private:
        int filas,columnas;
        ullint Freq_est;
        vector<vector<ullint>> C; //matriz de contadores
        pq_array *pq;

    public:
        CountMinCU(int d,int w);
        ~CountMinCU();
        void insertCMinCU(ullint n);
        ullint estimarFreq(ullint n);
        vector<uint32_t> topK();
        void print();
};
#endif

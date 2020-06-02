// #pragma once
#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include "Shape.h"
#include <vector>

class MonteCarlo
{
protected:
    Shape* shape;
    int step = 1;// 时间步长
    int time_max = 5;
    int center[3] = {10,10,10};
    std::vector<double> temp_1d;

public:
    MonteCarlo(Shape*);
    void setProbe(int x, int y, int z, double T);
    double iterate(); // 迭代一次，计算一个step产生的热传导，返回温度更新值
    void run(int verbose=0); // 不停止地直接迭代算法，直至到达稳态或到达最大迭代次数

    void output1d();
    void output3d();
};


#endif







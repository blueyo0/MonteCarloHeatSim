#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "Shape.h"
#include <vector>

enum SimMode {
	ONE_DIM, RANDOM_WALK
};


class MonteCarlo
{
protected:
	SimMode mode = SimMode::RANDOM_WALK;
	//SimMode mode = SimMode::ONE_DIM;
    Shape* shape;
    int step = 1;// 时间步长
    int time_max = 2000;
	int curr_iterate_time = 0;
	double default_value = 37.0;
	int monteCarloNum = 100;
    int center[3];

    std::vector<double> temp_1d;
	double ***temp_3d;

public:
    MonteCarlo(Shape*);
	void setIteration(int num);
    void setProbe(int x, int y, int z, double T);
	void reset();
	void setDefaultValue(double v);
    double iterate(); // 迭代一次，计算一个step产生的热传导，返回温度更新值
    void run(int verbose=0); // 不停止地直接迭代算法，直至到达稳态或到达最大迭代次数

	void runWithOneDim(int verbose = 0);
	void runWithRandomWalk(int verbose = 0);

	double computeTempWithRandomWalk(int x, int y, int z, int n);//RandomWalk模式的单个位置计算
	std::vector<double> computeProbs(int x, int y, int z);

    void output1d();
    void output3d();
};

#endif

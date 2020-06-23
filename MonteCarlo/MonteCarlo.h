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
    int step = 1;// ʱ�䲽��
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
    double iterate(); // ����һ�Σ�����һ��step�������ȴ����������¶ȸ���ֵ
    void run(int verbose=0); // ��ֹͣ��ֱ�ӵ����㷨��ֱ��������̬�򵽴�����������

	void runWithOneDim(int verbose = 0);
	void runWithRandomWalk(int verbose = 0);

	double computeTempWithRandomWalk(int x, int y, int z, int n);//RandomWalkģʽ�ĵ���λ�ü���
	std::vector<double> computeProbs(int x, int y, int z);

    void output1d();
    void output3d();
};

#endif

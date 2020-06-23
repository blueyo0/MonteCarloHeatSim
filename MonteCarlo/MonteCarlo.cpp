#include "MonteCarlo.h"
#include "GlobalFun.h"

#include <iomanip>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>

MonteCarlo::MonteCarlo(Shape* in_shape) : shape(in_shape) {
	center[0] = 10;
	center[1] = 10;
	center[2] = 10;
    Vector3d bd = shape->getBoundingBox();
    double x_step = 1/shape->getStep(0);
    int length = (int)sqrt(bd.x*bd.x+bd.y*bd.y+bd.z*bd.z)*x_step;
    temp_1d.resize(length, shape->getTemp(0,0,0));
	Vector3d size = shape->getSize();
	this->temp_3d = GlobalFun::create3DArray(size.x, size.y, size.z, 37.0);
}

void MonteCarlo::setIteration(int num)
{
	time_max = num;
}

void MonteCarlo::setProbe(int x, int y, int z, double T)
{
    center[0] = x;
    center[1] = y;
    center[2] = z;
    shape->setTemp(x,y,z,T);
    temp_1d[0] = T;
	temp_3d[x][y][z] = T;
}

void MonteCarlo::reset()
{
	for (int i = 1; i < temp_1d.size(); i++)
	{
		temp_1d[i] = default_value;
	}
	temp_3d = GlobalFun::create3DArray(shape->getSize().x, 
									   shape->getSize().y, 
									   shape->getSize().z, 
									   default_value);
}

void MonteCarlo::setDefaultValue(double v)
{
	default_value = v;
}

double MonteCarlo::iterate()
{
	double error = 0.0;
	if (this->mode == SimMode::ONE_DIM) {
		/*一维模拟三维*/
		std::vector<double> next_temp_1d(temp_1d.size());
		double lambda = shape->getAlpha(0, 0, 0);
		for (int i = 1; i < temp_1d.size() - 1; ++i) {
			next_temp_1d[i] = lambda * temp_1d[i + 1] + lambda * temp_1d[i - 1] + (1 - 2 * lambda)*temp_1d[i];
		}
		next_temp_1d[0] = temp_1d[0];
		next_temp_1d[next_temp_1d.size() - 1] = next_temp_1d[next_temp_1d.size() - 2];

		temp_1d.assign(next_temp_1d.begin(), next_temp_1d.end());
		int count = 0;
		for (int i = 0; i < temp_1d.size() - 1; ++i) {
			error += next_temp_1d[i] - temp_1d[i];
			count++;
		}
		error /= count;
	}
	else if (this->mode == SimMode::RANDOM_WALK) {
		/*随机游走算法的MonteCarlo三维模拟*/
		Vector3d size = shape->getSize();
		// TO-DO: 边界点温度变化计算
		for (int ix = 1; ix < size.x - 1; ++ix) {
			for (int iy = 1; iy < size.y - 1; ++iy) {
				for (int iz = 1; iz < size.z - 1; ++iz) {
					// 对所有的内部点，进行随机游走计算
					double next_T = computeTempWithRandomWalk(ix, iy, iz, this->monteCarloNum);
					this->temp_3d[ix][iy][iz] = next_T;
				}
			}
		}
	}
	
    return error;
}


void MonteCarlo::runWithOneDim(int verbose)
{
	for (int t = 0; t < time_max; t += step) {
		if (verbose == 1) {
			output1d();
		}
		else if (verbose == 3) {
			output3d();
		}
		iterate();
		Vector3d size = shape->getSize();
		for (int i = 0; i < size.x; ++i) {
			for (int j = 0; j < size.y; ++j) {
				for (int k = 0; k < size.z; ++k) {
					double dist = sqrt((i - center[0])*(i - center[0]) +
						(j - center[1])*(j - center[1]) +
						(k - center[2])*(k - center[2]));
					int lb = round(dist), hb = ceil(dist);
					double length = hb - lb;
					double next_T = temp_1d[lb];
					if (length > 0) {
						next_T = (temp_1d[lb] * abs(hb - dist) + temp_1d[hb] * abs(dist - lb)) / length;
					}
					if (verbose == 2) {
						std::cout.setf(std::ios::fixed);
						std::cout << std::setprecision(2) << next_T;
						std::cout.unsetf(std::ios::fixed);
						std::cout << "[" << dist << "]";
						std::cout << "(" << i << "," << j << "," << k << ")" << std::endl;
					}
					shape->setTemp(i, j, k, next_T);
				}
			}
		}
		if (verbose == 2) {
			std::cout << "=====================================================" << std::endl;
		}
	}
}


std::vector<double> MonteCarlo::computeProbs(int x, int y, int z)
{
	double W = 0.0005;
	double beta = 0.5;

	std::vector<double> res;
	double self_p = 1-W*(1-beta)*this->step;
	res.push_back(self_p);
	double neigh_p = this->step/(this->shape->getStep(0)*this->shape->getStep(0));
	std::vector<std::vector<int>> neigh;
	GlobalFun::getNeighPos({x,y,z}, neigh, shape->getSize().toVectorInt());
	for(std::vector<int> n : neigh){
		res.push_back(neigh_p*this->shape->getAlpha(n[0],n[1],n[2]));
	}
	return res;
}


//计算(x,y,z)位置, n次随机游走的结果
double MonteCarlo::computeTempWithRandomWalk(int x, int y, int z, int n)
{
	// 计算各个方向的概率
	std::vector<double> probs = this->computeProbs(x,y,z);
	std::vector<double> interval = GlobalFun::getProbsInterval(probs);
	
	// 开始迭代
	int N = (n < 1) ? 1 : n;
	double result = 0.0;
	std::vector<int> cur_pos = {x, y, z};
	srand((unsigned)time(NULL));
	for (int i = 0; i < N; ++i) {
		// TO-DO: 使用RW计算某个位置的温度	
		double rand_prob = rand() / double(RAND_MAX);

	}
	result /= N;
	return result;
}

void MonteCarlo::runWithRandomWalk(int verbose)
{
	Vector3d size = shape->getSize();
	if (verbose == 2) {
		std::cout << "====================[start]==========================" << std::endl;
	}
	// TO-DO: 检查差分稳定性，然后计算RW的概率：F, 1-W(1-beta)delta_t
	for (int t = 0; t < time_max; t += step) {
		this->curr_iterate_time = t+1;
		if (verbose == 3) {
			output3d();
		}
		iterate();
		for (int i = 0; i < size.x; ++i) {
			for (int j = 0; j < size.y; ++j) {
				for (int k = 0; k < size.z; ++k) {
					this->shape->setTemp(i, j, k, this->temp_3d[i][j][k]);
				}
			}
		}
	}
	if (verbose == 2) {
		std::cout << "======================[end]===========================" << std::endl;
	}

}

void MonteCarlo::run(int verbose)
{
	switch (this->mode)
	{
	case(SimMode::ONE_DIM):
		runWithOneDim(verbose);
		break;
	case(SimMode::RANDOM_WALK):
		runWithRandomWalk(verbose);
		break;
	default:
		break;
	}
} 

void MonteCarlo::output1d()
{
    for(int i=0; i<temp_1d.size()-1; ++i) {
           std::cout << temp_1d[i] << " ";
       }
       std::cout << std::endl << std::endl;
       std::cout << std::endl << std::endl;
}

void MonteCarlo::output3d()
{
    Vector3d size = shape->getSize();
    for(int i=0; i<size.x; ++i){
        for(int j=0; j<size.y; ++j){
            for(int k=0; k<size.z; ++k){
                std::cout << shape->getTemp(i,j,k) << " ";
            }
            std::cout << std::endl << std::endl;
        }  
        std::cout << std::endl << std::endl;  
    }    
    std::cout << std::endl << std::endl;
}


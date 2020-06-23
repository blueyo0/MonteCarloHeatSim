#pragma once
#include <vector>
#include <iostream>
#include <limits.h>


namespace GlobalFun{
	/*创建3维数组， 大小为(x,y,z)，初始值为T*/
    static double*** create3DArray(int x, int y, int z, int T){
        int size_with_step[3] = {x, y, z};
		double ***array = new double**[x];
        
        for(int i=0; i<size_with_step[0]; ++i){
			array[i] = new double*[size_with_step[1]];
		}
		for (int i = 0; i < size_with_step[0]; ++i) {
			for (int j = 0; j < size_with_step[1]; ++j) {
				array[i][j] = new double[size_with_step[2]];
			}
		}

		for (int i = 0; i < size_with_step[0]; ++i) {
			for (int j = 0; j < size_with_step[1]; ++j) {
				for (int k = 0; k < size_with_step[2]; ++k) {
					array[i][j][k] = 37.0;
                }
            }
        }
		return array;
    }

	/* 将概率数组划分到[0,1]区间，并给出边界点
	 * 例如：对于[1,1,2]，应当输出[0.00, 0.25, 0.50, 1.00]
	 */
	static std::vector<double> getProbsInterval(std::vector<double> in_prob){
		double sum = 0.0;
		for(double val : in_prob) sum += val;

		std::vector<double> interval;
		double prob_sum = 0.0;
		for(double val : in_prob) {
			interval.push_back(prob_sum);
			prob_sum += val/sum;
		}
		interval.push_back(prob_sum);
		return interval;
	}

	/* 检查点是否有效 */
	bool isPointValid(std::vector<int> pos, std::vector<int> boundary={INT_MAX, INT_MAX, INT_MAX}){
		int dim = pos.size();
		if(boundary.size()<dim) boundary.resize(dim, INT_MAX);
		bool flag  = true;
		for(int j=0; j<dim; ++j){
			if(pos[j] < 0 || pos[j] >= boundary[j]){
				flag = false;
			}
		}
		return flag;
	}

	/**
	 * 获取一个pos的neighbors的坐标
	*/
	bool getNeighPos(std::vector<int> pos, std::vector<std::vector<int>> &result,
					 std::vector<int> boundary={INT_MAX, INT_MAX, INT_MAX}){
		int dim = pos.size();
		for(int i=0; i<dim; ++i){
			std::vector<int> n_neigh(pos), p_neigh(pos);
			n_neigh[i]--;
			p_neigh[i]++;
						
			if(isPointValid(n_neigh, boundary)) result.push_back(n_neigh);
			if(isPointValid(p_neigh, boundary)) result.push_back(p_neigh);
		}
		return true;
	}

};




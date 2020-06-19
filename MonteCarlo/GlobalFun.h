#pragma once

namespace GlobalFun{
    static double*** create3DArray(int x, int y, int z, int T){
		double ***array;
        int size_with_step[3] = {x, y, z};
        
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
};


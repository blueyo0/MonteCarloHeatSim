#include "Shape.h"
#include <math.h>

Vector3d::Vector3d(double x, double y, double z)
    : x(x), y(y), z(z) {}


UniformCube::UniformCube()
{
    for(int i=0; i<size_with_step[0]; ++i){
        this->temp[i] = new double*[size_with_step[1]];
    }
    for(int i=0; i<size_with_step[0]; ++i){
        for(int j=0; j<size_with_step[1]; ++j){
            this->temp[i][j] = new double[size_with_step[2]];
        }
    }

    for(int i=0; i<size_with_step[0]; ++i){
        for(int j=0; j<size_with_step[1]; ++j){
            for(int k=0; k<size_with_step[2]; ++k){
                this->temp[i][j][k] = 37.0;
            }
        }
    }
}

UniformCube::UniformCube(Vector3d size, double step, double alpha, double default_temp)
    :step(step), alpha_factor(alpha)
{
    this->size_with_step[0] = ceil(size.x / step);
    this->size_with_step[1] = ceil(size.y / step);
    this->size_with_step[2] = ceil(size.z / step);
    this->temp = new double**[size_with_step[0]];
    for(int i=0; i<size_with_step[0]; ++i){
        this->temp[i] = new double*[size_with_step[1]];
    }
    for(int i=0; i<size_with_step[0]; ++i){
        for(int j=0; j<size_with_step[1]; ++j){
            this->temp[i][j] = new double[size_with_step[2]];
        }
    }
    for(int i=0; i<size_with_step[0]; ++i){
        for(int j=0; j<size_with_step[1]; ++j){
            for(int k=0; k<size_with_step[2]; ++k){
                this->temp[i][j][k] = default_temp;
            }
        }
    }
}

double UniformCube::getAlpha(int x, int y, int z)
{
    return this->alpha_factor;
}

double UniformCube::getTemp(int x, int y, int z)
{
    return this->temp[x][y][z];
}

double UniformCube::setTemp(int x, int y, int z, double T)
{
    this->temp[x][y][z] = T;
    return this->temp[x][y][z];
}

double UniformCube::getStep(int dim)
{
    return this->step;
}

Vector3d UniformCube::getBoundingBox()
{
    return Vector3d(
        this->size_with_step[0]*step,
        this->size_with_step[1]*step,
        this->size_with_step[2]*step
    );
}

Vector3d UniformCube::getSize()
{
    return  Vector3d(
        this->size_with_step[0],
        this->size_with_step[1],
        this->size_with_step[2]
    );
}



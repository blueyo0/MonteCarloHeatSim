#include "MonteCarlo.h"
#include <math.h>
#include <iostream>
#include <iomanip>

MonteCarlo::MonteCarlo(Shape* in_shape) : shape(in_shape) { 
    Vector3d bd = shape->getBoundingBox();
    double x_step = 1/shape->getStep(0);
    int length = (int)sqrt(bd.x*bd.x+bd.y*bd.y+bd.z*bd.z)*x_step;
    temp_1d.resize(length, shape->getTemp(0,0,0));
}

void MonteCarlo::setProbe(int x, int y, int z, double T)
{
    center[0] = x;
    center[1] = y;
    center[2] = z;
    shape->setTemp(x,y,z,T);
    temp_1d[0] = T;
}

double MonteCarlo::iterate()
{
    std::vector<double> next_temp_1d(temp_1d.size());
    // double lambda = shape->getAlpha(0,0,0)*step/(shape->getStep(0)*shape->getStep(0));
    double lambda = shape->getAlpha(0,0,0);//*step/(shape->getStep(0)*shape->getStep(0));
    for(int i=1; i<temp_1d.size()-1; ++i) {
        next_temp_1d[i] = lambda*temp_1d[i+1]+lambda*temp_1d[i-1]+(1-2*lambda)*temp_1d[i];
    }
    next_temp_1d[0] = temp_1d[0];
    next_temp_1d[next_temp_1d.size()-1] = next_temp_1d[next_temp_1d.size()-2];
    temp_1d.assign(next_temp_1d.begin(), next_temp_1d.end());
    return 0.0;
}

void MonteCarlo::run(int verbose)
{
    for(int t=0; t<time_max; t+=step){
        if(verbose==1){
            output1d();
        } 
        else if(verbose==3){
            output3d();
        }
        iterate();
        Vector3d size = shape->getSize();
        for(int i=0; i<size.x; ++i){
            for(int j=0; j<size.y; ++j){
                for(int k=0; k<size.z; ++k){
                    double dist = sqrt((i-center[0])*(i-center[0])+
                                       (j-center[1])*(j-center[1])+
                                       (k-center[2])*(k-center[2]));
                    int lb = round(dist), hb = ceil(dist);
                    double length = hb-lb;
                    double next_T = temp_1d[lb];
                    if(length>0){
                        next_T = (temp_1d[lb]*abs(hb-dist) + temp_1d[hb]*abs(dist-lb))/length;
                    }
                    if(verbose==2){
                        std::cout.setf(std::ios::fixed);
                        std::cout << std::setprecision(2) << next_T;
                        std::cout.unsetf(std::ios::fixed);
                        std::cout << "[" << dist << "]";
                        std::cout << "(" << i << "," << j << "," << k << ")" << std::endl;
                    }
                    shape->setTemp(i,j,k,next_T);
                }
            }    
        }
        if(verbose==2){    
            std::cout << "=====================================================" << std::endl;
        }
    }
} 

void MonteCarlo::output1d()
{
    for(int i=0; i<temp_1d.size()-1; ++i) {
        std::cout.setf(std::ios::fixed);
        std::cout << std::setprecision(2) << temp_1d[i] << " ";
        std::cout.unsetf(std::ios::fixed);
    }
    std::cout << std::endl << std::endl;
}

void MonteCarlo::output3d()
{
    Vector3d size = shape->getSize();
    for(int i=0; i<size.x; ++i){
        for(int j=0; j<size.y; ++j){
            for(int k=0; k<size.z; ++k){
                std::cout.setf(std::ios::fixed);
                std::cout << std::setprecision(2) << shape->getTemp(i,j,k) << " ";
                std::cout.unsetf(std::ios::fixed);
            }
            std::cout << std::endl << std::endl;
        }  
        std::cout << std::endl << std::endl;  
    }    
    std::cout << "=====================================================" << std::endl;
    std::cout << std::endl << std::endl;
}




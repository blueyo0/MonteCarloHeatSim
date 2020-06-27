#include "GlobalFun.h"
#include "Shape.h"
#include <iostream>
#include <vector>
#include <iomanip>

inline void test_getProbsInterval()
{
    std::vector<double> p;
    p.push_back(1);
    p.push_back(1);
    p.push_back(2);
    for (double val : p)
        std::cout << val << std::endl;
    std::vector<double> res = GlobalFun::getProbsInterval(p);
    std::cout << "hello world" << std::endl;
    for (double val : res)
        std::cout << val << std::endl;
}

inline void test_getNeighPos()
{
    std::vector<int> pos = {2,2,1};
    std::vector<std::vector<int>> res;
    if(GlobalFun::getNeighPos(pos, res, std::vector<int>({3,3,3}))){
        for(std::vector<int> pt : res){
            for(int idx : pt){
                std::cout << idx << ", ";
            }
            std::cout << std::endl;
        }
    }
}

inline void test_isPointBoundary()
{
    std::vector<int> pos = {2,2,1};
    std::cout << GlobalFun::isPointBoundary(pos, {50,26,50}) << std::endl;
    pos = {0,1,1};
    std::cout << GlobalFun::isPointBoundary(pos, {50,26,50}) << std::endl;
    pos = {1,0,1};
    std::cout << GlobalFun::isPointBoundary(pos, {50,26,50}) << std::endl;
    pos = {2,2,0};
    std::cout << GlobalFun::isPointBoundary(pos, {50,26,50}) << std::endl;
    pos = {48,24,48};
    std::cout << GlobalFun::isPointBoundary(pos, {50,26,50}) << std::endl;
    pos = {2,25,0};
    std::cout << GlobalFun::isPointBoundary(pos, {50,26,50}) << std::endl;
}


int main()
{
    // test_getNeighPos();
    // test_isPointBoundary();
    int time_cost = 3735;
    std::cout << "Elapsed time: " << std::setfill('0')
                << std::setw(2) << time_cost/3600 << ":"
                << std::setw(2) << (time_cost%3600)/60 << ":"
                << std::setw(2) << time_cost%60 << " s" << std::endl;
}



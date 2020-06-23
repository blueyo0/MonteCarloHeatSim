#include "GlobalFun.h"
#include "Shape.h"
#include <iostream>
#include <vector>

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


int main()
{
    test_getNeighPos();
}



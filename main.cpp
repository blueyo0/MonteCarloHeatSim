#include <iostream>
#include <math.h>
#include "Shape.h"
#include "MonteCarlo.h"

typedef UniformCube Ucube;

int main()
{
    Shape* s = new Ucube(Vector3d(1.1,1.2,1.3));
    MonteCarlo* algo = new MonteCarlo(s); 
    algo->setProbe(0, 0, 0, 60);
    algo->run(3);
}


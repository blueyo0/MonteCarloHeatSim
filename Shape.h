#pragma once
/*3维double类型向量*/
class Vector3d
{
public:
    double x,y,z;
    Vector3d(double x, double y, double z);
};

/*Shape接口*/
class Shape
{
public:
    virtual ~Shape() {}; 
    // 获取某个位置的热扩散率   
    virtual double getAlpha(int x, int y, int z) = 0;
    // 获取某个位置的温度
    virtual double getTemp(int x, int y, int z) = 0;
    // 设置某个位置的温度
    virtual double setTemp(int x, int y, int z, double T) = 0;
    // 获取某个维度的步长
    virtual double getStep(int dim) = 0;
    // 获取形状的BoundingBox的三维大小
    virtual Vector3d getBoundingBox() = 0;
    // 获取形状的网格的三维大小
    virtual Vector3d getSize() = 0;
};

/*均匀立方体*/
class UniformCube : public Shape
{
protected:
    double alpha_factor = 0.432; //热扩散率
    double step = 0.1; //各个维度使用统一步长
    int size_with_step[3] = {10,10,20}; //各个维度上的长度 = size_with_step*step
    double ***temp;
public:
    UniformCube();
    UniformCube(Vector3d size, double step=0.1, double alpha=0.432, double default_temp=37.0);
    virtual double getAlpha(int x, int y, int z) override;
    virtual double getTemp(int x, int y, int z) override;
    virtual double setTemp(int x, int y, int z, double T) override;
    virtual double getStep(int dim) override;
    virtual Vector3d getBoundingBox() override;
    virtual Vector3d getSize() override;
};






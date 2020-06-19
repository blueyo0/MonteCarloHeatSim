#pragma once
/*3άdouble��������*/
class Vector3d
{
public:
    double x,y,z;
    Vector3d(double x, double y, double z);
};

/*Shape�ӿ�*/
class Shape
{
public:
    virtual ~Shape() {}; 
    // ��ȡĳ��λ�õ�����ɢ��   
    virtual double getAlpha(int x, int y, int z) = 0;
    // ��ȡĳ��λ�õ��¶�
    virtual double getTemp(int x, int y, int z) = 0;
    // ����ĳ��λ�õ��¶�
    virtual double setTemp(int x, int y, int z, double T) = 0;
    // ��ȡĳ��ά�ȵĲ���
    virtual double getStep(int dim) = 0;
    // ��ȡ��״��BoundingBox����ά��С
    virtual Vector3d getBoundingBox() = 0;
    // ��ȡ��״���������ά��С
    virtual Vector3d getSize() = 0;
};

/*����������*/
class UniformCube : public Shape
{
protected:
    double alpha_factor = 0.432; //����ɢ��
    double step = 1; //����ά��ʹ��ͳһ����
	int size_with_step[3]; //����ά���ϵĳ��� = size_with_step*step
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

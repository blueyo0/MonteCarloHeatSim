#include "MainWindow.h"

#include "QFileDialog"
#include "vtkLookupTable.h"

#include "MonteCarlo/Shape.h"
#include "MonteCarlo/MonteCarlo.h"


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	// 初始化窗口对象
	ui.setupUi(this);
	therapyControl = new TherapyControl(this);
	ui_therapy = therapyControl->GetUI();

	// 设置VTK窗口的一些属性
	Renwin_3D = ui.qvtkWidget_3D->GetRenderWindow();
	Ren_3D = vtkRenderer::New();
	Renwin_3D->AddRenderer(Ren_3D);
	Renwin_3D->SetMultiSamples(0);
	Renwin_3D->SetAlphaBitPlanes(1);
	Ren_3D->SetUseDepthPeeling(1);

	
	const float constTemp = 37.0;	// 常温参数
	const float sourceTemp = 80.0;	// 热源温度参数
	// 显示三维数据的VTK数据对象
	VolumeRegionData = vtkImageData::New();
	VolumeRegionData->SetSpacing(1, 1, 1);				// 数据间隔
	VolumeRegionData->SetOrigin(-25, -25, -25);			// 数据起始位置
	VolumeRegionData->SetDimensions(50, 26, 50);	// 数据维度大小
	VolumeRegionData->AllocateScalars(VTK_FLOAT, 1);			// 分配内存
	// auto volumeRegionDataPtr = (float*)VolumeRegionData->GetScalarPointer();	// 数据指针，可直接修改，修改后需要调用Modified()函数保证VTK对象被更新

	// 构建模拟热场计算对象
	VolumeShape = new UniformCube(Vector3d(50, 26, 50), 1);
	MCAlgro = new MonteCarlo(VolumeShape);
	MCAlgro->setProbe(25, 25, 49, sourceTemp);
	MCAlgro->setDefaultValue(constTemp);
	MCAlgro->reset();
	MCAlgro->setIteration(5);						// 设置迭代计算的次数
	MCAlgro->run();

	this->UpdateDataFromMC();

	// 配置VTK数据显示的pipeline
	VolumeRegion = svrVolumeRegionActor::New();
	VolumeRegion->SetScalarRange(constTemp, sourceTemp);
	VolumeRegion->UserDataEnableOn();
	VolumeRegion->SetUserImageData(VolumeRegionData);
	Ren_3D->AddVolume(VolumeRegion);

	// 绑定QT的事件响应
	connect(ui.actionTherapyControl, SIGNAL(triggered()), this, SLOT(onOpenTherapyWin()));
	connect(this, SIGNAL(update3DRender()), this, SLOT(onUpdate3DRender()), Qt::BlockingQueuedConnection);
	connect(ui_therapy->horizontalSlider_factor, SIGNAL(valueChanged(int)), this, SLOT(onUpdateNeedleRegionFactor(int)));
}

MainWindow::~MainWindow()
{
	// 析构
}

void MainWindow::onOpenTherapyWin()
{
	// 以非模态窗口方式显示参数控制窗口
	therapyControl->show();
}

void MainWindow::onUpdate3DRender()
{
	// 更新显示VTK窗口
	this->Renwin_3D->Render();
}

void MainWindow::onUpdateNeedleRegionFactor(int v)
{
	// 修改参数后重新计算
	this->MCAlgro->reset();
	this->MCAlgro->setIteration(v);
	this->MCAlgro->run();
	this->UpdateDataFromMC();
	// 获取界面中个参数数值
	ui_therapy->label_factorValue->setText(QString::number(v));
	
	this->VolumeRegion->Modified();
	// 刷新窗口显示
	onUpdate3DRender();
}

void MainWindow::UpdateDataFromMC()
{
	auto dims = VolumeRegionData->GetDimensions();
	auto volumeRegionDataPtr = (float*)VolumeRegionData->GetScalarPointer();
	for (int iz = 0; iz < dims[2]; iz++)
	{
		for (int iy = 0; iy < dims[1]; iy++)
		{
			for (int ix = 0; ix < dims[0]; ix++)
			{
				int idx = ix + iy*dims[0] + iz*dims[1] * dims[0];
				volumeRegionDataPtr[idx] = VolumeShape->getTemp(ix, iy, iz);
			}
		}
	}
	VolumeRegionData->Modified();
}

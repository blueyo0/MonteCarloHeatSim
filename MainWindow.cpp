#include "MainWindow.h"

#include "QFileDialog"
#include "vtkLookupTable.h"

#include "MonteCarlo/Shape.h"
#include "MonteCarlo/MonteCarlo.h"


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	// ��ʼ�����ڶ���
	ui.setupUi(this);
	therapyControl = new TherapyControl(this);
	ui_therapy = therapyControl->GetUI();

	// ����VTK���ڵ�һЩ����
	Renwin_3D = ui.qvtkWidget_3D->GetRenderWindow();
	Ren_3D = vtkRenderer::New();
	Renwin_3D->AddRenderer(Ren_3D);
	Renwin_3D->SetMultiSamples(0);
	Renwin_3D->SetAlphaBitPlanes(1);
	Ren_3D->SetUseDepthPeeling(1);

	
	const float constTemp = 37.0;	// ���²���
	const float sourceTemp = 80.0;	// ��Դ�¶Ȳ���
	// ��ʾ��ά���ݵ�VTK���ݶ���
	VolumeRegionData = vtkImageData::New();
	VolumeRegionData->SetSpacing(1, 1, 1);				// ���ݼ��
	VolumeRegionData->SetOrigin(-25, -25, -25);			// ������ʼλ��
	VolumeRegionData->SetDimensions(50, 26, 50);	// ����ά�ȴ�С
	VolumeRegionData->AllocateScalars(VTK_FLOAT, 1);			// �����ڴ�
	// auto volumeRegionDataPtr = (float*)VolumeRegionData->GetScalarPointer();	// ����ָ�룬��ֱ���޸ģ��޸ĺ���Ҫ����Modified()������֤VTK���󱻸���

	// ����ģ���ȳ��������
	VolumeShape = new UniformCube(Vector3d(50, 26, 50), 1);
	MCAlgro = new MonteCarlo(VolumeShape);
	MCAlgro->setProbe(25, 25, 49, sourceTemp);
	MCAlgro->setDefaultValue(constTemp);
	MCAlgro->reset();
	MCAlgro->setIteration(5);						// ���õ�������Ĵ���
	MCAlgro->run();

	this->UpdateDataFromMC();

	// ����VTK������ʾ��pipeline
	VolumeRegion = svrVolumeRegionActor::New();
	VolumeRegion->SetScalarRange(constTemp, sourceTemp);
	VolumeRegion->UserDataEnableOn();
	VolumeRegion->SetUserImageData(VolumeRegionData);
	Ren_3D->AddVolume(VolumeRegion);

	// ��QT���¼���Ӧ
	connect(ui.actionTherapyControl, SIGNAL(triggered()), this, SLOT(onOpenTherapyWin()));
	connect(this, SIGNAL(update3DRender()), this, SLOT(onUpdate3DRender()), Qt::BlockingQueuedConnection);
	connect(ui_therapy->horizontalSlider_factor, SIGNAL(valueChanged(int)), this, SLOT(onUpdateNeedleRegionFactor(int)));
}

MainWindow::~MainWindow()
{
	// ����
}

void MainWindow::onOpenTherapyWin()
{
	// �Է�ģ̬���ڷ�ʽ��ʾ�������ƴ���
	therapyControl->show();
}

void MainWindow::onUpdate3DRender()
{
	// ������ʾVTK����
	this->Renwin_3D->Render();
}

void MainWindow::onUpdateNeedleRegionFactor(int v)
{
	// �޸Ĳ��������¼���
	this->MCAlgro->reset();
	this->MCAlgro->setIteration(v);
	this->MCAlgro->run();
	this->UpdateDataFromMC();
	// ��ȡ�����и�������ֵ
	ui_therapy->label_factorValue->setText(QString::number(v));
	
	this->VolumeRegion->Modified();
	// ˢ�´�����ʾ
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

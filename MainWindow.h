#pragma once

#include "ui_MainWin.h"

// VTK��ʾ���ڵĺ�����
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkImageData.h"
#include "svrVolumeRegionActor.h"

// ���Ʋ����Ĵ���
#include "TherapyControl.h"
#include "ui_TherapyWin.h"

// monte carlo����
#include "MonteCarlo/Shape.h"
#include "MonteCarlo/MonteCarlo.h"


// QT������
class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();
	private slots:
	// �򿪲������ƴ���
	void onOpenTherapyWin();
	// ˢ��3D��ʾ����
	void onUpdate3DRender();
	// ˢ�²���ֵ���㲢��ʾ
	void onUpdateNeedleRegionFactor(int v);

protected:

private:
	TherapyControl* therapyControl;
	Ui::MainWinUI ui;
	Ui::TherapyWin* ui_therapy;

	vtkRenderWindow *Renwin_3D;
	vtkRenderer *Ren_3D;

	vtkImageData* VolumeRegionData;
	svrVolumeRegionActor* VolumeRegion;
	Shape* VolumeShape;
	MonteCarlo* MCAlgro;

	// ��monte carlo������������VTK���ݶ������Ա�ˢ����ʾ
	void UpdateDataFromMC();
	
signals:
	void update3DRender();
};


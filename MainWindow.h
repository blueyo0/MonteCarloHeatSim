#pragma once

#include "ui_MainWin.h"

// VTK显示窗口的核心类
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkImageData.h"
#include "svrVolumeRegionActor.h"

// 控制参数的窗口
#include "TherapyControl.h"
#include "ui_TherapyWin.h"

// monte carlo代码
#include "MonteCarlo/Shape.h"
#include "MonteCarlo/MonteCarlo.h"


// QT主窗口
class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0);
	~MainWindow();
	private slots:
	// 打开参数控制窗口
	void onOpenTherapyWin();
	// 刷新3D显示窗口
	void onUpdate3DRender();
	// 刷新参数值计算并显示
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

	// 将monte carlo计算结果拷贝到VTK数据对象中以便刷新显示
	void UpdateDataFromMC();
	
signals:
	void update3DRender();
};


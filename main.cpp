/*!
* \file 		main.cpp
* \brief
*
*本程序中包含了简单的VTK显示和Monte Carlo代码
*
*\author 		Wang Tao(wangtao@accu-med.cn)
*\version 		V1.0
*\date 			2020/06/11
*/

#include "MainWindow.h"
#include <QApplication>
#include <QTranslator>

#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkRenderingFreeType)
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	QTranslator qTranslator;
	qTranslator.load(":/language/language_zh.qm");
	a.installTranslator(&qTranslator);

	MainWindow w;
	w.show();
	return a.exec();
}
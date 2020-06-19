#include "TherapyControl.h"

TherapyControl::TherapyControl(QWidget* parent) :QDialog(parent)
{
	ui.setupUi(this);
}

TherapyControl::~TherapyControl()
{

}

Ui::TherapyWin* TherapyControl::GetUI()
{
	return &ui;
}

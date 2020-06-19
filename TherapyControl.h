
#pragma once

#include <QDialog>

#include "ui_TherapyWin.h"

class TherapyControl : public QDialog
{
	Q_OBJECT
public:
	TherapyControl(QWidget* parent = 0);
	~TherapyControl();

	Ui::TherapyWin* GetUI();
protected:

private:
	Ui::TherapyWin ui;
};
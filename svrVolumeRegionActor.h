#pragma once 

#include "vtkVolume.h"
#include <vtkImageData.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include "vtkLookupTable.h"
#include "vtkScalarBarActor.h"

class svrVolumeRegionActor : public vtkVolume
{
public:
	static svrVolumeRegionActor* New();
	vtkTypeMacro(svrVolumeRegionActor, vtkVolume);
	void PrintSelf(ostream& os, vtkIndent indent) override;

	double* GetBounds() override;
	int RenderVolumetricGeometry(vtkViewport *viewport) override;
	int RenderOpaqueGeometry(vtkViewport *viewport) override;
	int RenderOverlay(vtkViewport *) override;
	void ReleaseGraphicsResources(vtkWindow *win) override;
	// int HasTranslucentPolygonalGeometry() override;

	//void SetImageData(vtkImageData* args);
	vtkGetObjectMacro(ImageData, vtkImageData);

	void SetSpacing(double _arg1, double _arg2, double _arg3);
	void SetSpacing(double _arg[3]){
		this->SetSpacing(_arg[0], _arg[1], _arg[2]);
	}
	vtkGetVector3Macro(Spacing, double);

	void SetCenter(double _arg1, double _arg2, double _arg3);
	void SetCenter(double _arg[3]){
		this->SetCenter(_arg[0], _arg[1], _arg[2]);
	}
	vtkGetVector3Macro(Center, double);

	vtkGetMacro(Radius, double);
	vtkSetMacro(Radius, double);

	vtkGetMacro(LogFactor, int);
	vtkSetMacro(LogFactor, int);

	//void SetColors(vtkColorTransferFunction* args);
	vtkGetObjectMacro(Colors, vtkColorTransferFunction);

	//void SetOpacities(vtkPiecewiseFunction* args);
	vtkGetObjectMacro(Opacities, vtkPiecewiseFunction);

	vtkGetObjectMacro(LookupTable, vtkLookupTable);
	vtkGetObjectMacro(ScalarBar, vtkScalarBarActor);

	vtkSetObjectMacro(UserImageData, vtkImageData);
	
	vtkGetMacro(UserDataEnable, int);
	vtkSetMacro(UserDataEnable, int);
	vtkBooleanMacro(UserDataEnable, int);

	vtkGetVector2Macro(ScalarRange, double);
	vtkSetVector2Macro(ScalarRange, double);
	
protected:
	svrVolumeRegionActor();
	~svrVolumeRegionActor();

	virtual void BuildSelf();

	vtkTimeStamp BuildTime;

	int UserDataEnable;
	vtkImageData* UserImageData;
	vtkImageData* ImageData;

	vtkLookupTable* LookupTable;
	vtkScalarBarActor* ScalarBar;

	double ScalarRange[2];
	double Spacing[3];
	double Radius;
	int LogFactor;

	vtkColorTransferFunction* Colors;
	vtkPiecewiseFunction* Opacities;
private:
	// not used function
	svrVolumeRegionActor(const svrVolumeRegionActor&) VTK_DELETE_FUNCTION;
	void operator=(const svrVolumeRegionActor&)VTK_DELETE_FUNCTION;
};

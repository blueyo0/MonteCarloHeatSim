#include "svrVolumeRegionActor.h"
#include "vtkObjectFactory.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkVolumeProperty.h"
#include "vtkGPUVolumeRayCastMapper.h"
#include "vtkMath.h"
#include "vtkLookupTable.h"
#include "vtkTransform.h"

vtkStandardNewMacro(svrVolumeRegionActor);

void svrVolumeRegionActor::PrintSelf(ostream& os, vtkIndent indent)
{
	Superclass::PrintSelf(os, indent);
}

double* svrVolumeRegionActor::GetBounds()
{
	this->BuildSelf();
	if (UserDataEnable)
	{
		if (this->UserImageData)
		{
			this->UserImageData->GetBounds(Bounds);
		}
		else
		{

		}
	}
	else
	{	
		if (this->ImageData)
		{
			this->ImageData->GetBounds(Bounds);
		}
		else
		{
		
		}
	}
	return Bounds;
}

int svrVolumeRegionActor::RenderVolumetricGeometry(vtkViewport* viewport)
{
	this->BuildSelf();

	int renderNum = 0;
	if (this->Visibility)
		renderNum += Superclass::RenderVolumetricGeometry(viewport);

	return renderNum;
}

int svrVolumeRegionActor::RenderOpaqueGeometry(vtkViewport *viewport)
{
	this->BuildSelf();

	int renderNum = 0;
	if (this->Visibility){
		renderNum += Superclass::RenderOpaqueGeometry(viewport);
		renderNum += this->ScalarBar->RenderOpaqueGeometry(viewport);
	}

	return renderNum;
}

int svrVolumeRegionActor::RenderOverlay(vtkViewport* viewport)
{
	int num = Superclass::RenderOverlay(viewport);
	num += this->ScalarBar->RenderOverlay(viewport);
	return num;
}

void svrVolumeRegionActor::ReleaseGraphicsResources(vtkWindow* win)
{
	Superclass::ReleaseGraphicsResources(win);
	this->ScalarBar->ReleaseGraphicsResources(win);
}

svrVolumeRegionActor::svrVolumeRegionActor()
{
	this->ImageData = vtkImageData::New();
	this->UserImageData = vtkImageData::New();
	
	this->Mapper = vtkGPUVolumeRayCastMapper::New();
	vtkGPUVolumeRayCastMapper* mapper = (vtkGPUVolumeRayCastMapper*)(this->Mapper);
	this->Mapper->SetInputDataObject(this->ImageData);

	this->GetProperty();
	LookupTable = vtkLookupTable::New();
	ScalarBar = vtkScalarBarActor::New();
	Colors = vtkColorTransferFunction::New();
	Opacities = vtkPiecewiseFunction::New();
	//this->SetImageData(this->ImageData);
	this->ScalarRange[0] = 37.0;
	this->ScalarRange[1] = 80.0;
	LookupTable->SetHueRange(0.7, 0);
	LookupTable->SetAlphaRange(0.3, 1);
	LookupTable->SetValueRange(1.0, 1.0);
	LookupTable->SetSaturationRange(1.0, 1.0);
	LookupTable->SetNumberOfTableValues(ScalarRange[1] - ScalarRange[0]);
	LookupTable->SetRange(ScalarRange);
	LookupTable->Build();
	LookupTable->SetTableValue(0, 0, 0, 0, 0);

	ScalarBar->SetLookupTable(LookupTable);
	ScalarBar->SetTitle("C");
	ScalarBar->SetNumberOfLabels(6); //ÉèÖÃ5¸ö±êÇ© 
	ScalarBar->SetUnconstrainedFontSize(1);
	ScalarBar->SetPosition(0.95, 0.0);
	ScalarBar->SetPosition2(0.05, 0.2);

	auto rangeSize = ScalarRange[1] - ScalarRange[0];
	for (int i = 0; i < rangeSize; i++)
	{
		auto idx = ScalarRange[0] + i / rangeSize;
		Colors->AddRGBPoint(idx, LookupTable->GetTableValue(i)[0], LookupTable->GetTableValue(i)[1], LookupTable->GetTableValue(i)[2]);
		Opacities->AddPoint(idx, LookupTable->GetTableValue(i)[3]);
	}
	
	this->Property->SetColor(Colors);
	this->Property->SetScalarOpacity(Opacities);
	this->Property->ShadeOn();
	this->Property->SetInterpolationTypeToLinear();
	this->Property->SetAmbient(0.2);
	this->Property->SetDiffuse(0.8);
	this->Property->SetSpecular(0.1);
	this->Property->SetSpecularPower(1);

	this->Center[0] = 0;
	this->Center[1] = 0;
	this->Center[2] = 0;

	this->Radius = 50;
	this->LogFactor = 1;
	this->Origin[0] = 0;
	this->Origin[1] = 0;
	this->Origin[2] = 0;
	this->Spacing[0] = 1;
	this->Spacing[1] = 1;
	this->Spacing[2] = 1;

	this->BuildTime.Modified();
	this->Modified();

	this->GetProperty();
}

svrVolumeRegionActor::~svrVolumeRegionActor()
{
	if (this->ImageData)
		this->ImageData->UnRegister(this);

}

void svrVolumeRegionActor::BuildSelf()
{
	if (this->BuildTime < this->MTime){

		LookupTable->SetNumberOfTableValues(ScalarRange[1] - ScalarRange[0]);
		LookupTable->SetRange(ScalarRange);
		LookupTable->Build();
		LookupTable->SetTableValue(0, 0, 0, 0, 0);
		auto rangeSize = ScalarRange[1] - ScalarRange[0];
		for (int i = 0; i < rangeSize; i++)
		{
			auto idx = ScalarRange[0] + i / rangeSize;
			Colors->AddRGBPoint(idx, LookupTable->GetTableValue(i)[0], LookupTable->GetTableValue(i)[1], LookupTable->GetTableValue(i)[2]);
			Opacities->AddPoint(idx, LookupTable->GetTableValue(i)[3]);
		}
		this->Property->Modified();
		this->ScalarBar->Modified();
		
		if (UserDataEnable)
		{
			this->Mapper->SetInputDataObject(this->UserImageData);
		}
		else{
			this->Mapper->SetInputDataObject(this->ImageData);

			
			// compute origins
			this->Origin[0] = Center[0] - Radius;
			this->Origin[1] = Center[1] - Radius;
			this->Origin[2] = Center[2] - Radius;
			// compute extent
			int dims[3] = {
				(Radius * 2) / Spacing[0],
				(Radius * 2) / Spacing[1],
				(Radius * 2) / Spacing[2]
			};

			this->ImageData->SetOrigin(this->Origin);
			this->ImageData->SetSpacing(this->Spacing);
			this->ImageData->SetDimensions(dims);
			this->ImageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
			auto dataPtr = (unsigned char*)ImageData->GetScalarPointer();

			memset(dataPtr, 0, sizeof(unsigned char)*dims[0] * dims[1] * dims[2]);

			if (dims[0] < 0 || dims[1] < 0 || dims[2] < 0){
				vtkErrorMacro(<< "Bad extent set for image data");
			}

			// generate image volume data
			this->ImageData->ComputeBounds();
			double bounds[6];
			ImageData->GetBounds(bounds);
			double center[3] = {
				(bounds[0] + bounds[1]) / 2,
				(bounds[2] + bounds[3]) / 2,
				(bounds[4] + bounds[5]) / 2,
			};

			for (double iz = center[2] - Radius; iz <= center[2] + Radius; iz += Spacing[2])
			{
				for (double iy = center[1] - Radius; iy <= center[1]/* + radius*/; iy += Spacing[1])
				{
					for (double ix = center[0] - Radius; ix <= center[0] + Radius; ix += Spacing[0])
					{
						if (ix<bounds[0] || ix>bounds[1] ||
							iy<bounds[2] || iy>bounds[3] ||
							iz<bounds[4] || iz>bounds[5])
							continue;

						// sphere shape
						double distance = sqrt(pow(center[0] - ix, 2) + pow(center[1] - iy, 2) + pow(center[2] - iz, 2));
						if (distance >= Radius)
							continue;

						const int idxX = (ix - Origin[0] + Spacing[0] / 2) / Spacing[0];
						const int idxY = (iy - Origin[1] + Spacing[1] / 2) / Spacing[1];
						const int idxZ = (iz - Origin[2] + Spacing[2] / 2) / Spacing[2];
						const int idx = idxX + idxY*dims[0] + idxZ*dims[0] * dims[1];
						double data = (1 - (distance / Radius));
						double dataScaled = std::log(1 + LogFactor*data) / std::log(1 + LogFactor);
						dataPtr[idx] = dataScaled * 255;
					}
				}
			}
			this->ImageData->Modified();
		}

		this->BuildTime = this->MTime;
	}
}

void svrVolumeRegionActor::SetSpacing(double _arg1, double _arg2, double _arg3)
{
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting Spacing to (" << _arg1 << "," << _arg2 << "," << _arg3 << ")");
	if ((this->Spacing[0] != _arg1) || (this->Spacing[1] != _arg2) || (this->Spacing[2] != _arg3))
	{
		this->Spacing[0] = _arg1;
		this->Spacing[1] = _arg2;
		this->Spacing[2] = _arg3;
		this->Modified();
	}
}

void svrVolumeRegionActor::SetCenter(double _arg1, double _arg2, double _arg3)
{
	vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting Center to (" << _arg1 << "," << _arg2 << "," << _arg3 << ")");
	if ((this->Center[0] != _arg1) || (this->Center[1] != _arg2) || (this->Center[2] != _arg3))
	{
		this->Center[0] = _arg1;
		this->Center[1] = _arg2;
		this->Center[2] = _arg3;
		this->Modified();
	}
}
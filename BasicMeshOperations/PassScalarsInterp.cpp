/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Copies scalars from one mesh to another. When the source is sparse, for every target point it finds a corresponding cell and interpolates the point scalars
*/
#include "vtkPolyDataReader.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkShortArray.h"
#include "vtkType.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkCell.h"
#include "VTKCommonTools.h"
#include "vtkFloatArray.h"


void PassScalarsFloat( int argc, char *argv[] )
{	

	std::cout << " reading source " << argv[1] << std::endl;
	vtkSmartPointer<vtkPolyData> sourcePd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[1]) );

	std::cout << " reading target " << argv[2] << std::endl;
	vtkSmartPointer<vtkPolyData> targetPd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[2]) );


	char* property_name = argv[3];
	std::cout << " property name: " << property_name << std::endl;


	std::cout<<"Interpolating pointdata"<<std::endl;

	vtkFloatArray *target_scalars = (vtkFloatArray*) targetPd -> GetPointData() -> GetArray(property_name);
	vtkFloatArray *source_scalars = (vtkFloatArray*) sourcePd -> GetPointData() -> GetArray(property_name);

	vtkSmartPointer<vtkCellLocator> in_source_locator = vtkSmartPointer<vtkCellLocator>::New();
	in_source_locator->SetDataSet(sourcePd);
	in_source_locator->BuildLocator();

	for (int i=0; i<targetPd -> GetNumberOfPoints(); i++) 
	{
		double *pt = targetPd -> GetPoint(i);

		double closestpoint[3];
		vtkIdType cellid;
		int subid;
		double dist2;
		in_source_locator->FindClosestPoint(pt,closestpoint,cellid,subid, dist2);

		//get the point weights for it's position inside the found cell
		double pcoords[3];
		double weights[3];
		double rubbish[3];
		sourcePd->GetCell(cellid)->EvaluatePosition(closestpoint,rubbish,subid,pcoords,dist2,weights);

		int pt0 = sourcePd->GetCell(cellid)->GetPointId(0);
		int pt1 = sourcePd->GetCell(cellid)->GetPointId(1);
		int pt2 = sourcePd->GetCell(cellid)->GetPointId(2);

		float v0 = source_scalars->GetValue(pt0);
		float v1 = source_scalars->GetValue(pt1);
		float v2 = source_scalars->GetValue(pt2);

		target_scalars->SetValue(i, weights[0]*v0 + weights[1]*v1 + weights[2]*v2);
	}



	char outputFileName[255];
	strcpy(outputFileName,argv[4]);
	std::cout << " save output to " << outputFileName << std::endl;

	CommonTools::SaveShapeToFile(targetPd,outputFileName);

}



int main( int argc, char *argv[] )
{
	
	if( argc < 5 )
	{
		std::cout << "When the source is sparse, for every target point it finds a corresponding cell and interpolates the point scalars." << std::endl;
		std::cout << "Usage (float point arrays only!): " << std::endl;
		std::cout << argv[0] << " <source> <target> <property name> <output>" << std::endl;
		return EXIT_FAILURE;
	}


	PassScalarsFloat( argc, argv );

	return 0;

}

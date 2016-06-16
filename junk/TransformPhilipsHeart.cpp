/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Transforms Philips mesh to suit EM simulations. Don't remember WTF is that
*/
#include "CommonTools.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkCleanPolyData.h>

//------------------------------------------------------------------
int main( int argc, char *argv[] )
{	

	if( argc < 2 )
	{
		std::cout << "Transforms Philips mesh to suit EM simulations v1.0 " << std::endl;
		return -1;
	}

	int c = 1;
	char* fullheart_filename = argv[c++];
	char* rv_filename = argv[c++];
	char* myo_filename = argv[c++];
	char* lv_filename = argv[c++];
	char* outshape_filename = argv[c++];
	
	vtkSmartPointer<vtkPolyData> fullheart = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(fullheart_filename) );
	vtkSmartPointer<vtkPolyData> rv = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(rv_filename) );
	vtkSmartPointer<vtkPolyData> myo = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(myo_filename) );
	vtkSmartPointer<vtkPolyData> lv = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(lv_filename) );
		
	///////////////////////////////////////////////////////////////////////////
	// 
	// body
	//
	
	vtkSmartPointer<vtkCleanPolyData> clean1 = vtkSmartPointer<vtkCleanPolyData>::New();
	clean1->SetInput( myo );
	clean1->Update();
		
	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booloper = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
	booloper->SetOperationToUnion();
	booloper->SetInput(0, clean1->GetOutput());
	booloper->SetInput(1, lv);
	booloper->SetTolerance(1e-1);
	booloper->Update();
	
	
	/////////////////////////////////////////////////////////////////////////////	
	CommonTools::SaveShapeToFile( booloper->GetOutput(), outshape_filename );
	
	return 0;
}

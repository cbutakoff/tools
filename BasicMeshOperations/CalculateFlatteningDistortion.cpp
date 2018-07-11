/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Given mesh and flattening, calculates distortion
*/
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkType.h"
#include "vtkSmartPointer.h"
#include "vtkCell.h"
#include "VTKCommonTools.h"

int main( int argc, char *argv[] )
{	

	if( argc < 5 )
	{
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0] << " <original_mesh> <flattened_mesh> <output_flattened_mesh>" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << " reading source " << argv[1] << std::endl;
	vtkSmartPointer<vtkPolyData> sourcePd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[1]) );

	std::cout << " reading target " << argv[2] << std::endl;
	vtkSmartPointer<vtkPolyData> targetPd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[2]) );


	CommonTools::SaveShapeToFile(targetPd,outputFileName);

	return 0;
}

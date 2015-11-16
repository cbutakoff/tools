/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief adds scalar array ro polydata and sets its value to a constant
*/
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include "vtkCellData.h"
#include "CommonTools.h"
#include "vtkSmartPointer.h"
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include "CommonTools.h"


template <typename vtk_array_type>
void AddArray( vtkPolyData* shapePt, bool use_points, double value, const char* property_name, vtk_array_type* fakeparam )
{
	cout << " property name: " << property_name << endl;
	cout << " scalar value: " << value << endl;

		
	vtkSmartPointer<vtk_array_type> regions = vtkSmartPointer<vtk_array_type>::New(); 
	
	
	regions -> SetName(property_name);

	if(use_points)
	{
		unsigned int npoints = shapePt->GetNumberOfPoints();
		regions->SetNumberOfValues(npoints);

		for (unsigned int i=0; i<npoints; i++) 
		{
			regions -> SetValue(i,value);
		}

		shapePt -> GetPointData() -> AddArray(regions);
	}
	else
	{
		unsigned int ncells = shapePt->GetNumberOfCells();
		regions->SetNumberOfValues(ncells);

		for (unsigned int i=0; i<ncells; i++) 
		{
			regions -> SetValue(i,value);
		}

		shapePt -> GetCellData() -> AddArray(regions);
	}

	//shapePt -> Update();
	
} 


//------------------------------------------------------------------
int main( int argc, char *argv[] )
{	

	if( argc < 5 )
	{
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0] << "SetScalars <shape> <property name> <value> <outputshape> <points|cells> <short|float>" << std::endl;
		return EXIT_FAILURE;
	}


	int c=1;
	const char *filenamein = argv[c++];
	const char *property = argv[c++];
	const char *c_value = argv[c++];
	const char *filenameout = argv[c++];
	const char *c_pointcell = argv[c++];
	const char *datatype = argv[c++];


	vtkSmartPointer<vtkPolyData> shape = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(filenamein) );


	cout << " opening " << filenamein << endl;

	bool use_points = strcmp(c_pointcell,"points")==0;
	double value = atof(c_value);

	if( strcmp(datatype,"short")==0 )
		AddArray( shape, use_points, value, property, (vtkShortArray*)NULL );
	else if( strcmp(datatype,"float")==0 )
		AddArray( shape, use_points, value, property, (vtkFloatArray*)NULL );
	else
	{
		std::cout<<"Unknown data type: "<<datatype<<std::endl;
		exit(-1);
	}

	CommonTools::SaveShapeToFile(shape,filenameout);

	return EXIT_SUCCESS;
}

/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Something to do with biventricular model generation
*/
#include <vtkPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>

#include <vtkImageData.h>
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkCellArray.h>

#include <vtkSmartPointer.h>

#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "CommonTools.h"


#include <vtkDelaunay3D.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkThreshold.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCellData.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

vtkPolyData* ExtractSurface( vtkUnstructuredGrid *volmesh, int id );

int main(int argc, char **argv)
{
	const int EPI_ID = 2;
	const int RV_ID = 3;
	const int LV_ID = 1;
	
	
	std::cout << "Version 0.99"<< std::endl;

	if( argc<2 ) 
	{
		std::cout << "Params: "<< std::endl;
		std::cout << "-shape <shape.vtk> \t\t- shape"<< std::endl;
		std::cout << "-image <image.vtk> \t\t- image"<< std::endl;
		std::cout << "-out <shape.vtk> \t\t- resulting image"<< std::endl;
		
		return -1;
	}

	char *inshapefile=NULL;
	char *inimagefile=NULL;
	char *outshapefile=NULL;

	
	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-shape")==0 )
		{
			inshapefile = argv[++c];
		}
		else if( strcmp(argv[c],"-image")==0 )
		{
			inimagefile = argv[++c];
		}
		else if( strcmp(argv[c],"-out")==0 )
		{
			outshapefile = argv[++c];
		}

	}


	std::cout<<"Shape: "<<inshapefile<<std::endl;
	std::cout<<"Image: "<<inimagefile<<std::endl;
	std::cout<<"Output: "<<outshapefile<<std::endl;
	

	
	vtkSmartPointer<vtkPolyData> inshape = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(inshapefile));

	vtkSmartPointer<vtkImageData> inimage = vtkSmartPointer<vtkImageData>::Take(
		CommonTools::LoadImage( inimagefile ));


	//extract LV and RV
	vtkSmartPointer<vtkDelaunay3D> del = vtkSmartPointer<vtkDelaunay3D>::New();
	del->SetInputData( inshape );
	del->BoundingTriangulationOff();
	del->Update();

	//delete cells with different pointids
	vtkCellArray* cells = del->GetOutput()->GetCells();
	vtkSmartPointer<vtkCellArray> newcells = vtkSmartPointer<vtkCellArray>::New();
	std::vector<int> celltypes;
	
	vtkShortArray *scalars = dynamic_cast<vtkShortArray*>(del->GetOutput()->GetPointData()->GetArray("regionID"));
	
	
	for( int i=0; i<del->GetOutput()->GetNumberOfCells(); i++ )
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		del->GetOutput()->GetCellPoints(i, ptIds);
						
		int id0 = ptIds->GetId(0);
		int id1 = ptIds->GetId(1);
		int id2 = ptIds->GetId(2);
		int id3 = ptIds->GetId(3);

		bool allequal = scalars->GetTuple1(id0) == scalars->GetTuple1(id1) && 
						scalars->GetTuple1(id0) == scalars->GetTuple1(id2) &&
						scalars->GetTuple1(id0) == scalars->GetTuple1(id3);
						
						
		if (allequal)
		{
			newcells->InsertNextCell( ptIds );
			celltypes.push_back( del->GetOutput()->GetCellType(i) );
		}
		
	}

	
	del->GetOutput()->SetCells( &celltypes[0], newcells );
	
	//extract LV surface
	del->GetOutput()->GetCellData()->SetActiveScalars("regionID");
	vtkSmartPointer<vtkPolyData> lv = vtkSmartPointer<vtkPolyData>::Take(
		ExtractSurface( del->GetOutput(), LV_ID ) );
	vtkSmartPointer<vtkPolyData> rv = vtkSmartPointer<vtkPolyData>::Take(
		ExtractSurface( del->GetOutput(), RV_ID ) );


	//generate epicardium
	inshape->GetPointData()->SetActiveScalars("regionID");
	vtkSmartPointer<vtkThreshold> thold = vtkSmartPointer<vtkThreshold>::New();
	thold->SetInputData(inshape);
	thold->ThresholdBetween( EPI_ID, EPI_ID );
	thold->AllScalarsOn();
	thold->Update();	
	
	vtkSmartPointer<vtkDelaunay3D> del_epi = vtkSmartPointer<vtkDelaunay3D>::New();
	del_epi->SetInputData( thold->GetOutput() );
	del_epi->BoundingTriangulationOff();
	del_epi->Update();

	vtkSmartPointer<vtkDataSetSurfaceFilter> surf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	surf->SetInputData( del_epi->GetOutput() );
	surf->Update();	
	
	vtkSmartPointer<vtkPolyData> epi = vtkSmartPointer<vtkPolyData>::New();
	epi->DeepCopy( surf->GetOutput() );
	
	
	//extract myocardium
	del->Update();

	//delete cells with different pointids
	cells = del->GetOutput()->GetCells();
	newcells = vtkSmartPointer<vtkCellArray>::New();
	celltypes.clear();
		
	for( int i=0; i<del->GetOutput()->GetNumberOfCells(); i++ )
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		del->GetOutput()->GetCellPoints(i, ptIds);
						
		int id0 = ptIds->GetId(0);
		int id1 = ptIds->GetId(1);
		int id2 = ptIds->GetId(2);
		int id3 = ptIds->GetId(3);

		bool allequal = scalars->GetTuple1(id0) == scalars->GetTuple1(id1) && 
						scalars->GetTuple1(id0) == scalars->GetTuple1(id2) &&
						scalars->GetTuple1(id0) == scalars->GetTuple1(id3);
						
		bool allepi = scalars->GetTuple1(id0) == EPI_ID &&
				scalars->GetTuple1(id1) == EPI_ID &&
				scalars->GetTuple1(id2) == EPI_ID &&
				scalars->GetTuple1(id3) == EPI_ID;
						
		if (!allequal || allepi)
		{
			newcells->InsertNextCell( ptIds );
			celltypes.push_back( del->GetOutput()->GetCellType(i) );
		}
		
	}
	
	del->GetOutput()->SetCells( &celltypes[0], newcells );
	surf->SetInputData( del->GetOutput() );
	surf->Update();		
	
	vtkSmartPointer<vtkPolyData> myo = vtkSmartPointer<vtkPolyData>::New();
	myo->DeepCopy( surf->GetOutput() );
	
	
	//CommonTools::SaveUnstructuredGrid(  del->GetOutput(), outshapefile);
	CommonTools::SaveShapeToFile( lv, "lv.vtk", NULL);
	CommonTools::SaveShapeToFile( rv, "rv.vtk", NULL);
	CommonTools::SaveShapeToFile( epi, "epi.vtk", NULL);
	CommonTools::SaveShapeToFile( myo, "myo.vtk", NULL);

	return 0;
}





vtkPolyData* ExtractSurface( vtkUnstructuredGrid *volmesh, int id )
{
	//volmesh->GetCellData()->SetActiveScalars("regionID");
	
	vtkSmartPointer<vtkThreshold> thold = vtkSmartPointer<vtkThreshold>::New();
	thold->SetInputData( volmesh );
	thold->ThresholdBetween( id, id );
	thold->AllScalarsOn();
	thold->Update();
	
	vtkSmartPointer<vtkDataSetSurfaceFilter> surf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	surf->SetInputData( thold->GetOutput() );
	surf->Update();
	
	vtkPolyData *result = vtkPolyData::New();
	result->DeepCopy(surf->GetOutput());
	
	return result;

}


/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Same as PassScalars, except that the search is done for every point/cell of target, not source!
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



void PassScalarsFloat( int argc, char *argv[], double radius )
{	

	std::cout << " reading source " << argv[1] << std::endl;
	vtkSmartPointer<vtkPolyData> sourcePd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[1]) );

	std::cout << " reading target " << argv[2] << std::endl;
	vtkSmartPointer<vtkPolyData> targetPd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[2]) );


	char* property_name = argv[3];
	std::cout << " property name: " << property_name << std::endl;

	bool use_point_data = sourcePd->GetPointData()->GetArray(property_name)!=NULL;

	if( use_point_data )
	{
		std::cout<<"Copying pointdata"<<std::endl;

		vtkFloatArray *target_scalars = (vtkFloatArray*) targetPd -> GetPointData() -> GetArray(property_name);
		vtkFloatArray *source_scalars = (vtkFloatArray*) sourcePd -> GetPointData() -> GetArray(property_name);

		vtkSmartPointer<vtkPointLocator> in_source_locator = vtkSmartPointer<vtkPointLocator>::New();
		in_source_locator->SetDataSet(sourcePd);
		in_source_locator->BuildLocator();

		for (int i=0; i<targetPd -> GetNumberOfPoints(); i++) 
		{
			double *p1 = targetPd -> GetPoint(i);
                        
                        vtkIdType pointId;
                        if(radius>=0)
                        {
                            double dist2;
                            pointId  = in_source_locator->FindClosestPointWithinRadius(radius, p1, dist2);
                        }
                        else
                            pointId  = in_source_locator->FindClosestPoint(p1);//(p1, p2, cellId, subId, dist);

                        
                        if(pointId>=0)
                            target_scalars->SetValue(i,source_scalars->GetValue(pointId));
                        else
                            target_scalars->SetValue(i,-1);
		}
	}
	else
	{
		std::cout<<"Copying celldata"<<std::endl;

		vtkFloatArray *target_scalars = (vtkFloatArray*) targetPd -> GetCellData() -> GetArray(property_name);
		vtkFloatArray *source_scalars = (vtkFloatArray*) sourcePd -> GetCellData() -> GetArray(property_name);

                std::cout<<"Building locator"<<std::endl;
		vtkSmartPointer<vtkCellLocator> in_source_locator = vtkSmartPointer<vtkCellLocator>::New();
		in_source_locator->SetDataSet(sourcePd);
		in_source_locator->BuildLocator();

		for( int i=0; i<targetPd->GetNumberOfCells(); i++)
		{
			double cell_center_p[3];
			double pt[3];
			double weights[3];
			int not_used;
			targetPd->GetCell(i)->GetParametricCenter(cell_center_p);
			targetPd->GetCell(i)->EvaluateLocation(not_used,cell_center_p,pt,weights);

			double closestpoint[3];
			vtkIdType cellid;
			int subid;
			double dist2;

                        vtkIdType pointId;
                        in_source_locator->FindClosestPoint(pt,closestpoint,cellid,subid, dist2);
                        
                        if(radius>=0)
                            if(dist2>radius*radius)
                                cellid = -1;

                        if(cellid>=0)
                            target_scalars->SetValue(i,source_scalars->GetValue(cellid));
                        else
                            target_scalars->SetValue(i,-1);		}
	}

	char outputFileName[255];
	strcpy(outputFileName,argv[4]);
	std::cout << " save output to " << outputFileName << std::endl;

	CommonTools::SaveShapeToFile(targetPd,outputFileName);

}



void PassScalarsShort( int argc, char *argv[], double radius )
{	


	std::cout << " reading source " << argv[1] << std::endl;
	vtkSmartPointer<vtkPolyData> sourcePd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[1]) );

	std::cout << " reading target " << argv[2] << std::endl;
	vtkSmartPointer<vtkPolyData> targetPd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(argv[2]) );


	char* property_name = argv[3];
	std::cout << " property name: " << property_name << std::endl;

	bool use_point_data = sourcePd->GetPointData()->GetArray(property_name)!=NULL;

	if( use_point_data )
	{
		std::cout<<"Copying pointdata"<<std::endl;

		vtkShortArray *target_scalars = (vtkShortArray*) targetPd -> GetPointData() -> GetArray(property_name);
		vtkShortArray *source_scalars = (vtkShortArray*) sourcePd -> GetPointData() -> GetArray(property_name);

		vtkSmartPointer<vtkPointLocator> in_source_locator = vtkSmartPointer<vtkPointLocator>::New();
		in_source_locator->SetDataSet(sourcePd);
		in_source_locator->BuildLocator();

		for (int i=0; i<targetPd -> GetNumberOfPoints(); i++) 
		{
			double *p1 = targetPd -> GetPoint(i);
                        
                        vtkIdType pointId;
                        if(radius>=0)
                        {
                            double dist2;
                            pointId  = in_source_locator->FindClosestPointWithinRadius(radius, p1, dist2);
                        }
                        else
                            pointId  = in_source_locator->FindClosestPoint(p1);//(p1, p2, cellId, subId, dist);

                        
                        if(pointId>=0)
                            target_scalars->SetValue(i,source_scalars->GetValue(pointId));
                        else
                            target_scalars->SetValue(i,-1);
		}
	}
	else
	{
		std::cout<<"Copying celldata"<<std::endl;

		vtkShortArray *target_scalars = (vtkShortArray*) targetPd -> GetCellData() -> GetArray(property_name);
		vtkShortArray *source_scalars = (vtkShortArray*) sourcePd -> GetCellData() -> GetArray(property_name);

		vtkSmartPointer<vtkCellLocator> in_source_locator = vtkSmartPointer<vtkCellLocator>::New();
		in_source_locator->SetDataSet(sourcePd);
		in_source_locator->BuildLocator();

		for( int i=0; i<targetPd->GetNumberOfCells(); i++)
		{
			double cell_center_p[3];
			double pt[3];
			double weights[3];
			int not_used;
			targetPd->GetCell(i)->GetParametricCenter(cell_center_p);
			targetPd->GetCell(i)->EvaluateLocation(not_used,cell_center_p,pt,weights);

			double closestpoint[3];
			vtkIdType cellid;
			int subid;
			double dist2;

                        vtkIdType pointId;
                        in_source_locator->FindClosestPoint(pt,closestpoint,cellid,subid, dist2);
                        
                        if(radius>=0)
                            if(dist2>radius*radius)
                                cellid = -1;

                        if(cellid>=0)
                            target_scalars->SetValue(i,source_scalars->GetValue(cellid));
                        else
                            target_scalars->SetValue(i,-1);
		}
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
		std::cout << "Same as PassScalars, except that the search is done for every point/cell of target, not source!" << std::endl;
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0] << " <source> <target> <property name> <output> <short|float> [radius]" << std::endl;
		return EXIT_FAILURE;
	}

	const char* datatype = argv[5];

        double radius = -1; //search radius
        if (argv[6]!=NULL)
            radius = atof(argv[6]);
        
	if( strcmp(datatype,"short")==0 )
		PassScalarsShort( argc, argv, radius );
	else if( strcmp(datatype,"float")==0 )
		PassScalarsFloat( argc, argv, radius );
	else
	{
		std::cout<<"Unknown data type: "<<datatype<<std::endl;
		exit(-1);
	}

	return 0;

}

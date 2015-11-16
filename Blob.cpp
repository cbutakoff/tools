/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file Blob.cpp  
    \brief Active surfaces.
  An active surfaces implementation that fits a 3D blob to a set of points. Under construction
*/
#include <vtkSmartPointer.h>
#include "CommonTools.h"
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkDoubleArray.h>
#include <vtkMassProperties.h>
#include <vtkCurvatures.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPoints.h>
#include <vtkPointLocator.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkType.h>

#include <vector>
#include <iostream>


//! A function that aproximates curvature.
void ApproximateCurvature( vtkPolyData* mesh, vtkDoubleArray* normals );
//! A function that gets point's neighbors.
void GetPointNeighbors( vtkPolyData* mesh, vtkIdType ptid, vtkIdList *ptIds);


int main(int argc, char** argv)
{
	//const char* inshapefilename = argv[1];
	//vtkSmartPointer<vtkPolyData> inshape = vtkSmartPointer<vtkPolyData>::Take(
	//	CommonTools::LoadShapeFromFile( inshapefilename ) );
		
	vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
	sphere->SetRadius(20);
	sphere->SetThetaResolution(15);
	sphere->SetPhiResolution(15);
	sphere->SetCenter(0,0,0);
	sphere->Update();
	//vtkSmartPointer<vtkDataSetSurfaceFilter> filter_extractsurface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	//filter_extractsurface->SetInputConnection( sphere->GetOutputPort() );
	//filter_extractsurface->Update();
		
	vtkSmartPointer<vtkPolyData> inshape = vtkSmartPointer<vtkPolyData>::New();
	inshape->DeepCopy( sphere->GetOutput() );
		
	vtkSmartPointer<vtkPolyDataNormals> filter_normals = vtkSmartPointer<vtkPolyDataNormals>::New();
	filter_normals->ComputePointNormalsOn();
	filter_normals->ComputeCellNormalsOff();
	filter_normals->SplittingOff();
	
	vtkSmartPointer<vtkCurvatures> filter_curvature = vtkSmartPointer<vtkCurvatures>::New();
	filter_curvature->SetCurvatureTypeToGaussian();
	
	vtkSmartPointer<vtkMassProperties> filter_volume = vtkSmartPointer<vtkMassProperties>::New();
	
	vtkSmartPointer<vtkPointLocator> locator = vtkSmartPointer<vtkPointLocator>::New();
	
	
	double dt = 0.001;
	double k_curvature = 0.5;
	double k_volume = 0.2;
	double k_pointattraction = 1;
	
	const double target_volume = 4*3.14*20*20*20/3; //volume of a 3cm ball
	
	//points for external forces, poles
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	pts->SetDataTypeToDouble();
	pts->SetNumberOfPoints(2);
	pts->SetPoint(0, 0,0,30);
	pts->SetPoint(1, 0,0,-30);
	
	
	bool converged = false;
	int iter;
	vtkSmartPointer<vtkPolyData> oldshape = vtkSmartPointer<vtkPolyData>::New();
	oldshape->DeepCopy(inshape);

	double volume = 0;
	double error  = 0;
	double old_error = 1000;
	const int max_iter=100000;

	for( iter = 0; iter<max_iter && !converged; iter++ )
	{
	
		
		//compute normals
		//std::cout<<"Compute normals"<<std::endl;
		filter_normals->SetInputData( inshape );
		filter_normals->Update();
		vtkDoubleArray* normals_points = static_cast<vtkDoubleArray*>(
			filter_normals->GetOutput()->GetPointData()->GetNormals() );
		
		//compute curvature
		//std::cout<<"Compute curvature"<<std::endl;
		//filter_curvature->SetInputConnection( filter_normals->GetOutputPort() );
		//filter_curvature->Update();
		ApproximateCurvature( inshape, normals_points );
		
		//vtkDoubleArray* curvature = static_cast<vtkDoubleArray*>(
		//	filter_curvature->GetOutput()->GetPointData()->GetArray("Gauss_Curvature"));
		vtkDoubleArray* curvature = static_cast<vtkDoubleArray*>(
			inshape->GetPointData()->GetArray("Curvature"));
		
		//compute volume
		//std::cout<<"Compute volume"<<std::endl;
		filter_volume->SetInputData( inshape );
		filter_volume->Update();
		volume = filter_volume->GetVolume();
		
		//displace vertices
		// -a* curvature* normal + b(VOL - volume)*normal
		//inshape->DeepCopy(filter_curvature->GetOutput());
		
		//std::cout<<"Displace vertices"<<std::endl;
		
		for( int i=0; i<inshape->GetNumberOfPoints(); i++)
		{
			double pt[3];
			inshape->GetPoint(i, pt);
			
			double n[3];
			normals_points->GetTuple(i, n);
			
			pt[0] += dt*(k_curvature*curvature->GetValue(i) + 
				k_volume*(target_volume-volume) )*n[0];
			pt[1] += dt*(k_curvature*curvature->GetValue(i) + 
				k_volume*(target_volume-volume) )*n[1];
			pt[2] += dt*(k_curvature*curvature->GetValue(i) + 
				k_volume*(target_volume-volume) )*n[2];
				
			inshape->GetPoints()->SetPoint(i, pt);
		}
		
		//add external force due to points
		locator->SetDataSet( inshape );
		locator->Update();
		
		for( int i=0; i<pts->GetNumberOfPoints(); i++)
		{
			int ptid = locator->FindClosestPoint( pts->GetPoint(i) );
			double pt_closest[3];
			double pt_original[3];
			inshape->GetPoint( ptid, pt_closest );
			pts->GetPoint( i, pt_original );
							
			double v[3];
			v[0] = pt_original[0] - pt_closest[0];
			v[1] = pt_original[1] - pt_closest[1];
			v[2] = pt_original[2] - pt_closest[2];

			double dist = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] ); 
			v[0] /= dist;
			v[1] /= dist;
			v[2] /= dist;

							
			pt_closest[0] += dt* k_pointattraction *v[0];
			pt_closest[1] += dt* k_pointattraction *v[1];
			pt_closest[2] += dt* k_pointattraction *v[2];
				
			inshape->GetPoints()->SetPoint(i, pt_closest);		
		}
		
		//inshape->Update();
		//check for conveence

		error  = 0;
		for( int i=0; i<inshape->GetNumberOfPoints(); i++)
		{
			double pt1[3];
			double pt2[3];
			inshape->GetPoint(i,pt1);
			oldshape->GetPoint(i,pt2);
			error += sqrtf((pt1[0]-pt2[0])*(pt1[0]-pt2[0]) +
					(pt1[0]-pt2[0])*(pt1[0]-pt2[0]) +
					(pt1[0]-pt2[0])*(pt1[0]-pt2[0]))/inshape->GetNumberOfPoints();
		}

		//if( fabs(old_error-error)<0.000001 ) converged=true;
		if( error<0.000001 ) converged=true;
		
		if( iter%1000 ==0 )
		{
			char filename[100];
			sprintf(filename, "shape%07d.vtk",iter);
			CommonTools::SaveShapeToFile( inshape, filename );
			printf("Iteration %07d/%07d. Volume: %0.2f/%0.2f. Error: %f\n",iter,max_iter, volume, target_volume, error);			
		}	

		
		oldshape->DeepCopy(inshape);	
		old_error=error;
		//printf("curv: %f \n",	curvature->GetValue(209));

		
	}
	
	
		
	char filename[100];
	sprintf(filename, "shape%07d.vtk",iter);
	CommonTools::SaveShapeToFile( inshape, filename );
	
    return 0;
}



void ApproximateCurvature( vtkPolyData* mesh, vtkDoubleArray* normals )
{
	//approximate curvature as length of the vector connecting the vertex with the center of mass of neighborhood
	const int npts = mesh->GetNumberOfPoints();
	
	vtkSmartPointer<vtkDoubleArray> curvature = vtkSmartPointer<vtkDoubleArray>::New();
	curvature->SetName("Curvature");
	curvature->SetNumberOfValues( npts );
	
	for( int i = 0; i<npts; i++ )
	{
		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		GetPointNeighbors( mesh, i, ptIds );
		
		double centr[] = {0,0,0};
		for( int j=0; j<ptIds->GetNumberOfIds(); j++)
		{
			double pt[3];
			mesh->GetPoint( ptIds->GetId(j), pt );
			centr[0] += pt[0];
			centr[1] += pt[1];
			centr[2] += pt[2];
		}
		centr[0] /= ptIds->GetNumberOfIds();
		centr[1] /= ptIds->GetNumberOfIds();
		centr[2] /= ptIds->GetNumberOfIds();
		
		double v[3];
		double vert[3];
		mesh->GetPoint( i, vert );
		v[0] = centr[0]-vert[0];
		v[1] = centr[1]-vert[1];
		v[2] = centr[2]-vert[2];
		
		double n[3];
		normals->GetTuple( i, n );
		const double dp = v[0]*n[0] + v[1]*n[1] + v[2]*n[2];
		
		//const double length = sqrtf( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
		curvature->SetValue( i, dp );
	}
	
	mesh->GetPointData()->RemoveArray( "Curvature" );
	mesh->GetPointData()->AddArray( curvature );
}


void GetPointNeighbors( vtkPolyData* mesh, vtkIdType ptid, vtkIdList *ptIds)
{
	ptIds->Reset();
	
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells( ptid,  cellIds );
	
	for(int j=0; j<cellIds->GetNumberOfIds(); j++ )
	{
		vtkSmartPointer<vtkIdList> ptIds_local = vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints( cellIds->GetId(j), ptIds_local );
		
		ptIds->InsertUniqueId( ptIds_local->GetId(0) );
		ptIds->InsertUniqueId( ptIds_local->GetId(1) );
		ptIds->InsertUniqueId( ptIds_local->GetId(2) );
	}
	
}

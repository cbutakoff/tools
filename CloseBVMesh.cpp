/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file CloseBVMesh.cpp  
    \brief Closes biventricular mesh by connecting endocardial edge to epicardial at the base. Does not close ventricles.

  It was made for Rafa's biventricular model. The mesh must have epicardium, rv endo and lv endo separable. No scalars are necessary. 
  The scalars must be vtkShortArray!!! (type short)
*/
#include "CommonTools.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkDelaunay2D.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkShortArray.h>
#include <vtkCellData.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkCellLocator.h>
#include <vtkPolygon.h>
#include <vtkStripper.h>

//! Copies cell scalars from source mesh to the target 
void CopyCellScalars(vtkPolyData* src, vtkPolyData* tgt, char* scalars_name, int fill_value);

//! Convert vtkPolyData polygon to vtkPolygon 
vtkPolygon* PolyData2Polygon( vtkPolyData* pd );

//! Verify if a point is in vtkPolygon 
int PointInPolygon( double x[3], vtkPolygon* pg );

//! Calculate area of a vtkPolygon 
double PolygonBoundaryArea( vtkPolygon* polygon );


//! Fills small holes in vtkPolyData by connecting vertices to the centroid
vtkPolyData* FillSmallHoles( vtkPolyData* pd );


//------------------------------------------------------------------
int main( int argc, char *argv[] )
{	

	if( argc < 2 )
	{
		std::cout << "Closes biventricular mesh v1.1" << std::endl;
		std::cout << "Usage: CloseBVMesh inshape outshape [cell_scalars]" << std::endl;
		std::cout << "cell_scalars - array of cell scalars to copy (e.g. simulation_scalars). Leave blank for no scalars" << std::endl;
		return -1;
	}

//	const char scalars_name[] = "simulation_scalars";
//	const char scalars_name[] = "RegionId";
	char *scalars_name = NULL;


	char* inshape = argv[1];
	char* outshape = argv[2];

	if( argc==4 )
	{
		scalars_name = argv[3];
	}
	


	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(inshape) );

	{
		std::cout<<"Clean up the mesh...";			
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner->SetInputData(pd);
//		cleaner->SetTolerance(0.0001);
//		cleaner->PointMergingOn();
		cleaner->Update();
		pd->DeepCopy( cleaner->GetOutput() );
		std::cout<<"done"<<std::endl;			
	}	
		
	std::cout<<"Connectivity filter...";	
	vtkSmartPointer<vtkPolyDataConnectivityFilter> split = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	split->SetInputData(pd);
	split->SetExtractionModeToAllRegions();
	split->ColorRegionsOn();
	split->Update();
	std::cout<<"done"<<std::endl;	
	
	//feature edges, all off , leave boudary only
	std::cout<<"Feature edges filter...";	
	vtkSmartPointer<vtkFeatureEdges> fe=vtkSmartPointer<vtkFeatureEdges>::New();
	fe->SetInputData( split->GetOutput() );
	fe->BoundaryEdgesOn ();
	fe->FeatureEdgesOff ();
	fe->NonManifoldEdgesOff ();
	fe->ManifoldEdgesOff ();
	fe->Update();
	std::cout<<"done"<<std::endl;	
	
	//delaunay 2d
	std::cout<<"Delaunay 2D filter...";	
	vtkSmartPointer<vtkDelaunay2D> del = vtkSmartPointer<vtkDelaunay2D>::New();
	del->SetInputData( fe->GetOutput() );
	del->SetTolerance(0);
	del->Update();
	std::cout<<"done"<<std::endl;	
		
	//CommonTools::SaveShapeToFile( del->GetOutput(), "delaunay.vtk");
		
		
	//clean triangles
	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
	std::cout<<"Removing extra triangles...";	

	//prepare the 3 polygons
/*	std::vector< vtkSmartPointer<vtkPolygon>  > polygons(3);
	double areas[3];
	for( int i=0; i<3; i++ )
	{
		vtkSmartPointer<vtkPolyData> sub = vtkSmartPointer<vtkPolyData>::Take(
			CommonTools::GetShapeSubSurface(split->GetOutput(), i, i) );
		
		polygons[i] = vtkSmartPointer<vtkPolygon>::Take(
			PolyData2Polygon( sub ) ); 
			
		areas[i] = PolygonBoundaryArea( polygons[i] );
	}

	int largest_poly_index = areas[1]>areas[0]?(areas[2]>areas[1]?2:1):(areas[0]>areas[2]?0:2);
	*/


	for( int i=0; i<del->GetOutput()->GetNumberOfCells(); i++)
	{
		vtkCell *cell = del->GetOutput()->GetCell(i);

		vtkIdType pts3[3];
		pts3[0] = cell->GetPointId(0);
		pts3[1] = cell->GetPointId(1);
		pts3[2] = cell->GetPointId(2);
		
		int s1, s2, s3;
		s1 = del->GetOutput()->GetPointData()->GetScalars()->GetTuple1(pts3[0]);
		s2 = del->GetOutput()->GetPointData()->GetScalars()->GetTuple1(pts3[1]);
		s3 = del->GetOutput()->GetPointData()->GetScalars()->GetTuple1(pts3[2]);
		
		//std::cout<<"Cell "<<i<<" scalar "<<del->GetOutput()->GetPointData()->GetScalars()->GetName()
		//	<<" : "<<s1<<" "<<s2<<" "<<s3<<std::endl;
		
		if( !(s1==s2 && s1==s3 && s2==s3) )			
			polys->InsertNextCell(3,pts3);	
	/*	else 
		{ 
			//candidate for elimination. 
			//get the cell centroid
			double pcoords[3];
			double x[3];
			double weights[3];
			int subId = cell->GetParametricCenter(pcoords);
			cell->EvaluateLocation(subId, pcoords, x, weights);
			
			//the centroid must be outside the smaller polys and inside the big poly
			int p[3];
			p[0] = PointInPolygon( x, polygons[0] );
			p[1] = PointInPolygon( x, polygons[1] );
			p[2] = PointInPolygon( x, polygons[2] );
			
			int the_condition = p[largest_poly_index];
			for( int j=0; j<3; j++ )
			{
				if( j!=largest_poly_index )
				{ 
					the_condition += (p[j]==0?1:0);
				}
			}
			
			if( the_condition==3 )
			{
				polys->InsertNextCell(3,pts3);	
			}
		}*/
	}

	vtkSmartPointer<vtkPolyData> cover = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->DeepCopy(del->GetOutput()->GetPoints());
	cover->SetPoints(points);
	cover->SetPolys(polys);
	//cover->Update();
	std::cout<<"done"<<std::endl;	


/*
	vtkSmartPointer<vtkShortArray> coverscalars = vtkSmartPointer<vtkShortArray>::New(); 
	coverscalars->SetName(scalars_name);
	coverscalars->SetNumberOfValues( cover->GetNumberOfCells() );
	for( int i=0; i<cover->GetNumberOfCells(); i++ )
		coverscalars->SetValue(i,100);

	cover->GetCellData()->AddArray( coverscalars );		
	cover->GetCellData()->CopyScalarsOn();
	cover->Update();
*/


	//CommonTools::SaveShapeToFile( cover, "cover.vtk");
	

	//final merging of the labelled stuff
	std::cout<<"Merge the surfaces with the cover...";				
	vtkSmartPointer<vtkAppendPolyData> final = vtkSmartPointer<vtkAppendPolyData>::New();
	final->AddInputData (pd);
	final->AddInputData (cover);
	final->Update();
	std::cout<<"done"<<std::endl;			
		

	//CommonTools::SaveShapeToFile( final->GetOutput(), "merged.vtk");
			
		
	std::cout<<"Clean up the mesh...";			
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(final->GetOutput());
	cleaner->Update();
	std::cout<<"done"<<std::endl;			

	vtkSmartPointer<vtkPolyData> fixedmesh = vtkSmartPointer<vtkPolyData>::Take(
		FillSmallHoles( cleaner->GetOutput() ) );

//	vtkSmartPointer<vtkPolyData> fixedmesh = vtkSmartPointer<vtkPolyData>::New();
//	fixedmesh->DeepCopy( cleaner->GetOutput() );


	if( scalars_name!=NULL )
	{
		CopyCellScalars(pd, fixedmesh, scalars_name, 100);
		fixedmesh->GetCellData()->SetActiveScalars(scalars_name);		
	}
		
	//run tesselate

	CommonTools::SaveShapeToFile( fixedmesh, outshape );
	
	return 0;
}


void CopyCellScalars(vtkPolyData* src, vtkPolyData* tgt, char* scalars_name, int fill_value)
{
	std::cout<<"Building locator...";
	vtkSmartPointer<vtkCellLocator> locator = vtkSmartPointer<vtkCellLocator>::New();
	locator->SetDataSet( src );
	locator->BuildLocator();
	locator->Update();
	std::cout<<"done"<<std::endl;


	std::cout<<"Reserving scalars...";	
	vtkSmartPointer<vtkShortArray> scalars = vtkSmartPointer<vtkShortArray>::New();
	scalars->SetName(scalars_name);
	scalars->SetNumberOfValues(tgt->GetNumberOfCells());
	 
	if( src->GetCellData()->GetArray(scalars_name)==NULL )
	{
		std::cout<<"Cell array "<<scalars_name<<" not found in the input mesh"<<std::endl;
	}


	for( int i=0; i<tgt->GetNumberOfCells(); i++ )
		scalars->SetValue(i, fill_value);
	std::cout<<"done"<<std::endl;	

	std::cout<<"Copying scalars...";		
	for( int i=0; i<tgt->GetNumberOfCells(); i++ )
	{
		//get the center of the cell and find the closest cell on the original meshes

		double pcoords[3];
		int subId = tgt->GetCell(i)->GetParametricCenter(pcoords);
		
		double x[3];
		double weights[3];
		tgt->GetCell(i)->EvaluateLocation(subId, pcoords, x, weights);
		
		double closestpoint[3];
		vtkIdType cellId;

		double dist2;
		locator->FindClosestPoint(x, closestpoint, cellId, subId, dist2 );

		int value = static_cast<vtkShortArray*>(src->GetCellData()->GetArray(scalars_name))->GetValue(cellId);

		if( dist2<1e-2)
		{
			scalars->SetValue(i, value);
		}		
	}
		
	tgt->GetCellData()->AddArray(scalars);	
}


vtkPolygon* PolyData2Polygon( vtkPolyData* pd )
{
	vtkSmartPointer<vtkStripper> strip = vtkSmartPointer<vtkStripper>::New();
	strip->SetInputData(pd);
	strip->Update();
	
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(strip->GetOutput());
	cleaner->Update();
	
	vtkPolygon* pg = vtkPolygon::New();
	for( int i = 0; i<cleaner->GetOutput()->GetNumberOfPoints(); i++ )
	{
		double x[3];
		cleaner->GetOutput()->GetPoint(i,x);
		x[2] = 0.0;
		pg->GetPoints()->InsertNextPoint(x );
	}
		
	return pg;
	
}

int PointInPolygon( double x[3], vtkPolygon* pg )
{
	double n[3];
	pg->ComputeNormal(pg->GetPoints()->GetNumberOfPoints(),
          static_cast<double*>(pg->GetPoints()->GetData()->GetVoidPointer(0)), n);	
	
	double bounds[6];
	pg->GetPoints()->GetBounds(bounds);
   
	return pg->PointInPolygon(x,
		  pg->GetPoints()->GetNumberOfPoints(), static_cast<double*>(
		  pg->GetPoints()->GetData()->GetVoidPointer(0)), bounds, n);
}

double PolygonBoundaryArea( vtkPolygon* polygon )
{
	double bounds[6];
	polygon->GetPoints()->GetBounds(bounds);
	
	const double a=bounds[1]-bounds[0];
	const double b=bounds[3]-bounds[2];
	const double c=bounds[5]-bounds[4];
	
	double result = (a?a:1)*(b?b:1)*(c?c:1);
	return result;
}

vtkPolyData* FillSmallHoles( vtkPolyData* pd )
{
	//feature edges, all off , leave boudary only
	vtkSmartPointer<vtkFeatureEdges> fe=vtkSmartPointer<vtkFeatureEdges>::New();
	fe->SetInputData( pd );
	fe->BoundaryEdgesOn ();
	fe->FeatureEdgesOff ();
	fe->NonManifoldEdgesOff ();
	fe->ManifoldEdgesOff ();
	fe->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> split = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	split->SetInputData(fe->GetOutput());
	split->SetExtractionModeToAllRegions();
	split->ColorRegionsOn();
	split->Update();
	
	vtkSmartPointer<vtkAppendPolyData> final = vtkSmartPointer<vtkAppendPolyData>::New();
	final->AddInputData (pd);

	for( int i=0; i<split->GetNumberOfExtractedRegions(); i++)
	{
		vtkSmartPointer<vtkPolyData> sub = vtkSmartPointer<vtkPolyData>::Take(
			CommonTools::GetShapeSubSurface(split->GetOutput(), i, i) );

//		char filename[100];
//		sprintf(filename,"sub %03d.vtk",i);
//		CommonTools::SaveShapeToFile( sub, filename, NULL);

		std::cout<<"Hole id :"<<i<<" pts: "<<sub->GetNumberOfPoints()<<std::endl;

/*
		vtkSmartPointer<vtkStripper> strip = vtkSmartPointer<vtkStripper>::New();
		strip->SetInputData(sub);
		strip->Update();
		
		vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner->SetInputData(strip->GetOutput());
		cleaner->Update();
	*/
		vtkSmartPointer<vtkPolyData> cover = vtkSmartPointer<vtkPolyData>::Take(
			CommonTools::GenerateHoleCover( sub ) );
		

		final->AddInputData (cover);

		
//		sprintf(filename,"delaunay_%03d.vtk",i);
//		CommonTools::SaveShapeToFile( cover, filename, NULL);
		
	}

	final->Update();

	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(final->GetOutput());
	cleaner->Update();


	vtkPolyData* result = vtkPolyData::New();
	result->DeepCopy( cleaner->GetOutput() );
	return result;
}


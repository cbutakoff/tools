/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Generate labels for Rafa's biventricular mesh
*/
#include <vtkPolyData.h>

#include <vtkDataSet.h>


#include <vtkAppendPolyData.h>

#include <vtkSmartPointer.h>

#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkPolyDataNormals.h"

#include "vtkPolyDataConnectivityFilter.h"
#include "vtkShortArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkCell.h"
#include "vtkPointData.h"

#include "CommonTools.h"

#include "vnl/vnl_vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//! add scalar array to the mesh. Don't remember why it is not in common tools
void SetScalars(vtkPolyData* mesh, int value, const char* array);


int main(int argc, char **argv)
{
	std::cout << "Version 2.0"<< std::endl;

	if( argc<2 ) 
	{
		std::cout << "Params: "<< std::endl;
		std::cout << "-i <shape.vtk> \t\t- intput shape"<< std::endl;
		std::cout << "-rvepi <shape.vtk> \t\t- shape"<< std::endl;
		std::cout << "-lvmyo <shape.vtk> \t\t- shape"<< std::endl;
		std::cout << "-o <shape.vtk> \t\t- output shape"<< std::endl;
		return -1;
	}

	char *inputshapefile=NULL;
	char *rvepishapefile=NULL;
	char *lvmyoshapefile=NULL;
	char *outputshapefile=NULL;

	const int SCALARS_LV_ENDO = 1;
	const int SCALARS_RV_ENDO = 2;
	const int SCALARS_LV_EPI = 3;
	const int SCALARS_RV_EPI = 4;
	const char scalars_name[] = "simulation_scalars";

	
	
	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-i")==0 )
		{
			inputshapefile = argv[++c];
		}
		else if( strcmp(argv[c],"-rvepi")==0 )
		{
			rvepishapefile = argv[++c];
		}
		else if( strcmp(argv[c],"-lvmyo")==0 )
		{
			lvmyoshapefile = argv[++c];
		}
		else if( strcmp(argv[c],"-o")==0 )
		{
			outputshapefile = argv[++c];
		}

	}


	std::cout<<"Labeling options:"<<std::endl;		
	std::cout<<"Array: "<<scalars_name<<std::endl;		
	std::cout<<"LV endo scalar: "<<SCALARS_LV_ENDO<<std::endl;		
	std::cout<<"RV endo scalar: "<<SCALARS_RV_ENDO<<std::endl;		
	std::cout<<"LV epi scalar: "<<SCALARS_LV_EPI<<std::endl;		
	std::cout<<"RV epi scalar: "<<SCALARS_RV_EPI<<std::endl<<std::endl;		
	

	std::cout<<"Loading shapes: "<<inputshapefile<<" | "<<rvepishapefile<<" | "<<lvmyoshapefile<<"...";		
	vtkSmartPointer<vtkPolyData> inshape = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(inputshapefile));
	vtkSmartPointer<vtkPolyData> rvepi = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(rvepishapefile));
	vtkSmartPointer<vtkPolyData> lvmyo = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(lvmyoshapefile));
	std::cout<<"done"<<std::endl;	
		
		
	//generate scalars for both inputs
	std::cout<<"Reset scalars for rvepi and lvmyo...";	
	SetScalars(rvepi, SCALARS_RV_EPI, scalars_name);
	SetScalars(lvmyo, SCALARS_LV_EPI, scalars_name);
	std::cout<<"done"<<std::endl;	


	//join the two meshes
	std::cout<<"Append rvepi and lvmyo...";	
	vtkSmartPointer<vtkAppendPolyData> original_labeled = vtkSmartPointer<vtkAppendPolyData>::New();
	original_labeled->AddInputData(rvepi);
	original_labeled->AddInputData(lvmyo);
	original_labeled->Update();
	std::cout<<"done"<<std::endl;	
	
	
	//split the input mesh
	std::cout<<"Connectivity filter...";	
	vtkSmartPointer<vtkPolyDataConnectivityFilter> split = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	split->SetInputData(inshape);
	split->SetExtractionModeToAllRegions();
	split->ColorRegionsOn();
	split->Update();
	std::cout<<"done"<<std::endl;	
	
	
	std::vector< vtkSmartPointer<vtkPolyData>  > sub(3);
	
	std::cout<<"Save connectivity output...";	
	for( int i=0; i<3; i++ )
		sub[i] = vtkSmartPointer<vtkPolyData>::Take(
			CommonTools::GetShapeSubSurface(split->GetOutput(), i, i) );
	std::cout<<"done"<<std::endl;	


	std::cout<<"Searching for epicardium...";	

	//find epi, we will have to split it into 2 parts
	int epi_index = 0; //largest number of points
	
	for( int i=1; i<3; i++ )
		if( sub[i]->GetNumberOfPoints() > sub[epi_index]->GetNumberOfPoints() )
			epi_index = i;
			
	std::cout<<"done"<<std::endl;	

			
	std::cout<<"Building locator for epicardium...";
	vtkSmartPointer<vtkCellLocator> locator = vtkSmartPointer<vtkCellLocator>::New();
	locator->SetDataSet( original_labeled->GetOutput() );
	locator->BuildLocator();
	locator->Update();
	std::cout<<"done"<<std::endl;


	std::cout<<"Reserving epi scalars...";	
	vtkSmartPointer<vtkShortArray> scalars = vtkSmartPointer<vtkShortArray>::New();
	scalars->SetName(scalars_name);
	scalars->SetNumberOfValues(sub[epi_index]->GetNumberOfCells());
	for( int i=0; i<sub[epi_index]->GetNumberOfCells(); i++ )
		scalars->SetValue(i, SCALARS_LV_EPI);
	std::cout<<"done"<<std::endl;	

	std::cout<<"Separating epi scalars into lv and rv...";		
	for( int i=0; i<sub[epi_index]->GetNumberOfCells(); i++ )
	{
		//get the center of the cell and find the closest cell on the original meshes

		double pcoords[3];
		int subId = sub[epi_index]->GetCell(i)->GetParametricCenter(pcoords);
		//std::cout<<"Pcenter "<<pcoords[0]<<" "<<pcoords[1]<<" "<<pcoords[2]<<std::endl;


		
		double x[3];
		double weights[3];
		sub[epi_index]->GetCell(i)->EvaluateLocation(subId, pcoords, x, weights);
		//std::cout<<"X "<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;

		
		double closestpoint[3];
		vtkIdType cellId;

		//std::cout<<"Getting closest point"<<std::endl;

		double dist2;
		locator->FindClosestPoint(x, closestpoint, cellId, subId, dist2 );

		//std::cout<<"Found scalar at cell "<<cellId<<std::endl;
		int value = static_cast<vtkShortArray*>(original_labeled->GetOutput()->GetCellData()->GetArray(scalars_name))->GetValue(cellId);
		//std::cout<<"with value "<<value<<std::endl;
		scalars->SetValue(i, value);
		
	}
		
	sub[epi_index]->GetCellData()->AddArray(scalars);
	std::cout<<"done"<<std::endl;	
	
	//////////////////////////////////////////////////
	//
	// Produce the endocardia
	//
	//
	//

	std::cout<<"Generate lvendo and rvendo scalars...";		
	//compute the centroids
	vnl_vector<double> rvcenter(3);
	vnl_vector<double> lvcenter(3);
	rvepi->GetCenter(rvcenter.begin());
	lvmyo->GetCenter(lvcenter.begin());
	
	int rvindex=-1;
	for(int i=0; i<3; i++)
	{
		vnl_vector<double> subcenter(3);
		
		if( i!=epi_index )
		{
			sub[i]->GetCenter( subcenter.begin() );
			double d_rv = (subcenter - rvcenter).squared_magnitude();
			double d_lv = (subcenter - lvcenter).squared_magnitude();
			
			if( d_rv > d_lv )
				SetScalars( sub[i], SCALARS_LV_ENDO, scalars_name);
			else
			{
				SetScalars( sub[i], SCALARS_RV_ENDO, scalars_name);				
				rvindex=i;
			}
				
		}
	
	}
	std::cout<<"done"<<std::endl;		

	
	std::cout<<"Build locator for lvmyo...";			
	//tweak rvendo to have also lvepi regionids
	locator->SetDataSet( lvmyo );
	locator->BuildLocator();
	locator->Update();
	std::cout<<"done"<<std::endl;			

	const double mark_thold = 0.5;
	std::cout<<"Adding lvepi to rvendo. Everything at dist2<"<<mark_thold<<" from lvmyo is marked...";			
	vtkShortArray* s = static_cast<vtkShortArray*>(sub[rvindex]->GetCellData()->GetArray(scalars_name));
	
	for( int i=0; i<sub[rvindex]->GetNumberOfCells(); i++ )
	{
		//get the center of the cell and find the closest cell on the original meshes

		double pcoords[3];
		int subId = sub[rvindex]->GetCell(i)->GetParametricCenter(pcoords);
		//std::cout<<"Pcenter "<<pcoords[0]<<" "<<pcoords[1]<<" "<<pcoords[2]<<std::endl;
		
		double x[3];
		double weights[3];
		sub[rvindex]->GetCell(i)->EvaluateLocation(subId, pcoords, x, weights);
		//std::cout<<"X "<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;

		
		double closestpoint[3];
		vtkIdType cellId;

		//std::cout<<"Getting closest point"<<std::endl;

		double dist2;
		locator->FindClosestPoint(x, closestpoint, cellId, subId, dist2 );

		if( dist2<mark_thold )
		{
			//std::cout<<"Fixing rvendo scalars"<<cellId<<std::endl;
			s->SetValue(i, SCALARS_LV_EPI);
		}
	}	
	//sub[rvindex]->Update();
	std::cout<<"done"<<std::endl;			

	
	std::cout<<"Merge the surfaces with new scalars...";			
	
	//final merging of the labelled stuff
	vtkSmartPointer<vtkAppendPolyData> final = vtkSmartPointer<vtkAppendPolyData>::New();
	final->AddInputData(sub[0]);
	final->AddInputData(sub[1]);
	final->AddInputData(sub[2]);
	final->Update();
	std::cout<<"done"<<std::endl;			

	
	std::cout<<"Saving "<<outputshapefile<<"...";			
	
	final->GetOutput()->GetPointData()->RemoveArray("RegionId");
	final->GetOutput()->GetCellData()->SetActiveScalars(scalars_name);
	
	CommonTools::SaveShapeToFile(final->GetOutput(), outputshapefile );
	std::cout<<"done"<<std::endl;			

	return 0;
}



void SetScalars(vtkPolyData* mesh, int value, const char* array)
{
	const int ncells = mesh->GetNumberOfCells();
	vtkSmartPointer<vtkShortArray> ar = vtkSmartPointer<vtkShortArray>::New();
	ar->SetNumberOfValues(ncells);
	ar->SetName(array);
	 
	for(int i=0; i<ncells; i++) ar->SetValue(i, value);
	
	mesh->GetCellData()->AddArray(ar);
}


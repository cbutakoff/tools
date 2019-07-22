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
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkCell.h>
#include <vtkExtractEdges.h>
#include "VTKCommonTools.h"


#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <limits>

#ifdef _OPENMP
# include <omp.h>
#endif


using namespace std;






//------------------------------------------------------------------
int main( int argc, char *argv[] )
{	

	if( argc <2  )
	{
		cout<<"Saves max edge length per cell into a text file"<<endl;
		cout<<"Supports openmp"<<endl<<endl;
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0] << "SetScalars <mesh.vtk> <outputfile.txt>" << std::endl;
		return EXIT_FAILURE;
	}


	int c=1;
	const char *filenamein = argv[c++];
	const char *filenameout = argv[c++];


	int ncpus = 1;
	#ifdef _OPENMP
	ncpus = omp_get_max_threads();
	cout<<"OpenMP enabled, using "<<ncpus<<" threads"<<endl<<endl;
	#endif
	

	cout << " opening " << filenamein << endl;
	vtkSmartPointer<vtkDataSetReader> rd = vtkSmartPointer<vtkDataSetReader>::New();
	rd->SetFileName(filenamein);
	CommonTools::AssociateProgressFunction(rd);
  	rd->Update();

	auto mesh = rd->GetOutput();

	vector<double> edgelength_global;
	edgelength_global.reserve(rd->GetOutput()->GetNumberOfCells());


	vtkIdType steps_completed = 0;
	const vtkIdType step_size   = 100;
	const vtkIdType total_steps = rd->GetOutput()->GetNumberOfCells() / step_size + 1;

	#pragma parallel 
	{
	vector<double> edgelength_local;
	edgelength_local.reserve( round( rd->GetOutput()->GetNumberOfCells()/ncpus ) );

	vtkIdType local_steps_count = 0;


	#pragma omp for nowait schedule(auto) 
	for( vtkIdType i=0; i<rd->GetOutput()->GetNumberOfCells(); i++ )
	{
            	if (local_steps_count++ % step_size == step_size-1)
            	{
			#pragma omp atomic
                	++steps_completed;

	                if (steps_completed % 100 == 1)
                	{
				#pragma omp critical
				cout<<"Done "<<steps_completed<<"/"<<total_steps<<"; "<< steps_completed*100/total_steps<< "%\r"<<flush;
                	}
            	}

		auto cell = mesh->GetCell(i);
		if(  cell->GetNumberOfPoints()!=4 )
		{
			cout<<endl<<"Only tetras supported. Cell "<< i <<" has "<<cell->GetNumberOfPoints()<<" vertices"<<endl;
			throw;
		}

		double coordinates[4][3];
		
		for(int kkk=0; kkk<4; kkk++)
			mesh->GetPoint( cell->GetPointId(kkk), coordinates[kkk] );
		

		int edge_ids[6][2] =
			    {
				{1,0},
				{2,1},
				{0,2},
				{3,0},
				{3,1},
				{3,2}
			    };

		double edge[3];
		double max_dist = 0;

		for( int kkk=0; kkk<6; kkk++)
		{
			for(int xyz=0; xyz<3; xyz++)
			{
				edge[xyz] = coordinates[ edge_ids[kkk][0] ][xyz] - coordinates[ edge_ids[kkk][1] ][xyz];
			}
			const double d = sqrt( edge[0]*edge[0] + edge[1]*edge[1] + edge[2]*edge[2] );

			if( d > max_dist ) max_dist = d;
		}

		edgelength_local.push_back(max_dist);

  	}
	//end parallel for
  	cout<<endl;
	

	#pragma omp critical
	{
	#ifdef _OPENMP
	//cout<<edgelength_local[0]<<endl;
	copy(edgelength_local.begin(), edgelength_local.end(), back_inserter(edgelength_global));
	#else
	edgelength_global.assign( edgelength_local.begin(), edgelength_local.end() );
	#endif
	}// end of critical

	}//end parallel

	cout<<"Saving the edge lengths"<<endl;


	ofstream ff( filenameout );
	for (vtkIdType k=0; k<edgelength_global.size(); ++k)
	{
		ff<<edgelength_global[k]<<endl;
	}


	return EXIT_SUCCESS;
}

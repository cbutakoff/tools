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
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkCell.h>
#include <vtkExtractEdges.h>
#include <vtkNew.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include "VTKCommonTools.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

//------------------------------------------------------------------
int main(int argc, char *argv[])
{

	if (argc < 2)
	{
		cout << "Saves max edge length per cell into a text file" << endl;
		cout << "Supports openmp" << endl
				 << endl;
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0] << "SetScalars <mesh.vtk|vtu> <outputfile.vtu>" << std::endl;
		return EXIT_FAILURE;
	}

	int c = 1;
	const char *filenamein = argv[c++];
	const char *filenameout = argv[c++];

	int ncpus = 1;
#ifdef _OPENMP
	ncpus = omp_get_max_threads();
	cout << "OpenMP enabled, using " << ncpus << " threads" << endl
			 << endl;
#endif

	vtkSmartPointer<vtkUnstructuredGrid> mesh;

	cout << " opening " << filenamein << endl;
	if (filenamein[strlen(filenamein) - 1] == 'u' || filenamein[strlen(filenamein) - 1] == 'U')
	{ //.vtu
		vtkNew<vtkXMLUnstructuredGridReader> rd;
		rd->SetFileName(filenamein);
		CommonTools::AssociateProgressFunction(rd);
		rd->Update();
		mesh = vtkUnstructuredGrid::SafeDownCast(rd->GetOutput());
	}
	else
	{
		vtkNew<vtkDataSetReader> rd;
		rd->SetFileName(filenamein);
		CommonTools::AssociateProgressFunction(rd);
		rd->Update();
		mesh = vtkUnstructuredGrid::SafeDownCast(rd->GetOutput());
	}

	vtkNew<vtkFloatArray> els_max;
	els_max->SetNumberOfComponents(1);
	els_max->SetNumberOfTuples(mesh->GetNumberOfCells());
	els_max->SetName("EdgeLengthMax");
	vtkNew<vtkFloatArray> els_min;
	els_min->SetNumberOfComponents(1);
	els_min->SetNumberOfTuples(mesh->GetNumberOfCells());
	els_min->SetName("EdgeLengthMin");



	vtkIdType steps_completed = 0;
	const vtkIdType step_size = 100;
	const vtkIdType total_steps = mesh->GetNumberOfCells() / step_size + 1;

#pragma parallel
	{
		vector<double> edgelength_max_local;
		vector<double> edgelength_min_local;
		vector<vtkIdType> cellIds_local;
		edgelength_max_local.reserve(round(mesh->GetNumberOfCells() / ncpus));
		edgelength_min_local.reserve(round(mesh->GetNumberOfCells() / ncpus));
		cellIds_local       .reserve(round(mesh->GetNumberOfCells() / ncpus));

		vtkIdType local_steps_count = 0;

#pragma omp for nowait schedule(auto)
		for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++)
		{
			if (local_steps_count++ % step_size == step_size - 1)
			{
#pragma omp atomic
				++steps_completed;

				if (steps_completed % 100 == 1)
				{
#pragma omp critical
					cout << "Done " << steps_completed << "/" << total_steps << "; " << steps_completed * 100 / total_steps << "%\r" << flush;
				}
			}

			auto cell = mesh->GetCell(i);
		  vtkIdType nEdges = cell->GetNumberOfEdges();

			double max_dist = 0;
			double min_dist = std::numeric_limits<double>::max();

			for( vtkIdType edgeid=0; edgeid<nEdges; edgeid++){
				vtkCell* edge = cell->GetEdge(edgeid);
				const vtkIdType pt0id = edge->GetPointId(0);
				const vtkIdType pt1id = edge->GetPointId(1);
				double pt0[3], pt1[3];
				mesh->GetPoint(pt0id, pt0);
				mesh->GetPoint(pt1id, pt1);
				const double d = sqrt( (pt1[0]-pt0[0])*(pt1[0]-pt0[0]) + 
      				 								 (pt1[1]-pt0[1])*(pt1[1]-pt0[1]) + 
      												 (pt1[2]-pt0[2])*(pt1[2]-pt0[2]) );
      
				if (d > max_dist)
					max_dist = d;
				if (d < min_dist)
					min_dist = d;

			}

			edgelength_max_local.push_back(max_dist);
			edgelength_min_local.push_back(min_dist);
			cellIds_local.push_back(i);
		}
		// end parallel for
		cout << endl;

#pragma omp critical
		{
			for( size_t i=0; i<cellIds_local.size(); i++ )
			{
    		els_max->SetTuple1( cellIds_local[i], edgelength_max_local[i] );
    		els_min->SetTuple1( cellIds_local[i], edgelength_min_local[i] );
			}
		} // end of critical

	} // end parallel

	cout << "Saving the edge lengths mesh" << endl;

	mesh->GetCellData()->AddArray(els_max);
	mesh->GetCellData()->AddArray(els_min);

	vtkNew<vtkXMLUnstructuredGridWriter> wr;
	wr->SetFileName(filenameout);
	wr->SetInputData(mesh);
	wr->EncodeAppendedDataOff();
	wr->SetCompressorTypeToNone();
	wr->SetHeaderTypeToUInt64();
	wr->Write();


	return EXIT_SUCCESS;
}

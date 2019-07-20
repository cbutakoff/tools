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

#include "VTKCommonTools.h"


#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>



using namespace std;






//------------------------------------------------------------------
int main( int argc, char *argv[] )
{	

	if( argc <2  )
	{
		std::cout << "Usage: " << std::endl;
		std::cout << argv[0] << "SetScalars <mesh.vtk> <outputfile.txt>" << std::endl;
		return EXIT_FAILURE;
	}


	int c=1;
	const char *filenamein = argv[c++];
	const char *filenameout = argv[c++];


	cout << " opening " << filenamein << endl;
	vtkSmartPointer<vtkDataSetReader> rd = vtkSmartPointer<vtkDataSetReader>::New();
	rd->SetFileName(filenamein);
  	rd->Update();

  	cout<<"Calculating edge lengths"<<endl;

	vtkIdType npoints = rd->GetOutput()->GetNumberOfPoints();	
	//SpMat edges( npoints, npoints );


	ofstream ff( filenameout );
	for( vtkIdType i=0; i<rd->GetOutput()->GetNumberOfCells(); i++ )
	{
      		if( i%10000 == 0 )
			cout<<"Done "<<i<<"/"<<rd->GetOutput()->GetNumberOfCells()<<"; "<< i*100/rd->GetOutput()->GetNumberOfCells()<< "%\r"<<flush;

		auto cell = rd->GetOutput()->GetCell(i);
		//cout<<"Cell "<<i<<": ";
		//for(int kk=0; kk<cell->GetNumberOfPoints();kk++)		
		//	cout<<cell->GetPointId(kk)<<",";
		//cout<<endl;


      		for( int j=0; j<cell->GetNumberOfEdges(); j++ )
      		{
          		auto edge = cell->GetEdge(j);
          		auto p1 = edge->GetPointId(0);
          		auto p2 = edge->GetPointId(1);
			//cout<<"edge "<<p1<<" "<<p2<<endl;
	
			//if( Get( edges, p1, p2)==0 )
			{
				double coords1[3];
				double coords2[3];
			  	rd->GetOutput()->GetPoint(p1, coords1);
				rd->GetOutput()->GetPoint(p2, coords2);	
				double d = sqrt((coords1[0]-coords2[0])*(coords1[0]-coords2[0]) +
						(coords1[1]-coords2[1])*(coords1[1]-coords2[1]) +
						(coords1[2]-coords2[2])*(coords1[2]-coords2[2]));
				//Insert( edges, p1, p2, d );
				ff<<d<<endl;
			}

      		}
  	}
  	cout<<endl;
/*
	cout<<"Saving the edge lengths"<<endl;


	for (vtkIdType k=0; k<edges.outerSize(); ++k)
		for (SpMat::InnerIterator it(edges,k); it; ++it)
		{
			ff<<it.value()<<endl;
		}

*/
	return EXIT_SUCCESS;
}

/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
	\brief Convert chaste .ele .node tetrahedral mesh to VTK unstructured grid (only for Oxford Rabbit)
	*/
	
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkType.h>
#include <vtkTetra.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <vnl/vnl_vector.h>

#include "CommonTools.h"

//! Reads binary chaste .node file and stores into vtkPoints
void ReadChasteNodes(const char* filename, vtkPoints* pts );
//! Reads binary chaste .ele file of tetrahedra and stores into vtkCells, it also retrieves the first scalar array
void ReadChasteElements(const char* filename, vtkCellArray* cells, vtkIntArray* scalars);



int main( int argc, char **argv)
{
	std::cout<<"Convert chaste .ele .node tetrahedral mesh to VTK unstructured grid. v1.0"<<std::endl;
	std::cout<<"Options: <file.ele> <file.node> <output.vtk>"<<std::endl;
	
	if( argc<4 ) 
	{
		std:cout<<"Insufficient parameters"<<std::endl;
		return EXIT_FAILURE;
	}
	
	int c=1;
	const char *ele_filename = argv[c++];
	const char *node_filename = argv[c++];
	const char *out_filename = argv[c++];
	
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkIntArray> cellid = vtkSmartPointer<vtkIntArray>::New();
	
	ReadChasteNodes( node_filename, pts );

	//CommonTools::SavePoints(pts, "points.vtk");

	ReadChasteElements( ele_filename, cells, cellid);
	
	/*
	pts->SetNumberOfPoints(4);
	//cells->SetNumberOfCells( 1 );
	cellid->SetNumberOfValues( 1 );

	pts->SetPoint(0, 0,0,0 );
	pts->SetPoint(1, 1,0,0 );
	pts->SetPoint(2, 0,1,0 );
	pts->SetPoint(3, 0,0,1 );
	
	vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
	tetra->GetPointIds()->SetId( 0, 3 );
	tetra->GetPointIds()->SetId( 1, 2 );
	tetra->GetPointIds()->SetId( 2, 1 );
	tetra->GetPointIds()->SetId( 3, 0 );

	
	cellid->SetValue( 0, 1 );
	
	cells->InitTraversal();
	cells->InsertNextCell(tetra);
	*/

	
	vtkSmartPointer<vtkUnstructuredGrid> output = vtkSmartPointer<vtkUnstructuredGrid>::New();
	output->SetPoints(pts);
	output->SetCells(VTK_TETRA,cells);
	output->GetCellData()->AddArray(cellid);
	//output->Update();

	std::cout<<"Ncells: "<<output->GetNumberOfCells()<<std::endl;

	vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetFileName(out_filename);
	writer->SetInputData(output);
	writer->SetFileTypeToBinary();
	writer->Update();
	
	
	return EXIT_SUCCESS;
}


//! Reads binary chaste .node file and stores into vtkPoints
void ReadChasteNodes(const char* filename, vtkPoints* pts )
{
	std::cout<<"1"<<std::endl;

	//first read the number of nodes
	const long int N = 4283195;
	const int dim = 3;
	const long int offset = 0x12;


	//then read the nodes
	FILE * pFile;
	pFile = fopen (filename,"rb");
	if (pFile==NULL)
	{
		std::cout<<"Open failed: "<<filename<<std::endl;
		throw("Open failed");
	}	

	vnl_vector<double> ptr(N*dim);
	
	fseek( pFile, offset, SEEK_SET );

	/*
	unsigned char abc[4];
	fread( abc, 1, 4, pFile );
	char buffer [33];
	float abcf;
	
    *((unsigned char*)(&abcf) + 3) = abc[0];
    *((unsigned char*)(&abcf) + 2) = abc[1];
    *((unsigned char*)(&abcf) + 1) = abc[2];
    *((unsigned char*)(&abcf) + 0) = abc[3];
	
	sprintf(buffer, "%x %x %x %x %f", abc[0], abc[1], abc[2], abc[3], abcf );
	std::cout<<"abc = "<<buffer<<std::endl;
	
	fseek( pFile, offset, SEEK_SET );
	*/
	
	std::cout<<"Attempting to read nodes"<<std::endl;
	fread( (double*)ptr.begin(), sizeof(double), N*dim, pFile );
	fclose(pFile);

	std::cout<<"Nodes loaded"<<std::endl;
	
	std::cout<<"1st node:"<<ptr[0]<<"; "<<ptr[1]<<"; "<<ptr[2]<<std::endl;
	std::cout<<"2nd node:"<<ptr[3]<<"; "<<ptr[4]<<"; "<<ptr[5]<<std::endl;
	
	pts->SetNumberOfPoints( N );
	
	for( int i=0; i<N; i++ )
	{
		const long int j = i*dim;
		pts->SetPoint(i, ptr[j], ptr[j+1], ptr[j+2] );
	}
	
}

//! Reads binary chaste .ele file of tetrahedra and stores into vtkCells, it also retrieves the first scalar array
void ReadChasteElements(const char* filename, vtkCellArray* cells, vtkIntArray* scalars)
{
	//first read the number of nodes
	const long int N = 24099051;
	const int dim = 5; //tetra,  4+1 for scalars
	const long int offset = 0x11;

	//then read the nodes
	FILE * pFile;
	pFile = fopen (filename,"rb");
	if (pFile==NULL)
	{
		std::cout<<"Open failed: "<<filename<<std::endl;
		throw("Open failed");
	}	

	vnl_vector<vtkTypeUInt32> ptr(N*dim);

	
	fseek( pFile, offset, SEEK_SET );
	fread( (vtkTypeUInt32*)ptr.begin(), sizeof(vtkTypeUInt32), N*dim, pFile );
	fclose(pFile);
	
	//cells->SetNumberOfCells( N );
	scalars->SetNumberOfValues( N );

	std::cout<<"1st element:"<<ptr[0]<<"; "<<ptr[1]<<"; "<<ptr[2]<<"; "<<ptr[3]<<"; "<<ptr[4]<<std::endl;
	std::cout<<"2nd element:"<<ptr[5]<<"; "<<ptr[6]<<"; "<<ptr[7]<<"; "<<ptr[8]<<"; "<<ptr[9]<<std::endl;

	cells->Reset();
	cells->InitTraversal();


	for( int i=0; i<N; i++ )
	{
		if( i%10000 == 0  )
		{
			printf("Elements added: %010d / %d \r", i, N );
		}
		
		const long int k = i*dim;
		vtkSmartPointer<vtkTetra> tetra = vtkSmartPointer<vtkTetra>::New();
		for(int j=0; j<4; j++)
		{
			tetra->GetPointIds()->SetId( j, ptr[k+j] );
		}
		scalars->SetValue( i, ptr[k+4] );
		cells->InsertNextCell(tetra);
	}
	std::cout<<std::endl;
	
}


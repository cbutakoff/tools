/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Make biventricular mesh for Rafa's model
    */

#include <vtkImageData.h>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include "vtkImageOpenClose3D.h"
#include <vtkImageContinuousDilate3D.h>

#include "CommonTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>




int main(int argc, char **argv)
{
	std::cout << "Version 1.1"<< std::endl;

	if( argc<2 ) 
	{
		std::cout << "Params: "<< std::endl;
		std::cout << "-i <image.vtk> \t\t- input"<< std::endl;
		std::cout << "-o <image.vtk> \t\t- output"<< std::endl;
		std::cout << "-vo <number> \t\t- open value"<< std::endl;
		std::cout << "-vc <number> \t\t- close value"<< std::endl;
		std::cout << "-ks <x> <y> <z> \t\t- kernel size, 3-vector"<< std::endl;
		return EXIT_FAILURE;
	}

	char *inputimagefile=NULL;
	char *outputimagefile=NULL;
	float openavalue = 0;
	float closevalue = 255;
	int kx = 3;
	int ky = 3;
	int kz = 3;

	
	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-i")==0 )
		{
			inputimagefile = argv[++c];
		}
		else if( strcmp(argv[c],"-o")==0 )
		{
			outputimagefile = argv[++c];
		}
		else if( strcmp(argv[c],"-vo")==0 ) //open value
		{
			openavalue = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-vc")==0 ) //close value
		{
			closevalue = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-ks")==0 ) //kernel size (3-vector)
		{
			kx = atoi(argv[++c]);
			ky = atoi(argv[++c]);
			kz = atoi(argv[++c]);
		}

	}


	vtkSmartPointer<vtkImageData> inimage = vtkSmartPointer<vtkImageData>::Take(
		CommonTools::LoadImage( inputimagefile ) );

	std::cout<<"Applying vtkOpenClose3D. Open v.: "<<openavalue<<"; close v.: "<<closevalue<<"; kernel: "<<kx<<", "<<ky<<", "<<kz<<std::endl;
	
	/*
	vtkSmartPointer<vtkImageOpenClose3D> openclose = vtkSmartPointer<vtkImageOpenClose3D>::New();
	openclose->SetInputData( inimage );
	openclose->SetKernelSize(kx,ky,kz);
	openclose->SetOpenValue(openavalue);
	openclose->SetCloseValue(closevalue);
	openclose->Update();
	*/
	
	vtkSmartPointer<vtkImageContinuousDilate3D> dilate = vtkSmartPointer<vtkImageContinuousDilate3D>::New();
	dilate->SetKernelSize( kx,ky,kz );
	dilate->SetInputData( inimage );
	dilate->Update();
		
	CommonTools::SaveImage( dilate->GetOutput(), outputimagefile );
		
	return EXIT_SUCCESS;
}

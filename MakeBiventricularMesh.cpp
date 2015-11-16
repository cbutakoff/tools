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
#include <vtkPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>

#include <vtkImageData.h>
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>

#include <vtkAppendPolyData.h>

#include <vtkSmartPointer.h>

#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkImageWeightedSum.h"
#include "vtkImageOpenClose3D.h"
#include "vtkImageMarchingCubes.h"
#include "vtkImageContinuousErode3D.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkSmoothPolyDataFilter.h"

#include "CommonTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>


//! creates an empty image
vtkImageData* CreateEmptyImage( double* bounds, double* spacing, double padding, float value );
//! creates a mask
vtkImageData* CreateMask( vtkImageData* image, vtkPolyData* shape, float value );
//! decimates mesh
vtkPolyData* DecimateMesh( vtkPolyData* mesh, int nFaces );


int main(int argc, char **argv)
{
	std::cout << "Version 1.1"<< std::endl;

	if( argc<2 ) 
	{
		std::cout << "Params: "<< std::endl;
		std::cout << "-lvendo <shape.vtk> \t\t- shape"<< std::endl;
		std::cout << "-rvepi <shape.vtk> \t\t- shape"<< std::endl;
		std::cout << "-lvmyo <shape.vtk> \t\t- shape"<< std::endl;
		std::cout << "-o <shape.vtk> \t\t- output shape"<< std::endl;
		return -1;
	}

	char *lvendoshapefile=NULL;
	char *rvepishapefile=NULL;
	char *lvmyoshapefile=NULL;
	char *outputshapefile=NULL;

	double vx=0.5,vy=0.5,vz=0.5; //voxel dimensions
	//double vx=1,vy=1,vz=1; //voxel dimensions
	
	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-lvendo")==0 )
		{
			lvendoshapefile = argv[++c];
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

	vtkSmartPointer<vtkPolyData> lvendo = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(lvendoshapefile));
	vtkSmartPointer<vtkPolyData> rvepi = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(rvepishapefile));
	vtkSmartPointer<vtkPolyData> lvmyo = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(lvmyoshapefile));
		
		
		

	double spacing[3];

	const int padding = 10; //space to add to shape extremes


	

	

	std::cout<<"Generating an empty image"<<std::endl;

	//combine all  the meshes to get the dimensions
	vtkSmartPointer<vtkAppendPolyData> app = vtkSmartPointer<vtkAppendPolyData>::New();
	app->AddInputData(lvendo);
	app->AddInputData(rvepi);
	app->AddInputData(lvmyo);
	app->Update();

	spacing[0]=vx;
	spacing[1]=vy;
	spacing[2]=vz;

	double bounds[6];
	app->GetOutput()->GetBounds(bounds);
	std::cout<<"Bounds: "<<bounds[0]<<", "<<bounds[1]<<", "<<bounds[2]<<", "<<bounds[3]<<", "<<bounds[4]<<", "<<bounds[5]<<std::endl;

	vtkSmartPointer<vtkImageData> image_final = vtkSmartPointer<vtkImageData>::Take( CreateEmptyImage( bounds, spacing, padding, 0 ) );
	CommonTools::SaveImage( image_final, outputshapefile );
	
	std::cout<<"Creating lvendo image"<<std::endl;
	vtkSmartPointer<vtkImageData> image_lvendo = vtkSmartPointer<vtkImageData>::Take( CreateMask( image_final, lvendo, 10 ) );
	std::cout<<"Creating lvmyo image"<<std::endl;
	vtkSmartPointer<vtkImageData> image_lvmyo = vtkSmartPointer<vtkImageData>::Take( CreateMask( image_final, lvmyo, 10 ) );
	std::cout<<"Creating rvepi image"<<std::endl;
	vtkSmartPointer<vtkImageData> image_rvepi = vtkSmartPointer<vtkImageData>::Take( CreateMask( image_final, rvepi, 10 ) );

	std::cout<<"Summing the masks"<<std::endl;
	vtkSmartPointer<vtkImageWeightedSum> sum = vtkSmartPointer<vtkImageWeightedSum>::New();
	sum->AddInputData( image_lvendo );
	sum->AddInputData( image_lvmyo );
	sum->AddInputData( image_rvepi );
	sum->SetWeight( 0, 1 );
	sum->SetWeight( 1, 1 );
	sum->SetWeight( 2, 1 );
	sum->NormalizeByWeightOff();
	sum->Update();
	
	std::cout<<"Applying morphological operations to generate RV endo"<<std::endl;
	vtkSmartPointer<vtkImageOpenClose3D> openclose = vtkSmartPointer<vtkImageOpenClose3D>::New();
	openclose->SetInputData( sum->GetOutput() );
	openclose->SetKernelSize(3,3,3);
	openclose->SetOpenValue(0);
	openclose->SetCloseValue(10);
	openclose->Update();
	
	vtkSmartPointer<vtkImageContinuousErode3D> erode = vtkSmartPointer<vtkImageContinuousErode3D>::New();
	erode->SetInputData( openclose->GetOutput() );
	erode->SetKernelSize(7,7,7);
	erode->Update();

	std::cout<<"Restoring the left ventricle"<<std::endl;

	vtkSmartPointer<vtkImageWeightedSum> sum1 = vtkSmartPointer<vtkImageWeightedSum>::New();
	sum1->AddInputData( openclose->GetOutput() );
	sum1->AddInputData( erode->GetOutput() );
	sum1->AddInputData( image_lvendo );
	sum1->AddInputData( image_lvmyo );
	sum1->SetWeight( 0, 1 );
	sum1->SetWeight( 1, -2 );
	sum1->SetWeight( 2, -4 );
	sum1->SetWeight( 3, 4 );
	sum1->NormalizeByWeightOff();
	sum1->Update();

	std::cout<<"Smooth the image"<<std::endl;
	vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
	smooth->SetInputData( sum1->GetOutput() );
	smooth->SetStandardDeviation( 1.5 );	
	smooth->Update();

	std::cout<<"Marching cubes"<<std::endl;

	vtkSmartPointer<vtkImageMarchingCubes> mc = vtkSmartPointer<vtkImageMarchingCubes>::New();
	mc->SetInputData( smooth->GetOutput() );
	mc->SetValue( 0, 1 );
	mc->Update();

	std::cout<<"Smooth mesh"<<std::endl;
	vtkSmartPointer<vtkSmoothPolyDataFilter> meshsmooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	meshsmooth->SetInputData( mc->GetOutput() );
	meshsmooth->FeatureEdgeSmoothingOn();
	meshsmooth->SetNumberOfIterations(1000); 
	meshsmooth->Update();

	vtkSmartPointer<vtkPolyData> finalmesh = vtkSmartPointer<vtkPolyData>::Take( DecimateMesh( meshsmooth->GetOutput(), 10000 ) );

	std::cout<<"Saving "<<outputshapefile<<std::endl;
	CommonTools::SavePolydata( finalmesh, outputshapefile, true );
//	CommonTools::SaveImage( openclose->GetOutput(), outputshapefile );

	return 0;
}


vtkPolyData* DecimateMesh( vtkPolyData* mesh, int nFaces )
{
		const char temp_filename[] = "decimation_temporary_file12894286.ply";
		CommonTools::SaveShapeToFile(mesh, temp_filename);
		CommonTools::GenerateDecimationScript("decimation_script_196239462.mlx",nFaces);
		char cmdline[1000];
		sprintf(cmdline,"meshlabserver -i %s -o %s -s decimation_script_196239462.mlx", temp_filename, temp_filename);
		int rubbish = system(cmdline);
		vtkPolyData *tempPD = CommonTools::LoadShapeFromFile( temp_filename );
		remove("decimation_script_196239462.mlx");
		remove(temp_filename);
		
		return tempPD;
}


vtkImageData* CreateMask( vtkImageData* image, vtkPolyData* shape, float value )
{

	vtkSmartPointer<vtkPolyDataToImageStencil> filt = 
		vtkSmartPointer<vtkPolyDataToImageStencil>::New();


	filt->SetOutputOrigin( image->GetOrigin() );
	filt->SetOutputSpacing( image->GetSpacing() );
	filt->SetOutputWholeExtent( image->GetExtent() );
	filt->SetInputData( shape );
	filt->Update();

	vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
	stencil->SetStencilData(filt->GetOutput());
	stencil->SetBackgroundValue(value);
	stencil->SetInputData(image);
	stencil->ReverseStencilOn();
	
	vtkImageData* imageout = vtkImageData::New();
	//stencil->GetOutput()->Update();
	imageout->DeepCopy( dynamic_cast<vtkImageData*>( stencil->GetOutput()) );
	//imageout->Update();
	
	
	return imageout;

}



vtkImageData* CreateEmptyImage( double* bounds, double* spacing, double padding, float value )
{
		vtkImageData* image = vtkImageData::New();
		double origin[3];
	
		double t[3];
		t[0] = -bounds[0]+padding*spacing[0];
		t[1] = -bounds[2]+padding*spacing[1];
		t[2] = -bounds[4]+padding*spacing[2];


		origin[0]=-t[0];
		origin[1]=-t[1];
		origin[2]=-t[2];
		int dims[3];

		dims[0]=floor((bounds[1]-bounds[0])/spacing[0])+padding*2;
		dims[1]=floor((bounds[3]-bounds[2])/spacing[1])+padding*2;
		dims[2]=floor((bounds[5]-bounds[4])/spacing[2])+padding*2;

		//image->SetScalarTypeToFloat();
		image->SetOrigin(origin);
		image->SetSpacing(spacing);
		image->SetDimensions(dims);
		//image->Update();
		image->AllocateScalars(VTK_FLOAT,1);
		memset(image->GetScalarPointer(),value,dims[0]*dims[1]*dims[2]*sizeof(float));	
		
		return image;
}

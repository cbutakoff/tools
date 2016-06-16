/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Resample image, smooth, reconstruct shape.
    */
#include <vtkImageData.h>

#include <vtkSmartPointer.h>
#include <vtkImageResample.h>
#include <vtkImageMarchingCubes.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkCleanPolyData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkImageCast.h>
#include <vtkImageConstantPad.h>
#include <vtkImageChangeInformation.h>
#include <vtkShortArray.h>
#include <vtkCellData.h>

#include "CommonTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>


int main(int argc, char **argv)
{
	std::cout << "Resample image, smooth, reconstruct shape. Version 2.0"<< std::endl;

	if( argc<2 ) 
	{
		std::cout << "Params: "<< std::endl;
		std::cout << "-image <image.vtk> \t\t- image (mask)"<< std::endl;
		std::cout << "-out_shape <shape.vtk> \t\t- resulting shape"<< std::endl;
		std::cout << "-sigma <float> \t\t- Gaussian smoothing. Default=3."<< std::endl;
		std::cout << "-region <int> \t\t- regionID to assign (scalars, unimportant, default=0)"<< std::endl;
		std::cout << "-faces <int> \t\t- number of faces to leave after simplification (def=500)"<< std::endl;
		std::cout << "-flip \t\t- flip the normals"<< std::endl;
		std::cout << "-fx <float> \t\t- resampling factor x (def=1)"<< std::endl;
		std::cout << "-fy <float> \t\t- resampling factor y (def=1)"<< std::endl;
		std::cout << "-fz <float> \t\t- resampling factor z (def=6, don't remember why)"<< std::endl;
		std::cout << "-vx <float> \t\t- spacing for x (def=0.2)"<< std::endl;
		std::cout << "-vy <float> \t\t- spacing for y (def=0.2)"<< std::endl;
		std::cout << "-vz <float> \t\t- spacing for z (def=0.2)"<< std::endl;
		std::cout << "-ri <image.vtk> \t\t- optional, to save resampled image"<< std::endl;		
		return -1;
	}

	char *inimagefile=NULL;
	char *outshapefile=NULL;
	char *resampledimagefile=NULL;
		
	float factorx = 1.0;
	float factory = 1.0;
	float factorz = 6.0;
	float vx = 0.2;
	float vy = 0.2;
	float vz = 0.2;
	float GaussSigma = 3.0;

	short int regionID = 0;
	int nFaces = 500;

	bool flipnormals=false;

	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-image")==0 )
		{
			inimagefile = argv[++c];
		}
		else if( strcmp(argv[c],"-out_shape")==0 )
		{
			outshapefile = argv[++c];
		}
		else if( strcmp(argv[c],"-ri")==0 )
		{
			resampledimagefile = argv[++c];
		}
		else if( strcmp(argv[c],"-sigma")==0 )
		{
			GaussSigma = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-region")==0 )
		{
			regionID = atoi(argv[++c]);
		}
		else if( strcmp(argv[c],"-faces")==0 )
		{
			nFaces = atoi(argv[++c]);
		}
		else if( strcmp(argv[c],"-flip")==0 )
		{
			flipnormals = true;
		}
		else if( strcmp(argv[c],"-fx")==0 )
		{
			factorx = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-fy")==0 )
		{
			factory = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-fz")==0 )
		{
			factorz = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-vx")==0 )
		{
			vx = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-vy")==0 )
		{
			vy = atof(argv[++c]);
		}
		else if( strcmp(argv[c],"-vz")==0 )
		{
			vz = atof(argv[++c]);
		}

		
	}

	

	vtkSmartPointer<vtkImageData> imagein = vtkSmartPointer<vtkImageData>::Take(
		CommonTools::LoadImage( inimagefile ) );


        if( outshapefile==NULL && resampledimagefile!=NULL )
        {
            //just resample
            std::cout<<"Resampling"<<std::endl;
            std::cout<<"Cast to float"<<std::endl;
            vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
            cast->SetInputData(imagein);
            cast->SetOutputScalarTypeToFloat();
            cast->Update();


            vtkSmartPointer<vtkImageResample> resamp = vtkSmartPointer<vtkImageResample>::New();
            resamp->SetInputData( cast->GetOutput() );
            resamp->SetAxisMagnificationFactor(0,factorx);
            resamp->SetAxisMagnificationFactor(1,factory);
            resamp->SetAxisMagnificationFactor(2,factorz);
            resamp->SetAxisOutputSpacing(0, vx);
            resamp->SetAxisOutputSpacing(1, vy);
            resamp->SetAxisOutputSpacing(2, vz);
            resamp->SetInterpolationModeToCubic();
            resamp->Update();

            CommonTools::SaveImage(resamp->GetOutput(), resampledimagefile);
            return 0;
        }
        
        
        
        
	std::cout<<"Cast to float"<<std::endl;
	vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
	cast->SetInputData(imagein);
	cast->SetOutputScalarTypeToFloat();
	cast->Update();

	std::cout<<"Pad along Z axis"<<std::endl;

	double spacing[3];
	cast->GetOutput()->GetSpacing(spacing);
	
	int padx=0;
	int pady=0;
	int padz=2;

	//trick to pad on the top and left
	vtkSmartPointer<vtkImageChangeInformation> translator = vtkSmartPointer<vtkImageChangeInformation>::New();
	translator->SetInputData(cast->GetOutput());
	translator->SetExtentTranslation(padx, pady, padz);
	translator->SetOriginTranslation(-padx*spacing[0], -pady*spacing[1], -padz*spacing[2]);
	translator->Update();

	int extent[6];
	translator->GetOutput()->GetExtent(extent);
	extent[4] = 0;
	//extent[5] += padz*spacing[2];
	
	vtkSmartPointer<vtkImageConstantPad> pad = vtkSmartPointer<vtkImageConstantPad>::New();
	pad->SetInputData(translator->GetOutput());
	pad->SetConstant(0);
	pad->SetOutputWholeExtent(extent);
	pad->Update();


	std::cout<<"Resampling"<<std::endl;
	vtkSmartPointer<vtkImageResample> resamp = vtkSmartPointer<vtkImageResample>::New();
	resamp->SetInputData( pad->GetOutput() );
	resamp->SetAxisMagnificationFactor(0,factorx);
	resamp->SetAxisMagnificationFactor(1,factory);
	resamp->SetAxisMagnificationFactor(2,factorz);
	resamp->SetAxisOutputSpacing(0, vx);
	resamp->SetAxisOutputSpacing(1, vy);
	resamp->SetAxisOutputSpacing(2, vz);
	resamp->SetInterpolationModeToCubic();
	resamp->Update();
	
	if( resampledimagefile!=NULL )
		CommonTools::SaveImage(resamp->GetOutput(), resampledimagefile);
	
	if( outshapefile!=NULL )
	{
		std::cout<<"Smoothing: "<<GaussSigma<<std::endl;
	
		vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
		smooth->SetStandardDeviation(GaussSigma);
		smooth->SetInputData(resamp->GetOutput());
		smooth->Update();

		//CommonTools::SaveImage(smooth->GetOutput(), "smoothed.vtk");

		std::cout<<"Running marching cubes"<<std::endl;

		vtkSmartPointer<vtkImageMarchingCubes> marchingcubes = vtkSmartPointer<vtkImageMarchingCubes>::New();
		marchingcubes->SetInputData(smooth->GetOutput());
		marchingcubes->ComputeNormalsOn();
		marchingcubes->SetValue(0,0.5);
		marchingcubes->Update();	


		const char *temp_filename = "cb_jn32ghkjdhq3hdkj.ply";

		vtkSmartPointer<vtkPolyData> tempPD;
		CommonTools::SaveShapeToFile(marchingcubes->GetOutput(),temp_filename);
		CommonTools::GenerateDecimationScript("decimation_script_196239462.mlx",nFaces);
		char cmdline[1000];
		sprintf(cmdline,"meshlabserver -i %s -o %s -s decimation_script_196239462.mlx", temp_filename, temp_filename);
		int rubbish = system(cmdline);
		tempPD = vtkSmartPointer< vtkPolyData > ::Take( CommonTools::LoadShapeFromFile( temp_filename) );
		remove("decimation_script_196239462.mlx");
		remove(temp_filename);

		//cleanup
		std::cout<<"Cleaning the mesh"<<std::endl;

		vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
		clean->SetInputData(tempPD);
		clean->Update();


		//generate normals
		if( flipnormals ) std::cout<<"Flipping normals"<<std::endl;
		vtkSmartPointer<vtkPolyDataNormals> normalgen = vtkSmartPointer<vtkPolyDataNormals>::New();
		normalgen->SetInputData(clean->GetOutput());
		normalgen->SetFlipNormals( flipnormals );
		normalgen->SplittingOff();
		normalgen->Update();

		vtkSmartPointer<vtkDataSetSurfaceFilter> extrsurf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
		extrsurf->SetInputData(normalgen->GetOutput());
		extrsurf->Update();

		//generate scalars
		vtkSmartPointer<vtkShortArray> regions = vtkSmartPointer<vtkShortArray>::New();
		regions->SetNumberOfValues(extrsurf->GetOutput()->GetNumberOfCells());

		for(int i=0; i<extrsurf->GetOutput()->GetNumberOfCells(); i++)
			regions->SetValue(i, regionID);

		regions->SetName("regionID");
		extrsurf->GetOutput()->GetCellData()->AddArray(regions);

		CommonTools::SaveShapeToFile(extrsurf->GetOutput(),outshapefile);
	}

	return 0;
}


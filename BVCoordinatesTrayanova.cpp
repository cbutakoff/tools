/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Generate layers of Left atrium. 
*/
#include <vtkCleanPolyData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkIVWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include "CommonTools.h"
#include "vtkSmartPointer.h"
#include <vector>
#include "vtkCellLocator.h"
#include "vtkPolyDataNormals.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_cross.h"
#include "vtkGenericCell.h"
#include "vtkWedge.h"
#include "vtkIdList.h"
#include "vtkCell.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkCellType.h"
#include "vtkCellData.h"
#include "vtkImageEuclideanDistance.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkImageGradient.h"
#include "vtkImageContinuousErode3D.h"
#include "vtkImageContinuousDilate3D.h"
#include "vtkShortArray.h"
#include "vtkCellArray.h"
#include "vtkStreamTracer.h"
#include "vtkImageMathematics.h"
#include "vtkImageCast.h"
#include "vtkAssignAttribute.h"
#include "vtkSplineFilter.h"
#include "vtkImageContinuousDilate3D.h"
#include "vtkImageContinuousErode3D.h"
#include "vtkImageResample.h"
#include "vtkImageMarchingCubes.h"
#include "vtkImageClip.h"
#include "vtkMath.h"
#include "vtkImageMathematics.h"
#include "vtkImageCast.h"
#include "vtkImageGaussianSmooth.h"
#include "vtkImageMarchingCubes.h"
#include "vtkCleanPolyData.h"
#include "vtkImageNormalize.h"
#include "vtkImageThreshold.h"
#include "vtkPointLocator.h"
#include <vtkvmtkPolyDataSurfaceRemeshing.h>

//#define SAVE_DEBUG_INFO


//! generate mask given epi and endo
void GenerateImageMask(vtkImageData* res_image, float fg, float bg, vtkPolyData* endo, 
			vtkPolyData* epi, float VoxelSize=0.5);
void GenerateTemperatureField( vtkPolyData* epi, vtkPolyData* lvendo, vtkPolyData* rvendo, 
        float VoxelSize/*=0.5*/, int laplace_iterations/*=100*/, 
        char* outfile_lv, char *outfile_rv, char* outfile_epi );


 
void usage(char *exe)
{
	std::cout<< "Generate scalar fields as in Trayanova's paper" << std::endl;
	std::cout<<std::endl;
	std::cout<< "Options:" << std::endl;
	std::cout<<"-rvendo <mesh.vtk> \t - endocardium mask (0-bg)"<<std::endl;
	std::cout<<"-lvendo <mesh.vtk> \t - epicardium mask (0-bg)"<<std::endl;
	std::cout<<"-epi <mesh.vtk> \t - epicardium mask (0-bg)"<<std::endl;
	std::cout<<"-voxel_lap <float> \t - Voxel size for laplace eq. solution (isotropic), the smaller-the better (e.g. 0.2)"<<std::endl;
	std::cout<<"-iter <int> \t - Number of iterations for Laplace field, default 100"<<std::endl;
	std::cout<<"-o1 <image.vtk> \t - Number of iterations for Laplace field, default 100"<<std::endl;
	std::cout<<"-o2 <image.vtk> \t - Number of iterations for Laplace field, default 100"<<std::endl;
	std::cout<<"-o3 <image.vtk> \t - Number of iterations for Laplace field, default 100"<<std::endl;

	exit(0);
}



int main(int argc, char **argv)
{
	std::cout<<"Version 1.0"<<std::endl;
	if (argc<4) usage(argv[0]);
	
	char* lvendo_filename;
	char* rvendo_filename;
	char* epi_filename;
        float voxel_lap = 0.2;
	int laplace_iter = 100;
        char* outimage1;
        char* outimage2;
        char* outimage3;
        
	
	//kernel for dilating interior boundary
	//for defining exterior boundary
	//new ext mask = dilated interior U old ext. mask
	int ext_dilate_kernelsize = 5; 
	
	
	
        std::cout<<"Parsing parameters"<<std::endl;
	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-lvendo")==0 )
		{
			lvendo_filename = argv[++c];
		}
		if( strcmp(argv[c],"-rvendo")==0 )
		{
			rvendo_filename = argv[++c];
		}
		else if ( strcmp(argv[c],"-epi")==0 )
		{
			epi_filename = argv[++c];
		}
		else if ( strcmp(argv[c],"-voxel_lap")==0 )
		{
			voxel_lap = atof(argv[++c]);
		}
		else if ( strcmp(argv[c],"-iter")==0 )
		{
			laplace_iter = atoi(argv[++c]);
		}
		else if ( strcmp(argv[c],"-o1")==0)
		{
			outimage1 = argv[++c];
		}
		else if ( strcmp(argv[c],"-o2")==0)
		{
			outimage2 = argv[++c];
		}
		else if ( strcmp(argv[c],"-o3")==0)
		{
			outimage3 = argv[++c];
		}
	}
	
	
        std::cout<<"Parameters parsed"<<std::endl;	

	
        std::cout<<"Checking for file existence"<<std::endl;	
	CommonTools::FileExists(lvendo_filename);
	CommonTools::FileExists(rvendo_filename);
	CommonTools::FileExists(epi_filename);

        std::cout<<"Loading "<<lvendo_filename<<std::endl;
	vtkSmartPointer<vtkPolyData> lvendo = vtkSmartPointer<vtkPolyData>::Take(
                CommonTools::LoadShapeFromFile(lvendo_filename) );
        std::cout<<"Loading "<<rvendo_filename<<std::endl;
	vtkSmartPointer<vtkPolyData> rvendo = vtkSmartPointer<vtkPolyData>::Take(
                CommonTools::LoadShapeFromFile(rvendo_filename) );
        std::cout<<"Loading "<<epi_filename<<std::endl;
	vtkSmartPointer<vtkPolyData> epi = vtkSmartPointer<vtkPolyData>::Take(
                CommonTools::LoadShapeFromFile(epi_filename) );


              		
	GenerateTemperatureField(epi, lvendo, rvendo, voxel_lap, laplace_iter, outimage1, outimage2, outimage3);


	return 0;
}



void GenerateTemperatureField( vtkPolyData* epi, vtkPolyData* lvendo, 
        vtkPolyData* rvendo, float VoxelSize/*=0.5*/, 
        int laplace_iterations/*=100*/,
        char* outfile_lv, char *outfile_rv, char* outfile_epi )
{
	vtkSmartPointer<vtkImageData> mask_lvendo = vtkSmartPointer<vtkImageData>::New();
	GenerateImageMask(mask_lvendo,0,-1, lvendo, epi, VoxelSize);

       	vtkSmartPointer<vtkImageData> mask_rvendo = vtkSmartPointer<vtkImageData>::New();
	GenerateImageMask(mask_rvendo,0,-1, rvendo, epi, VoxelSize);

	vtkSmartPointer<vtkImageData> mask_epi = vtkSmartPointer<vtkImageData>::New();
	GenerateImageMask(mask_epi,0,-1,epi,epi, VoxelSize);

#ifdef SAVE_DEBUG_INFO
	CommonTools::SaveImage(mask_lvendo,"mask_lvendo.vtk");
	CommonTools::SaveImage(mask_rvendo,"mask_rvendo.vtk");
	CommonTools::SaveImage(mask_epi,"mask_epi.vtk");
#endif
        
	//erode mask
        std::cout<<"Eroding endos"<<std::endl;
	vtkSmartPointer<vtkImageContinuousErode3D> erode = vtkSmartPointer<vtkImageContinuousErode3D>::New();
	erode->SetInputData (mask_lvendo);
	erode->SetKernelSize(5,5,5);
	erode->Update();
        mask_lvendo->DeepCopy(erode->GetOutput());

	erode->SetInputData (mask_rvendo);
	erode->Update();
        mask_rvendo->DeepCopy(erode->GetOutput());

	vtkSmartPointer<vtkImageData> field;


        //***********************************************************************
        //*
        //*
        //*
        //*
        //*
	std::cout<<"Computing Laplace field: LV endo vs EPI"<<std::endl;
	//combine masks
	vtkSmartPointer<vtkImageData> combinedMask = vtkSmartPointer<vtkImageData>::New();
	combinedMask->DeepCopy(mask_epi);

	char* mask_lvendo_data = static_cast<char*>(mask_lvendo->GetScalarPointer());
	char* mask_epi_data = static_cast<char*>(mask_epi->GetScalarPointer());
	char* combinedmask_data = static_cast<char*>(combinedMask->GetScalarPointer());
	for(int i=0; i<combinedMask->GetNumberOfPoints(); i++)
		combinedmask_data[i] = -( (mask_epi_data[i]+2)*mask_lvendo_data[i] );

	for(int i=0; i<combinedMask->GetNumberOfPoints(); i++)
		if ( combinedmask_data[i] == 0 ) 
                    combinedmask_data[i] = 1;
                else if( combinedmask_data[i] == 1 )
                    combinedmask_data[i] = 0;
                    
        
        //CommonTools::SaveImage(combinedMask,"combined_mask.vtk");

        field = vtkSmartPointer<vtkImageData>::New();
        CommonTools::laplace3D_voxelsize(combinedMask, field,laplace_iterations);

	CommonTools::SaveImage(field, outfile_lv);


        //***********************************************************************
        //*
        //*
        //*
        //*
        //*

	std::cout<<"Computing Laplace field: RV endo vs EPI"<<std::endl;
	//combine masks
	combinedMask->DeepCopy(mask_epi);

	mask_epi_data = static_cast<char*>(mask_epi->GetScalarPointer());
	char* mask_rvendo_data = static_cast<char*>(mask_rvendo->GetScalarPointer());
	combinedmask_data = static_cast<char*>(combinedMask->GetScalarPointer());
	for(int i=0; i<combinedMask->GetNumberOfPoints(); i++)
		combinedmask_data[i] = -( (mask_epi_data[i]+2)*mask_rvendo_data[i] );

	for(int i=0; i<combinedMask->GetNumberOfPoints(); i++)
		if ( combinedmask_data[i] == 0 ) 
                    combinedmask_data[i] = 1;
                else if( combinedmask_data[i] == 1 )
                    combinedmask_data[i] = 0;
                    

        //CommonTools::SaveImage(combinedMask,"combined_mask_correct.vtk");
        field = vtkSmartPointer<vtkImageData>::New();
        CommonTools::laplace3D_voxelsize(combinedMask, field,laplace_iterations);

	CommonTools::SaveImage(field, outfile_rv);

	
        
        
        //***********************************************************************
        //*
        //*
        //*
        //*
        //*

	std::cout<<"Computing Laplace field: RV+LV endo vs EPI"<<std::endl;
	//combine masks
	combinedMask->DeepCopy(mask_epi);

       	mask_lvendo_data = static_cast<char*>(mask_lvendo->GetScalarPointer());
	mask_rvendo_data = static_cast<char*>(mask_rvendo->GetScalarPointer());
	mask_epi_data = static_cast<char*>(mask_epi->GetScalarPointer());
	combinedmask_data = static_cast<char*>(combinedMask->GetScalarPointer());
	for(int i=0; i<combinedMask->GetNumberOfPoints(); i++)
            combinedmask_data[i] = -( (mask_epi_data[i]+2)*(mask_rvendo_data[i]+mask_lvendo_data[i]+1));


	for(int i=0; i<combinedMask->GetNumberOfPoints(); i++)
		if ( combinedmask_data[i] == 0 ) 
                    combinedmask_data[i] = 1;
                else if( combinedmask_data[i] == 1 )
                    combinedmask_data[i] = 0;
                    
        //CommonTools::SaveImage(combinedMask,"combined_mask.vtk");
        field = vtkSmartPointer<vtkImageData>::New();
        CommonTools::laplace3D_voxelsize(combinedMask, field,laplace_iterations);

	CommonTools::SaveImage(field, outfile_epi);
}






void GenerateImageMask( vtkImageData* res_image, float fg, float bg, vtkPolyData* endo, vtkPolyData* epi, float VoxelSize )
{
	//close the LV endo
	std::cout<<"Filling the holes"<<std::endl;

	vtkSmartPointer<vtkPolyData> endo_closed = vtkSmartPointer<vtkPolyData>::Take( CommonTools::CloseSurface(endo) );


	const int padding = 10; //space to add to shape extremes

	double origin[3];
	double spacing[3];
	int wextent[6];


	vtkSmartPointer<vtkImageData> inimage = vtkSmartPointer<vtkImageData>::New();

	std::cout<<"Generating an image"<<std::endl;

	spacing[0]=VoxelSize;
	spacing[1]=VoxelSize;
	spacing[2]=VoxelSize;

	double bounds[6];
	epi->GetBounds(bounds);
	double t[3];
	t[0] = -bounds[0]+padding*spacing[0];
	t[1] = -bounds[2]+padding*spacing[1];
	t[2] = -bounds[4]+padding*spacing[2];


	origin[0]=-t[0];
	origin[1]=-t[1];
	origin[2]=-t[2];
	int dims[3];
	epi->GetBounds(bounds);
	dims[0]=floor((bounds[1]-bounds[0])/spacing[0])+padding*2;
	dims[1]=floor((bounds[3]-bounds[2])/spacing[1])+padding*2;
	dims[2]=floor((bounds[5]-bounds[4])/spacing[2])+padding*2;

	inimage->SetOrigin(origin);
	inimage->SetSpacing(spacing);
	inimage->SetDimensions(dims);
//	inimage->SetScalarTypeToChar();
//	inimage->Update();
	inimage->AllocateScalars(VTK_CHAR,1);
	inimage->GetExtent(wextent);

	char* voxels = static_cast<char*>(inimage->GetScalarPointer());
	for(int i = 0; i<dims[0]*dims[1]*dims[2]; i++)
	{
		voxels[i] = fg;
	}

	//generate mask
	vtkSmartPointer<vtkPolyDataToImageStencil> filt = 
		vtkSmartPointer<vtkPolyDataToImageStencil>::New();


	filt->SetOutputOrigin( origin );
	filt->SetOutputSpacing( spacing );
	filt->SetOutputWholeExtent( wextent );
	filt->SetInputData ( endo_closed );
	filt->Update();

	vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
	stencil->SetStencilData(filt->GetOutput());
	stencil->SetBackgroundValue(bg);
	stencil->SetInputData (inimage);
	stencil->SetOutput(res_image);
	stencil->Update();
}

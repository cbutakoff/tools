/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Mesh segmentation using distance to skeleton along the gradiwent of the solution of the laplacian equation
*/
#include "CommonTools.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkType.h>
#include <vtkImageGradient.h>
#include <vtkAssignAttribute.h>
#include <vtkStreamTracer.h>
#include <vtkImageContinuousDilate3D.h>


//#include <itkImage.h>
//#include <itkImageFileReader.h>
//#include <itkImageFileWriter.h>


//custom made
//#include <itkImageToVTKImageFilter.h>
//#include <itkBinaryThinningImageFilter3D.h>


#include <iostream>




//------------------------------------------------------------------
int main( int argc, char *argv[] )
{	
	std::cout<<"MeshSegmentationLaplace v1.999"<<std::endl;
	std::cout<<"MeshSegmentationLaplace <image_mask.vtk> <mesh.vtk> <out_shape.vtk>"<<std::endl;
	std::cout<<"Image must have isotropic small voxels and is of type char"<<std::endl;

	if( argc<4 )
	{
		return EXIT_SUCCESS;
	}

	const char* filename_skeleton = "34875974_skeleton.vtk";
	const char* filename_skeleton_labeled = "34875974_skeleton_labeled.vtk";
	const char* filename_mask_pgm = "34875974_mask.pgm";
	const char* filename_mask_vtk = "34875974_mask.vtk";
	
	char cmd_label_skeleton[100];
	sprintf(cmd_label_skeleton, "LabelBranches3D %s %s", filename_skeleton, filename_skeleton_labeled);
	
	int c=1;
	const char *mask_filename = argv[c++];
	const char *mesh_filename = argv[c++];
	const char *out_filename = argv[c++];
	
	//-----------------------------------------------
	//
	//    reading files
	//
	//------------------------------------------------
	vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(mesh_filename) );
	
	//typedef itk::Image< char,  3 >  ImageTypeChar;
	//typedef itk::ImageFileReader< ImageTypeChar >  ImageReaderTypeChar;
	//typedef itk::ImageFileWriter< ImageTypeChar >  ImageWriterTypeChar;
	//typedef itk::ImageToVTKImageFilter<ImageTypeChar> ITKImageToVTKImageFilterTypeChar;
	//typedef itk::BinaryThinningImageFilter3D< ImageTypeChar, ImageTypeChar > ThinningFilterTypeChar;

	
	//-----------------------------------------------
	
	//-----------------------------------------------
	//
	//    Skeletonization
	//
	//------------------------------------------------	
	// 
	//ImageTypeChar::Pointer skeleton = ImageTypeChar::New();
	
	//skeletonize here
	// Define the thinning filter
	std::cout <<  "Thinning... " <<std::endl; 
//	ThinningFilterTypeChar::Pointer thinningFilter = ThinningFilterTypeChar::New();
//	thinningFilter->SetInputData( mask );
//	thinningFilter->Update();	
//	skeleton = thinningFilter->GetOutput();

	char command_vtk2pgm[100]; 
	char command_pgm2vtk[100]; 
	char command_thinning[100]; 
	sprintf( command_vtk2pgm, "itkvol2pgm %s %s 1", mask_filename, filename_mask_pgm );
	sprintf( command_pgm2vtk, "pgm2itkvol %s %s", filename_mask_pgm, filename_skeleton );
	sprintf( command_thinning, "pink.skelpar3d %s 9 50 %s", filename_mask_pgm, filename_mask_pgm );
	
	//convert mask to pgm
	std::cout<<"Running VTK2PGM..."<<std::endl;
	system(command_vtk2pgm);
	//run thinning
	std::cout<<"Running thinning..."<<std::endl;
	system(command_thinning);
	//convert thinned mask to vtk
	std::cout<<"Runnint pgm2vtk..."<<std::endl;
	system(command_pgm2vtk);
	//load thinned mask from vtk
	

	//load mask
	vtkSmartPointer<vtkImageData> mask = vtkSmartPointer<vtkImageData>::Take( 
		CommonTools::LoadImage(mask_filename));
	


	//-----------------------------------------------
	//
	//    Label the skeleton
	//
	//------------------------------------------------	
	// run skeletonization
	std::cout <<  "Label skeleton... " <<std::endl; 
	int res = system( cmd_label_skeleton );
	if( res == EXIT_FAILURE )
	{
		std::cout <<  "Error labeling skeleton. Exiting... " <<std::endl; 
		return EXIT_FAILURE;
	}
	

	
	// load the result from filename_skeleton_labeled
	std::cout <<  "Load Skeleton... " <<std::endl; 
	vtkSmartPointer<vtkImageData> skeleton = vtkSmartPointer<vtkImageData>::Take( 
		CommonTools::LoadImage(filename_skeleton_labeled));	
	skeleton->SetOrigin( mask->GetOrigin() );
	skeleton->SetSpacing( mask->GetSpacing() );
	//skeleton->Update();
	CommonTools::SaveImage( skeleton, filename_skeleton_labeled );
	//-----------------------------------------------

	//-----------------------------------------------
	//
	//    Solve laplace equation on the mask
	//
	//------------------------------------------------	
	//Solve Laplace equation
	vtkSmartPointer<vtkImageData> field = vtkSmartPointer<vtkImageData>::New();
	
	std::cout <<  "Mask itk->vtk... " <<std::endl; 
	
	
	//dilate the mask to make sure the mesh is inside the generated vector field
	std::cout<<"Dilating the mask..."<<std::endl;
	vtkSmartPointer<vtkImageContinuousDilate3D> dilate = vtkSmartPointer<vtkImageContinuousDilate3D>::New();
	dilate->SetInputData( mask );
	dilate->SetKernelSize(5,5,5);
	dilate->Update();

	//merge mask and the skeleton, relabel stuff
	/* Input values (c) Ruben Cardenes + Constantine Butakoff
	3 Outside domain
	1 Exterior boundary 
	0 Interior boundary
	2 Inside domain 
	*/ 


	int dims_mask[3];
	int dims_skeleton[3];
	dilate->GetOutput()->GetDimensions(dims_mask);
	skeleton->GetDimensions(dims_skeleton);
	if( dims_mask[0]!=dims_skeleton[0] || dims_mask[1]!=dims_skeleton[1] || dims_mask[2]!=dims_skeleton[2] )
	{
		std::cout<<"Dimensions of the skeleton and mask do not match. Aborting."<<std::endl;
		return EXIT_FAILURE;
	}

	
	char* voxels_mask = static_cast<char*>(dilate->GetOutput()->GetScalarPointer());
	char* voxels_skeleton = static_cast<char*>(skeleton->GetScalarPointer());
	for(int i = 0; i<dims_mask[0]*dims_mask[1]*dims_mask[2]; i++)
	{
		if( voxels_mask[i] == 0 ) //outside
			voxels_mask[i] = 3;
		else
		{
			if( voxels_skeleton[i] != 0 )
				voxels_mask[i] = 1;
			else
				voxels_mask[i] = 2;
		}
	}
	
	//CommonTools::SaveImage( skeleton, "skeleton.vtk"); 
	CommonTools::SaveImage( dilate->GetOutput(), "merged_mask.vtk"); 
	
	std::cout <<  "Solving Laplace eq.... " <<std::endl; 
	CommonTools::laplace3D_voxelsize(dilate->GetOutput(), field, 1000);

	CommonTools::SaveImage( field, "field.vtk"); 

	
	//compute gradient of DT
	std::cout<<"Computing gradient of the Laplacian solution"<<std::endl;
	vtkSmartPointer<vtkImageGradient> grad = vtkSmartPointer<vtkImageGradient>::New();
	grad->SetInputData(field);
	grad->SetDimensionality(3);
	grad->Update();
		
	//don't remember what the heck it is, but it's important
	vtkSmartPointer<vtkAssignAttribute> aa = vtkSmartPointer<vtkAssignAttribute>::New();
	aa->SetInputData( grad->GetOutput() );
	aa->Assign( vtkDataSetAttributes::SCALARS,
			vtkDataSetAttributes::VECTORS,
			vtkAssignAttribute::POINT_DATA );
	aa->Update();

	//generate streamlines for mesh vertices
	std::cout<<"Generate streamlines."<<std::endl;
	vtkSmartPointer<vtkStreamTracer> tracer = vtkSmartPointer<vtkStreamTracer>::New();
	tracer->SetInputData( aa->GetOutput() );
	tracer->SetSourceData( mesh );
	tracer->SetIntegrationDirectionToForward();
	tracer->SetIntegratorTypeToRungeKutta45();
	tracer->SetMaximumPropagation(20);
	tracer->Update();

	CommonTools::SavePolydata(tracer->GetOutput(),"tracer.vtk", true);

	//For every mesh point find the ending voxel in the labeled 
	//skeleton, get the label and assign to the vertex
	
	//------------------------------------------------	
	
	return EXIT_SUCCESS;
}



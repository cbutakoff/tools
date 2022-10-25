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


//#define USE_VMTK 

#include <vtkCleanPolyData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkIVWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
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


#include "VTKCommonTools.h"
#include "CommonTools.h"


#ifdef USE_VMTK
#include <vtkvmtkPolyDataSurfaceRemeshing.h>
#endif

#define FIELD_DT 0
#define FIELD_LAPLACE 1
//#define SAVE_DISTANCE_TFORM
//#define SAVE_DEBUG_DATA

//! generate points between epi and endo along vector field
void GenerateLayersAlongField(vtkPoints* layers, vtkPolyData* epi, vtkPolyData* endo, int nLayers, 
			int field_type=0, float VoxelSize=0.5, int laplace_iterations=1000);


//! generate mask given epi and endo
void GenerateImageMask(vtkImageData* res_image, float fg, float bg, vtkPolyData* endo, 
			vtkPolyData* epi, float VoxelSize=0.5);

//! Uniform remeshing. Requires admesh in path for mesh fixing
void UniformRemesh(vtkPolyData* mesh, float targetarea=1.0);

 
void usage(char *exe)
{
	std::cout<< "Generate layers of the atrium. " << std::endl;
	std::cout<<std::endl;
	std::cout<< "Options:" << std::endl;
	std::cout<<"-endo <mask.vtk> \t - endocardium mask (0-bg)"<<std::endl;
	std::cout<<"-epi <mask.vtk> \t - epicardium mask (0-bg)"<<std::endl;
	std::cout<<"-layers <int> \t - Number of layers of points including epi and endo. Min:2."<<std::endl;
	std::cout<<"-o <prefix> \t - prefix for the output filename with layers(files will be prefix00000.vtk, prefix00001.vtk,...)"<<std::endl;
	std::cout<<"-voxel <float> \t - Voxel size (isotropic), the smaller-the better (e.g. 0.2)"<<std::endl;
	std::cout<<"-iter <int> \t - Number of iterations for Laplace field, default 100"<<std::endl;
	std::cout<<"-nfaces <int> \t - Number of faces to retain during decimation (default 4000, negative - disables decimation)"<<std::endl;
	std::cout<<"-smooth_rad <float> \t - Smoothing std for mesh reconstruction (default 5)"<<std::endl;
	std::cout<<"-volmesh <mesh.vtk> \t - filename to store the volumetric mesh (optional)"<<std::endl;
	std::cout<<"-ta <float> \t - target area for vmtk remesh (only for layers, default 1, <0 - disable, requires admesh)"<<std::endl;

	exit(0);
}



int main(int argc, char **argv)
{
	std::cout<<"Version 2.0"<<std::endl;
	if (argc<4) usage(argv[0]);
	
	char* endo_filename;
	char* epi_filename;
	int nLayers=0;
	char* outshape_prefix;
	float VoxelSize = 0.5;
	int laplace_iter = 100;
	int nfaces = 4000;
	float smooth_radius = 5.0;
	char* volmesh_filename = NULL;
	float targetarea = 1.0;
        
	int hf_kernelsize = 3; //hole filling kernel
	
	//kernel for dilating interior boundary
	//for defining exterior boundary
	//new ext mask = dilated interior U old ext. mask
	int ext_dilate_kernelsize = 5; 
	
	
	int points_per_level = -1;
	
        std::cout<<"Parsing parameters"<<std::endl;
	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-endo")==0 )
		{
			endo_filename = argv[++c];
		}
		else if ( strcmp(argv[c],"-epi")==0 )
		{
			epi_filename = argv[++c];
		}
		else if ( strcmp(argv[c],"-layers")==0 )
		{
			nLayers = atoi(argv[++c]);
		}
		else if ( strcmp(argv[c],"-voxel")==0 )
		{
			VoxelSize = atof(argv[++c]);
		}
		else if ( strcmp(argv[c],"-iter")==0 )
		{
			laplace_iter = atoi(argv[++c]);
		}
		else if ( strcmp(argv[c],"-o")==0)
		{
			outshape_prefix = argv[++c];
		}
		else if ( strcmp(argv[c],"-nfaces")==0)
		{
			nfaces = atoi(argv[++c]);
		}
		else if ( strcmp(argv[c],"-smooth_rad")==0)
		{
			smooth_radius = atof(argv[++c]);
		}
		else if ( strcmp(argv[c],"-volmesh")==0)
		{
			volmesh_filename = argv[++c];
		}
		else if ( strcmp(argv[c],"-ta")==0)
		{
			 targetarea = atof(argv[++c]);
		}
	}
	
	
        std::cout<<"Parameters parsed"<<std::endl;	
	if(nLayers<2)
	{
		std::cout<<"Number of layers must be >=2"<<std::endl;
	}
	
        std::cout<<"Checking for file existence"<<std::endl;	
	CommonTools::FileExists(endo_filename);
	CommonTools::FileExists(epi_filename);

        std::cout<<"Loading "<<endo_filename<<std::endl;
	vtkSmartPointer<vtkImageData> mask_endo = vtkSmartPointer<vtkImageData>::Take(	CommonTools::LoadImage(endo_filename) );
        std::cout<<"Loading "<<epi_filename;
	vtkSmartPointer<vtkImageData> mask_epi = vtkSmartPointer<vtkImageData>::Take( CommonTools::LoadImage(epi_filename) );

	double vx, vy, vz;
	mask_endo->GetSpacing(vx, vy, vz);

	double ox, oy, oz;
	mask_endo->GetOrigin(ox, oy, oz);

	std::cout<<"Cropping the masks..."<<std::endl;
	{
		//run marching cubes to get epi mesh
		vtkSmartPointer<vtkImageMarchingCubes> mc = vtkSmartPointer<vtkImageMarchingCubes>::New();
		mc->SetInputData(mask_epi);
		mc->Update();
		//CommonTools::SaveShapeToFile(mc->GetOutput(),"epi_mc.vtk");
		
		int image_extent[6];
		mask_epi->GetExtent( image_extent );
		
		double 	bounds[6];
		mc->GetOutput()->GetBounds( bounds );
		bounds[0] = bounds[0]-5*vx;
		bounds[1] = bounds[1]+5*vx;
		bounds[2] = bounds[2]-5*vy;
		bounds[3] = bounds[3]+5*vy;
		bounds[4] = bounds[4]-5*vz;
		bounds[5] = bounds[5]+5*vz;
		

		int minX = vtkMath::Floor(fabs(bounds[0]-ox)/vx);
		int minY = vtkMath::Floor(fabs(bounds[2]-oy)/vy);
		int minZ = vtkMath::Floor(fabs(bounds[4]-oz)/vz);
		int maxX = vtkMath::Ceil(fabs(bounds[1]-ox)/vx);
		int maxY = vtkMath::Ceil(fabs(bounds[3]-oy)/vy);
		int maxZ = vtkMath::Ceil(fabs(bounds[5]-oz)/vz);

		if ( minX<image_extent[0] ) minX=image_extent[0];
		if ( minY<image_extent[2] ) minY=image_extent[2];
		if ( minZ<image_extent[4] ) minZ=image_extent[4];
		if ( maxX>image_extent[1] ) maxX=image_extent[1];
		if ( maxY>image_extent[3] ) maxY=image_extent[3];
		if ( maxZ>image_extent[5] ) maxZ=image_extent[5];

		std::cout<<"Clipping: "<<minX<<" "<<maxX<<" "<<minY
			<<" "<<maxY<<" "<<minZ<<" "<<maxZ<<" "<<std::endl;

		vtkSmartPointer<vtkImageClip> clip = vtkSmartPointer<vtkImageClip>::New();
		clip->SetInputData( mask_epi );
		clip->SetOutputWholeExtent( minX, maxX, minY, maxY, minZ, maxZ );
		clip->ClipDataOn();
		clip->Update();
		
		vtkSmartPointer<vtkImageNormalize> mynorm = vtkSmartPointer<vtkImageNormalize>::New();
		mynorm->SetInputData(clip->GetOutput());
		mynorm->Update();

		mask_epi->DeepCopy( mynorm->GetOutput() );

		clip->SetInputData( mask_endo );
		clip->Update();
		mynorm->Update();

		mask_endo->DeepCopy( mynorm->GetOutput() );

#ifdef SAVE_DEBUG_DATA
		CommonTools::SaveImage( mask_epi, "epi_clipped.vtk");
		CommonTools::SaveImage( mask_endo, "endo_clipped.vtk");
#endif 
	}

	std::cout<<"Hole filling in the masks..."<<std::endl;
	{
		vtkSmartPointer<vtkImageContinuousDilate3D> dilate = vtkSmartPointer<vtkImageContinuousDilate3D>::New();
		dilate->SetInputData(mask_epi);
		dilate->SetKernelSize(hf_kernelsize,hf_kernelsize,hf_kernelsize);
		dilate->Update();
		
		vtkSmartPointer<vtkImageContinuousErode3D> erode = vtkSmartPointer<vtkImageContinuousErode3D>::New();
		erode->SetInputData(dilate->GetOutput());
		erode->SetKernelSize(hf_kernelsize,hf_kernelsize,hf_kernelsize);
		erode->Update();
		
		mask_epi->DeepCopy(erode->GetOutput());

		
		dilate->SetInputData(mask_endo);
		dilate->Update();
		erode->Update();
		mask_endo->DeepCopy(erode->GetOutput());
		
#ifdef SAVE_DEBUG_DATA
		CommonTools::SaveImage( mask_epi, "epi_hf.vtk");
		CommonTools::SaveImage( mask_endo, "endo_hf.vtk");
#endif
	}

	std::cout<<"Resampling the masks: "<<vx/VoxelSize<<" "<<vy/VoxelSize
	<<" "<<vz/VoxelSize<<std::endl;
	{
		vtkSmartPointer<vtkImageResample> resamp = vtkSmartPointer<vtkImageResample>::New();
		resamp->SetInputData( mask_endo );
		resamp->SetAxisMagnificationFactor(0,vx/VoxelSize);
		resamp->SetAxisMagnificationFactor(1,vy/VoxelSize);
		resamp->SetAxisMagnificationFactor(2,vz/VoxelSize);
		resamp->SetAxisOutputSpacing(0, VoxelSize);
		resamp->SetAxisOutputSpacing(1, VoxelSize);
		resamp->SetAxisOutputSpacing(2, VoxelSize);
		resamp->SetInterpolationModeToCubic();
	
		resamp->Update();
		
		//remove negative pixels
		vtkSmartPointer<vtkImageThreshold> thold = vtkSmartPointer<vtkImageThreshold>::New();
		thold->SetInputData(resamp->GetOutput());
		thold->ThresholdByLower(0);
		thold->ReplaceInOn();
		thold->ReplaceOutOff();
		thold->SetInValue(0);
		thold->Update();

		vtkSmartPointer<vtkImageThreshold> thold1 = vtkSmartPointer<vtkImageThreshold>::New();
		thold1->SetInputData(thold->GetOutput());
		thold1->ThresholdByUpper(1);
		thold1->ReplaceInOn();
		thold1->ReplaceOutOff();
		thold1->SetInValue(1);
		thold1->Update();
		
		mask_endo->DeepCopy(thold1->GetOutput());
		
		resamp->SetInputData( mask_epi );
		resamp->Update();
		thold->Update();
		thold1->Update();
		mask_epi->DeepCopy(thold1->GetOutput());

#ifdef SAVE_DEBUG_DATA
		CommonTools::SaveImage( mask_epi, "epi_res.vtk");
		CommonTools::SaveImage( mask_endo, "endo_res.vtk");
#endif
	}

	{
		std::cout<<"Making epi bigger than endo..."<<std::endl;
		vtkSmartPointer<vtkImageContinuousDilate3D> dilate = vtkSmartPointer<vtkImageContinuousDilate3D>::New();
		dilate->SetInputData(mask_endo);
		dilate->SetKernelSize(ext_dilate_kernelsize,ext_dilate_kernelsize,ext_dilate_kernelsize);
		dilate->Update();
		
		vtkSmartPointer<vtkImageMathematics> math = vtkSmartPointer<vtkImageMathematics>::New();
		math->SetOperationToMax();
		math->SetInput1Data( mask_epi );
		math->SetInput2Data( dilate->GetOutput() );
		math->Update();
		mask_epi->DeepCopy(math->GetOutput());
	}
	
#ifdef SAVE_DEBUG_DATA
	CommonTools::SaveImage(mask_epi,"epi_big.vtk");
#endif

	vtkSmartPointer<vtkPolyData> endo_mesh = vtkSmartPointer<vtkPolyData>::New();
	{
		std::cout<<"Smooth endo mask and get endo surface..."<<std::endl;
		
		std::cout<<"Smoothing..."<<std::endl;
		vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
		smooth->SetInputData(mask_endo);
		smooth->SetStandardDeviation( smooth_radius );
		smooth->Update();
		
		std::cout<<"Running marching cubes..."<<std::endl;
		vtkSmartPointer<vtkImageMarchingCubes> marchingcubes = vtkSmartPointer<vtkImageMarchingCubes>::New();
		marchingcubes->SetInputData(smooth->GetOutput());
		marchingcubes->ComputeNormalsOn();
		marchingcubes->SetValue(0,0.5);
		marchingcubes->Update();
		
		vtkSmartPointer<vtkPolyData> tempPD;
		if( nfaces>0 )
		{
			std::cout<<"Decimating... (Meshlab)"<<std::endl;
			char temp_filename[100];
			strcpy( temp_filename, "costa_nhfh38hf87.ply" );
		
			CommonTools::SaveShapeToFile(marchingcubes->GetOutput(),temp_filename);
			CommonTools::GenerateDecimationScript("decimation_script_196239462.mlx", nfaces);
			char cmdline[1000];
			sprintf(cmdline,"meshlabserver -i %s -o %s -s decimation_script_196239462.mlx", temp_filename, temp_filename);
			int rubbish = system(cmdline);
			tempPD = vtkSmartPointer< vtkPolyData > ::Take( CommonTools::LoadShapeFromFile( temp_filename) );
			remove("decimation_script_196239462.mlx");
			remove(temp_filename);		
		}
		else
			tempPD = marchingcubes->GetOutput();
						
		std::cout<<"Cleaning the mesh..."<<std::endl;

		vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
		clean->SetInputData(tempPD);
		clean->Update();
		endo_mesh->SetPoints(clean->GetOutput()->GetPoints());
		endo_mesh->SetPolys(clean->GetOutput()->GetPolys());
		//endo_mesh->Update();

		vtkSmartPointer<vtkPolyDataNormals> normalgen = vtkSmartPointer<vtkPolyDataNormals>::New();
		normalgen->SetInputData(endo_mesh);
		normalgen->SplittingOff();
		normalgen->FlipNormalsOn();
		normalgen->Update();

		endo_mesh->DeepCopy(normalgen->GetOutput());

		//generate normals
        if(targetarea>0)
        {
            UniformRemesh(endo_mesh, targetarea);
        }
                
                
	}
	
	
	vtkSmartPointer<vtkPolyData> epi_mesh = vtkSmartPointer<vtkPolyData>::New();
	{
		std::cout<<"Smooth epi mask and get epi surface..."<<std::endl;
		
		std::cout<<"Smoothing..."<<std::endl;
		vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
		smooth->SetInputData(mask_epi);
		smooth->SetStandardDeviation(smooth_radius);
		smooth->Update();
		
		std::cout<<"Running marching cubes..."<<std::endl;
		vtkSmartPointer<vtkImageMarchingCubes> marchingcubes = vtkSmartPointer<vtkImageMarchingCubes>::New();
		marchingcubes->SetInputData(smooth->GetOutput());
		marchingcubes->ComputeNormalsOn();
		marchingcubes->SetValue(0,0.5);
		marchingcubes->Update();

		vtkSmartPointer<vtkPolyData> tempPD;

		if(nfaces>0)
		{
			std::cout<<"Decimating... (Meshlab)"<<std::endl;
			char temp_filename[100];
			strcpy( temp_filename, "costa_nhfh38hf87.ply" );
		
			CommonTools::SaveShapeToFile(marchingcubes->GetOutput(),temp_filename);
			CommonTools::GenerateDecimationScript("decimation_script_196239462.mlx",1000);
			char cmdline[1000];
			sprintf(cmdline,"meshlabserver -i %s -o %s -s decimation_script_196239462.mlx", temp_filename, temp_filename);
			int rubbish = system(cmdline);
			tempPD = vtkSmartPointer< vtkPolyData > ::Take( CommonTools::LoadShapeFromFile( temp_filename) );
			remove("decimation_script_196239462.mlx");
			remove(temp_filename);		
		}
		else
			tempPD = marchingcubes->GetOutput();

		std::cout<<"Cleaning the mesh..."<<std::endl;

		vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
		clean->SetInputData(tempPD);
		clean->Update();
		epi_mesh->SetPoints(clean->GetOutput()->GetPoints());
		epi_mesh->SetPolys(clean->GetOutput()->GetPolys());
		//epi_mesh->Update();

		//generate normals
		vtkSmartPointer<vtkPolyDataNormals> normalgen = vtkSmartPointer<vtkPolyDataNormals>::New();
		normalgen->SetInputData(epi_mesh);
		normalgen->SplittingOff();
		normalgen->Update();

		epi_mesh->DeepCopy(normalgen->GetOutput());
	}
		

#ifdef SAVE_DEBUG_DATA
	CommonTools::SaveShapeToFile(endo_mesh,"endo_mesh.vtk");	
	CommonTools::SaveShapeToFile(epi_mesh,"epi_mesh.vtk");	
#endif
		
		
	/////////////////////////////////////////////////////////
	//
	//  At this point we have epi_mask bigger than endo_mask
	//  both are float with pixels in [0,1]
	//	endo_mesh - is the decimated smoothed mesh of the endocardium
	//
	//////////////////////////////////////////////////////////////////
		
	//free images
	mask_endo = vtkSmartPointer<vtkImageData>::New();
	mask_epi = vtkSmartPointer<vtkImageData>::New();

	const int nPointsPerLayer = endo_mesh->GetNumberOfPoints();
	vtkSmartPointer<vtkPoints> layers =  vtkSmartPointer<vtkPoints>::New();
	layers->SetNumberOfPoints(nLayers*nPointsPerLayer);
	
	GenerateLayersAlongField(layers, epi_mesh, endo_mesh, nLayers, FIELD_LAPLACE, VoxelSize, laplace_iter);




	std::cout<<"Generate layers"<<std::endl;
	for(int i=0; i<nLayers; i++)
	{
		char filename[100];
		
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		pts->SetNumberOfPoints(nPointsPerLayer);

		for(int j=0; j<nPointsPerLayer; j++)
		{
			pts->SetPoint(j, layers->GetPoint(i*nPointsPerLayer+j));
		}
		
		vtkSmartPointer<vtkPolyData> layer = vtkSmartPointer<vtkPolyData>::New();
		layer->SetPoints(pts);
		layer->SetPolys(endo_mesh->GetPolys());
		//layer->Update();

		

        if(targetarea>0)
        {
            UniformRemesh(layer, targetarea);
        }
        
		vtkSmartPointer<vtkPolyDataNormals> normalgen = vtkSmartPointer<vtkPolyDataNormals>::New();
		normalgen->SetInputData(layer);
		normalgen->SplittingOff();
		normalgen->Update();

		//layer->DeepCopy(normalgen->GetOutput());

        std::cout<<"Smoothing layer"<<std::endl;
        sprintf( filename, "%s%05d.vtk", outshape_prefix, i );
		CommonTools::SaveShapeToFile(normalgen->GetOutput(), filename);
                
	}


	if(volmesh_filename!=NULL)
	{
		//generate elements
		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkShortArray> celllayerid = vtkSmartPointer<vtkShortArray>::New();
		vtkSmartPointer<vtkShortArray> pointlayerid = vtkSmartPointer<vtkShortArray>::New();
		vtkSmartPointer<vtkShortArray> pointSubpartId = vtkSmartPointer<vtkShortArray>::New();
		celllayerid->SetName("CellLayerID");
		pointlayerid->SetName("PointLayerID");
		pointSubpartId->SetName("subpartID");
		
		std::cout<<"Generating wedge elements"<<std::endl;
		for(int i=0; i<endo_mesh->GetNumberOfCells(); i++ )
		{
			vtkCell *endocell = endo_mesh->GetCell(i);
			if (endocell->GetNumberOfPoints()!=3)
			{
				std::cout<<"Nontriangular faces in endo mesh. Do something!"<<std::endl;
				exit(0);
			}

			vtkIdList* ptIds = endocell->GetPointIds();
			

			for (int layer=0; layer<nLayers-1; layer++)
			{
				vtkSmartPointer<vtkWedge> wedge = vtkSmartPointer<vtkWedge>::New();
				for(int j=0; j<3; j++)
				{
					wedge->GetPointIds()->SetId(j, layer*nPointsPerLayer+ptIds->GetId(j));
					wedge->GetPointIds()->SetId(j+3, (layer+1)*nPointsPerLayer+ptIds->GetId(j));
				}
				cells->InsertNextCell(wedge);
				celllayerid->InsertNextValue(layer);
			}

		}
		
		endo_mesh->SetPolys(cells);
		

		std::cout<<"Saving"<<std::endl;
		vtkSmartPointer<vtkUnstructuredGrid> output = vtkSmartPointer<vtkUnstructuredGrid>::New();
	//	output->Allocate(cells->GetNumberOfCells());
		output->SetPoints(layers);

		//generate pointlayerid
		pointlayerid->SetNumberOfValues(output->GetNumberOfPoints());
		for(int i=0; i<nPointsPerLayer; i++)
		{
			for (int j=0; j<nLayers; j++)
			{
				pointlayerid->SetValue(j*nPointsPerLayer+i,j);
			}
		}


		//generate subpartID
		pointSubpartId->SetNumberOfValues(pointlayerid->GetNumberOfTuples());
		for(int i=0;i<pointlayerid->GetNumberOfTuples(); i++)
		{
			vtkIdType id = pointlayerid->GetValue(i);
			if (id==0)
			{
				pointSubpartId->SetValue(i,1);
			} 
			else if (id==nLayers-1)
			{
				pointSubpartId->SetValue(i,2);
			}
			else
			{
				pointSubpartId->SetValue(i,0);
			}
		}
			

		//save everything
		output->SetCells(VTK_WEDGE,cells);
		output->GetCellData()->AddArray(celllayerid);
		output->GetPointData()->AddArray(pointlayerid);
		output->GetPointData()->AddArray(pointSubpartId);
		
		
		vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		writer->SetFileName(volmesh_filename);
		writer->SetInputData(output);
		writer->SetFileTypeToBinary();
		writer->Update();
	}

	return 0;
}



void GenerateLayersAlongField( vtkPoints* layers, vtkPolyData* epi, vtkPolyData* endo, int nLayers, int field_type/*=0*/, float VoxelSize/*=0.5*/, int laplace_iterations/*=100*/ )
{
	const int nPointsPerLayer = endo->GetNumberOfPoints();

	vtkSmartPointer<vtkImageData> mask_endo = vtkSmartPointer<vtkImageData>::New();
	GenerateImageMask(mask_endo,0,-1,endo,epi, VoxelSize);

	//blVTKHelperTools::SaveImage(mask_endo,"maskendo.vtk");


	vtkSmartPointer<vtkImageData> mask_epi = vtkSmartPointer<vtkImageData>::New();
	GenerateImageMask(mask_epi,0,-1,epi,epi, VoxelSize);

	//erode mask
	vtkSmartPointer<vtkImageContinuousErode3D> erode = vtkSmartPointer<vtkImageContinuousErode3D>::New();
	erode->SetInputData(mask_endo);
	erode->SetKernelSize(5,5,5);
	erode->Update();



	//blVTKHelperTools::SaveImage(erode->GetOutput(),"stencil.vtk");

	vtkSmartPointer<vtkImageData> field;
	vtkSmartPointer<vtkImageEuclideanDistance> dist = vtkSmartPointer<vtkImageEuclideanDistance>::New();


	//std::cout<<"Dilating epi mask..."<<std::endl;
	//vtkSmartPointer<vtkImageContinuousDilate3D> dilate = vtkSmartPointer<vtkImageContinuousDilate3D>::New();
	//dilate->SetInputData(mask_epi);
	//dilate->SetKernelSize(5,5,5);
	//dilate->Update();

	std::cout<<"Computing Laplace field"<<std::endl;
	//combine masks
	vtkSmartPointer<vtkImageData> combinedMask = vtkSmartPointer<vtkImageData>::New();
	combinedMask->DeepCopy(mask_epi);

	char* mask_endo_data = static_cast<char*>(erode->GetOutput()->GetScalarPointer());
	char* mask_epi_data = static_cast<char*>(mask_epi->GetScalarPointer());
	char* combinedmask_data = static_cast<char*>(combinedMask->GetScalarPointer());
	for(int i=0; i<combinedMask->GetNumberOfPoints(); i++)
		combinedmask_data[i] = -( (mask_epi_data[i]+2)*mask_endo_data[i] );

	vtkSmartPointer<vtkImageCast> field_caster = vtkSmartPointer<vtkImageCast>::New();

	double* fielddata;

	switch(field_type)
	{
	case FIELD_DT:
		std::cout<<"Computing distance transform"<<std::endl;
		dist->SetInputData(erode->GetOutput());
		dist->InitializeOn();
		dist->ConsiderAnisotropyOff();
		dist->Update();

		fielddata = static_cast<double*>(dist->GetOutput()->GetScalarPointer());
		for(int i=0; i<dist->GetOutput()->GetNumberOfPoints(); i++)
			if( mask_epi_data[i]<0 )
				fielddata[i] = 1;
			else
				fielddata[i] = fielddata[i]/100000;
				
		
		field_caster->SetInputData( dist->GetOutput() );
		field_caster->SetOutputScalarTypeToFloat();
		field_caster->Update();

		dist = NULL;

		field = field_caster->GetOutput();
		
		
		
		break;
	case FIELD_LAPLACE:
		{
			//blVTKHelperTools::SaveImage(combinedMask,"combined_mask.vtk");
			field = vtkSmartPointer<vtkImageData>::New();
			CommonTools::laplace3D_voxelsize(combinedMask, field,laplace_iterations);
		}
		break;
	default:
		std::cout<<"Unknown field type"<<std::endl;
		return;
	}

	//CommonTools::SaveImage(field,"field.vtk");

	//field = vtkSmartPointer<vtkImageData>::Take( CommonTools::LoadImage("distance.vtk") );

	std::cout<<"Cleaning up the mask"<<std::endl;
	
	//transform the mask to have 1 inside and 0 outside

	vtkSmartPointer<vtkImageMathematics> epi_shifted = vtkSmartPointer<vtkImageMathematics>::New();
	epi_shifted->SetInput1Data( mask_epi );
	epi_shifted->SetOperationToAddConstant();
	epi_shifted->SetConstantC(1);
	epi_shifted->Update();


	vtkSmartPointer<vtkImageCast> caster = vtkSmartPointer<vtkImageCast>::New();
	caster->SetInputData( epi_shifted->GetOutput() );
	caster->SetOutputScalarTypeToFloat();
	caster->Update();

	std::cout<<"vtkImageMathematics - epi_inverted"<<std::endl;
	vtkSmartPointer<vtkImageMathematics> epi_inverted = vtkSmartPointer<vtkImageMathematics>::New();
	epi_inverted->SetInput1Data( caster->GetOutput() );
	epi_inverted->SetOperationToAddConstant();
	epi_inverted->SetConstantC(-1);
	epi_inverted->Update();



	std::cout<<"vtkImageMathematics - im_math_mask"<<std::endl;
	vtkSmartPointer<vtkImageMathematics> im_math_mask = vtkSmartPointer<vtkImageMathematics>::New();
	im_math_mask->SetInput1Data( field );
	im_math_mask->SetInput2Data( caster->GetOutput() );
	im_math_mask->SetOperationToMultiply();
	im_math_mask->Update();

	std::cout<<"vtkImageMathematics - im_math_setbg"<<std::endl;
	vtkSmartPointer<vtkImageMathematics> im_math_setbg = vtkSmartPointer<vtkImageMathematics>::New();
	im_math_setbg->SetInput1Data( im_math_mask->GetOutput() );
	im_math_setbg->SetInput2Data( epi_inverted->GetOutput() );
	im_math_setbg->SetOperationToSubtract();
	im_math_setbg->Update();

	//CommonTools::SaveImage(im_math_setbg->GetOutput(),"product.vtk");


	std::cout<<"Computing gradient of the distance transform"<<std::endl;
	vtkSmartPointer<vtkImageGradient> grad = vtkSmartPointer<vtkImageGradient>::New();
	grad->SetInputData(field);
	grad->SetDimensionality(3);
	grad->Update();

	//vtkImageData* grad_image = grad->GetOutput();
		

	vtkSmartPointer<vtkAssignAttribute> aa = vtkSmartPointer<vtkAssignAttribute>::New();
	aa->SetInputData( grad->GetOutput() );
	aa->Assign( vtkDataSetAttributes::SCALARS,
			vtkDataSetAttributes::VECTORS,
			vtkAssignAttribute::POINT_DATA );
	aa->Update();


	std::cout<<"Generate streamlines."<<std::endl;
	vtkSmartPointer<vtkStreamTracer> tracer = vtkSmartPointer<vtkStreamTracer>::New();
	tracer->SetInputData( aa->GetOutput() );
	tracer->SetSourceData( endo );
	tracer->SetIntegrationDirectionToForward();
	tracer->SetIntegratorTypeToRungeKutta45();
	tracer->SetMaximumPropagation(200);
	tracer->Update();


	//CommonTools::SavePolydata(tracer->GetOutput(),"tracer.vtk");

	//resample streamlines
	vtkSmartPointer<vtkSplineFilter> spline = vtkSmartPointer<vtkSplineFilter>::New();
	spline->SetInputData(tracer->GetOutput());
	spline->SetSubdivideToSpecified();
	spline->SetNumberOfSubdivisions(nLayers-1);
	spline->Update();

	//CommonTools::SavePolydata(spline->GetOutput(),"tracer_resampled.vtk");

	std::cout<<"Streamlines: "<<spline->GetOutput()->GetLines()->GetNumberOfCells()<<std::endl;


	//move every point from endo to epi along the gradient

	std::cout<<"Generating point layers"<<std::endl;

	for (int i=0; i<spline->GetOutput()->GetLines()->GetNumberOfCells(); i++)
	{		
		vtkCell* line = spline->GetOutput()->GetCell(i);		

		for(int j=0;j<nLayers; j++)
		{
			double pt[3];

			line->GetPoints()->GetPoint(j, pt); 

			layers->SetPoint(i+j*nPointsPerLayer,pt);
		}

	}
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
	inimage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
//	inimage->Update();
//	inimage->AllocateScalars();
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
	filt->SetInputData( endo_closed );
	filt->Update();

	vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
	stencil->SetStencilData(filt->GetOutput());
	stencil->SetBackgroundValue(bg);
	stencil->SetInputData(inimage);
	stencil->SetOutput(res_image);
	stencil->Update();
}


void UniformRemesh(vtkPolyData* mesh, float targetarea)
{
    char filename[100];
    //post process: fix and uniformly remesh
    sprintf( filename, "tempshape3274899816.stl" );
    CommonTools::SaveShapeToFile(mesh, filename);

    char fixcmd[500];
    std::cout<<"Running admesh.."<<std::endl;
    sprintf(fixcmd,"admesh -n -u -f -d %s -b%s -t0.1 -i50",filename,filename);
    int rubbish = system(fixcmd);

    vtkSmartPointer<vtkPolyData> temp_mesh = vtkSmartPointer<vtkPolyData>::Take( 
            CommonTools::LoadShapeFromFile(filename) );

    remove(filename);

    //remesh
#ifdef USE_VMTK
    vtkSmartPointer<vtkvmtkPolyDataSurfaceRemeshing> remesh = 
            vtkSmartPointer<vtkvmtkPolyDataSurfaceRemeshing>::New();
    remesh->SetInputData(temp_mesh);
    remesh->SetNumberOfIterations(10);
    remesh->SetElementSizeModeToTargetArea();
    remesh->SetTargetArea(targetarea);
    remesh->Update(); 
    
    mesh->DeepCopy(remesh->GetOutput());
#endif
}

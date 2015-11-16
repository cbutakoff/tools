/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Generate volumetric mesh of Left Ventricle. 
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


#define FIELD_DT 0
#define FIELD_LAPLACE 1
//#define SAVE_DISTANCE_TFORM

//! generate points between epi and endo along surface normals
void GenerateLayersAlongNormals(vtkPoints* layers, vtkPolyData* epi, vtkPolyData* endo, int nLayers);
//! generate points between epi and endo along vector field
void GenerateLayersAlongField(vtkPoints* layers, vtkPolyData* epi, vtkPolyData* endo, int nLayers, int field_type=0, float VoxelSize=0.5, int laplace_iterations=100);
//! generate mask given epi and endo
void GenerateImageMask(vtkImageData* res_image, float fg, float bg, vtkPolyData* endo, vtkPolyData* epi, float VoxelSize=0.5);

//! generate local coordinates 
void GenerateLocalCoordinateCircLongit(vtkUnstructuredGrid* volmesh, int nLayers, int nPointsPerLevel);
//! generate local coordinates 
void GenerateLocalCoordinateRadial(vtkUnstructuredGrid* volmesh, int nLayers);

 
void usage(char *exe)
{
	std::cout<< "The LV mesh must have \"subpartID\" scalar array of type short with 2 for epicardium. Everything else should be endocardium." << std::endl;
	std::cout<<std::endl;
	std::cout<< "Options:" << std::endl;
	std::cout<<"-endo <polydata.vtk> \t - Endocardium mesh"<<std::endl;
	std::cout<<"-endo <polydata.vtk> \t - Endocardium mesh"<<std::endl;
	std::cout<<"-lv <polydata.vtk> \t - Full left ventricular mesh with subpartID"<<std::endl;
	std::cout<<"-layers <int> \t - Number of layers of points including epi and endo. Min:2."<<std::endl;
	std::cout<<"-o <polydata.vtk> \t - Output filename"<<std::endl;
	std::cout<<"-method normals|distance|laplace \t - first is based on surface normals, the second-on distance image transform"<<std::endl;
	std::cout<<"-coordinates <points_per_level> \t - if specified, coordinate system for every point will be generated (normally points_per_level=48)"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"If method = distance|laplace specify:"<<std::endl;
	std::cout<<"-voxel <float> \t - Voxel size (isotropic), the smaller-the better"<<std::endl;
	std::cout<<"-iter <int> \t - Number of iterations for Laplace field, default 100"<<std::endl;
	std::cout<<"Saves distance transform to distance.vtk unless disabled at compile time"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"If method = normals make sure the endo normals point inwards"<<std::endl;

	exit(0);
}



int main(int argc, char **argv)
{
	std::cout<<"Version 1.4"<<std::endl;
	if (argc<4) usage(argv[0]);
	
	char* endo_filename;
	char* fulllv_filename;
	int nLayers=0;
	char* outshape_filename;
	float VoxelSize = 0.5;
	int laplace_iter = 100;
	
	int points_per_level = -1;
	
	typedef enum {methodNormals, methodDistance, methodLaplace} VolumetrizationMethod;

	VolumetrizationMethod method=methodDistance;

	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-endo")==0 )
		{
			endo_filename = argv[++c];
		}
		else if ( strcmp(argv[c],"-lv")==0 )
		{
			fulllv_filename = argv[++c];
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
			outshape_filename = argv[++c];
		}
		else if ( strcmp(argv[c],"-coordinates")==0)
		{
			points_per_level = atoi(argv[++c]);
		}
		else if ( strcmp(argv[c],"-method")==0)
		{
			const char *str_type = argv[++c];
			if( strcmp(str_type , "normals")==0 ) 
				method = methodNormals;
			else if( strcmp(str_type , "distance")==0 ) 
				method = methodDistance;
			else if( strcmp(str_type , "laplace")==0 ) 
				method = methodLaplace;
		}
	}
	
	if(nLayers<2)
	{
		std::cout<<"Number of layers must be >=2"<<std::endl;
	}
	
	CommonTools::FileExists(endo_filename);
	CommonTools::FileExists(fulllv_filename);

	vtkSmartPointer<vtkPolyData> endo = vtkSmartPointer<vtkPolyData>::Take(	CommonTools::LoadShapeFromFile(endo_filename) );
	vtkSmartPointer<vtkPolyData> fulllv = vtkSmartPointer<vtkPolyData>::Take( CommonTools::LoadShapeFromFile(fulllv_filename) );

	vtkSmartPointer<vtkPolyDataNormals> normalgen = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalgen->SetInputData(endo);
	normalgen->SplittingOff();
	normalgen->ComputeCellNormalsOff();
	normalgen->ComputePointNormalsOn();
	normalgen->Update();

		
	fulllv->GetPointData()->SetActiveScalars("subpartID");
	vtkSmartPointer<vtkPolyData> epi = vtkSmartPointer<vtkPolyData>::Take( CommonTools::GetShapeSubSurface(fulllv,2) );


	const int nPointPerLayer = endo->GetNumberOfPoints();
	vtkSmartPointer<vtkPoints> layers =  vtkSmartPointer<vtkPoints>::New();
	layers->SetNumberOfPoints(nLayers*nPointPerLayer);

	//generate point layers
	switch (method)
	{
	case methodNormals:
		GenerateLayersAlongNormals(layers, epi, normalgen->GetOutput(), nLayers);
		break;
	case methodDistance:
		GenerateLayersAlongField(layers, epi, endo, nLayers, FIELD_DT, VoxelSize, laplace_iter);
		break;
	case methodLaplace:
		GenerateLayersAlongField(layers, epi, endo, nLayers, FIELD_LAPLACE, VoxelSize, laplace_iter);
		break;
	default:
		std::cout<<"Unknown method"<<std::endl;
	}


	//generate elements
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkShortArray> AHA17 = vtkSmartPointer<vtkShortArray>::New();
	vtkSmartPointer<vtkShortArray> celllayerid = vtkSmartPointer<vtkShortArray>::New();
	vtkSmartPointer<vtkShortArray> pointlayerid = vtkSmartPointer<vtkShortArray>::New();
	vtkSmartPointer<vtkShortArray> pointSubpartId = vtkSmartPointer<vtkShortArray>::New();
	AHA17->SetName("AHA17");
	celllayerid->SetName("CellLayerID");
	pointlayerid->SetName("PointLayerID");
	pointSubpartId->SetName("subpartID");

	bool aha17_present = endo->GetCellData()->HasArray("AHA17");
	vtkShortArray* aha17_old = NULL;
	if(aha17_present) 
		aha17_old = dynamic_cast<vtkShortArray*>(endo->GetCellData()->GetArray("AHA17"));
	else
		std::cout<<"AHA17 cell array not found. Ignoring."<<std::endl;



	std::cout<<"Generating wedge elements"<<std::endl;
	for(int i=0; i<endo->GetNumberOfCells(); i++ )
	{
		vtkCell *endocell = endo->GetCell(i);
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
				wedge->GetPointIds()->SetId(j, layer*nPointPerLayer+ptIds->GetId(j));
				wedge->GetPointIds()->SetId(j+3, (layer+1)*nPointPerLayer+ptIds->GetId(j));
			}
			cells->InsertNextCell(wedge);
			if(aha17_present) AHA17->InsertNextValue( aha17_old->GetValue(i) );
			celllayerid->InsertNextValue(layer);
		}

	}
	
	endo->SetPolys(cells);
	

	std::cout<<"Saving"<<std::endl;
	vtkSmartPointer<vtkUnstructuredGrid> output = vtkSmartPointer<vtkUnstructuredGrid>::New();
//	output->Allocate(cells->GetNumberOfCells());
	output->SetPoints(layers);

	//generate pointlayerid
	pointlayerid->SetNumberOfValues(output->GetNumberOfPoints());
	for(int i=0; i<nPointPerLayer; i++)
	{
		for (int j=0; j<nLayers; j++)
		{
			pointlayerid->SetValue(j*nPointPerLayer+i,j);
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
	if(aha17_present) output->GetCellData()->AddArray(AHA17);
	output->GetCellData()->AddArray(celllayerid);
	output->GetPointData()->AddArray(pointlayerid);
	output->GetPointData()->AddArray(pointSubpartId);

//	blShapeUtils::ShapeUtils::SaveShapeToFile(output,outshape_filename);
	
	if( points_per_level>0 )
	{
		std::cout<<"Generating coordinate system with points per level = "<<points_per_level<<std::endl;
		GenerateLocalCoordinateRadial(output, nLayers);
		GenerateLocalCoordinateCircLongit(output, nLayers, points_per_level ); //add here projection onto perpendicular plane to the Radial vector
	}
	
	vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetFileName(outshape_filename);
	writer->SetInputData(output);
	writer->SetFileTypeToBinary();
	writer->Update();

	return 0;
}


void GenerateLocalCoordinateRadial(vtkUnstructuredGrid* volmesh, int nLayers)
{
	vtkSmartPointer<vtkFloatArray> radial = vtkSmartPointer<vtkFloatArray>::New();

	const int nPoints = volmesh->GetNumberOfPoints();

	radial->SetNumberOfComponents(3);
	radial->SetNumberOfTuples(nPoints);
	radial->SetName("RadialAxis");

	/*---------------------------------------------------
	 * 
	 *      Generate circumferential vectors
	 * 
	 */
	
	const int nPointsPerLayer = nPoints/nLayers;

	//layers go from endo to epi
	//do all the layers except the last one (epi)
	vnl_vector<double> v(3);
	vnl_vector<double> p1(3);
	vnl_vector<double> p2(3);
	
	for( int nLayer = 0; nLayer < nLayers-1; nLayer++)
	{
		for( int ind =0; ind<nPointsPerLayer; ind++ )
		{
			const int pt1Index = ind + nLayer*nPointsPerLayer;
			const int pt2Index = ind + (nLayer+1)*nPointsPerLayer;

			volmesh->GetPoint(pt1Index, p1.data_block());
			volmesh->GetPoint(pt2Index, p2.data_block());
			v = p2-p1;
			v.normalize();

			radial->SetTuple(pt1Index, v.data_block());				
		}
	}

	//make special treatment for the last layer
	for( int ind =0; ind<nPointsPerLayer; ind++ )
	{
		const int pt1Index = ind + (nLayers-1)*nPointsPerLayer;
		const int pt2Index = ind + (nLayers-2)*nPointsPerLayer;

		volmesh->GetPoint(pt1Index, p1.data_block());
		volmesh->GetPoint(pt2Index, p2.data_block());
		v = p1-p2;
		v.normalize();

		radial->SetTuple(pt1Index, v.data_block());				
	}

	volmesh->GetPointData()->AddArray( radial );

}


//nPointsPerLevel - number of points on one Z-level in one layer (one circumference)
void GenerateLocalCoordinateCircLongit(vtkUnstructuredGrid* volmesh, int nLayers, int nPointsPerLevel)
{
	vtkSmartPointer<vtkFloatArray> circ = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkFloatArray> longit = vtkSmartPointer<vtkFloatArray>::New();

	const int nPoints = volmesh->GetNumberOfPoints();

	circ->SetNumberOfComponents(3);
	circ->SetNumberOfTuples(nPoints);
	circ->SetName("CircumferentialAxis");


	longit->SetNumberOfComponents(3);
	longit->SetNumberOfTuples(nPoints);
	longit->SetName("LongitudinalAxis");


	vtkFloatArray *radial = dynamic_cast<vtkFloatArray*>(volmesh->GetPointData()->GetArray("RadialAxis"));


	/*---------------------------------------------------
	 * 
	 *      Generate circumferential vectors
	 * 
	 */
	
	const int nPointsPerLayer = nPoints/nLayers;
	const int nLevelsPerLayer = (nPointsPerLayer-1)/nPointsPerLevel; //0s level is just one point
	
	//std::cout<<"Points per layer: "<<nPointsPerLayer<<std::endl;
	//std::cout<<"Levels per layer: "<<nLevelsPerLayer<<std::endl;
	
	//extract the coordinates of the circumference
	vnl_matrix<double> coords(nPointsPerLevel+2,3); //pad with 2 more values for cyclicity
	vnl_vector<double> v(3);
	vnl_vector<double> p1(3);
	vnl_vector<double> p2(3);
	
	for( int nLayer = 0; nLayer < nLayers; nLayer++)
	{
		//add apex, it's exceptional
		v.fill(0);
		const int ptZeroIndex = 0 + nLayer*nPointsPerLayer;
		circ->SetTuple(ptZeroIndex, v.data_block());				
		longit->SetTuple(ptZeroIndex, v.data_block());				


		for( int nLevel=0; nLevel < nLevelsPerLayer; nLevel++)
		{
			double pt[3];
			
			const int ptFirstIndex = 0 + nLevel*nPointsPerLevel + nLayer*nPointsPerLayer  +1; //+1 to account for apex
			const int ptLastIndex = nPointsPerLevel-1 + nLevel*nPointsPerLevel + nLayer*nPointsPerLayer;
			
			volmesh->GetPoint(ptLastIndex, pt);
			coords.set_row( 0, pt );
			volmesh->GetPoint(ptFirstIndex, pt);
			coords.set_row( coords.rows()-1, pt );
			
			for( int ind=0; ind<nPointsPerLevel; ind++ )
			{
				const int ptIndex = ind + nLevel*nPointsPerLevel + nLayer*nPointsPerLayer +1;//+1 to account for apex
				volmesh->GetPoint(ptIndex, pt);
				coords.set_row( ind+1, pt );
			}
			
			//generate the vectors as central differences
			for( unsigned int i=1; i<coords.rows()-1; i++ )
			{
				p1 = coords.get_row(i-1);
				p2 = coords.get_row(i+1);
				v = p2-p1;
				v.normalize();
				
				

				const int ptIndex = i-1 + nLevel*nPointsPerLevel + nLayer*nPointsPerLayer +1;	//+1 to account for apex			

				vnl_vector<double> v_radial(3);
				radial->GetTuple( ptIndex, v_radial.data_block() );

				//make sure it is orthogonal to radial direction
				v = v- dot_product(v,v_radial)*v_radial;
				v.normalize();

				circ->SetTuple(ptIndex, v.data_block());				
				
				vnl_vector<double> v_longit = vnl_cross_3d(v, v_radial);
				v_longit.normalize();
				longit->SetTuple(ptIndex, v_longit.data_block());				
			}
					
		}
	}
	
	
	 volmesh->GetPointData()->AddArray( circ );
	 volmesh->GetPointData()->AddArray( longit );
}






void GenerateLayersAlongNormals(vtkPoints* layers, vtkPolyData* epi, vtkPolyData* endo, int nLayers)
{
	const int nPointsPerLayer = endo->GetNumberOfPoints();

	vnl_matrix<double> endo_points(endo->GetNumberOfPoints(),3);
	vnl_matrix<double> located_epi_points(endo->GetNumberOfPoints(),3);

	std::cout<<"Generate layers along normals"<<std::endl;
	vtkSmartPointer<vtkCellLocator> loc = vtkSmartPointer<vtkCellLocator>::New();
	loc->SetDataSet(epi);
	loc->BuildLocator();

	//CommonTools::SavePolydata(epi,"epi111.vtk");
	//CommonTools::SavePolydata(endo,"endo1111.vtk");

	vtkSmartPointer<vtkPolyDataNormals> endo_normalgen = vtkSmartPointer<vtkPolyDataNormals>::New();
	endo_normalgen->SetInputData(endo);
	endo_normalgen->SplittingOff();
	endo_normalgen->Update();
	
	//vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	//pts->SetNumberOfPoints(endo->GetNumberOfPoints());

	//locate the points on epi lying on the intersection with the normal from endo
	for (int i=0; i<endo->GetNumberOfPoints(); i++)
	{
		vnl_vector<double> normal(3);
		endo_normalgen->GetOutput()->GetPointData()->GetNormals()->GetTuple(i, normal.begin());
		normal.normalize();

		vnl_vector<double> pt1(3);
		vnl_vector<double> pt2(3);
		endo->GetPoint(i,pt1.begin());
		pt2 = pt1 - 70.0*normal;

		endo_points.set_row(i, pt1);

		int subId;
		double t;
		vnl_vector<double> x(3);
		double pcoords[3];
		loc->IntersectWithLine(pt1.begin(),pt2.begin(),0,t,x.begin(),pcoords,subId);
		located_epi_points.set_row(i, x);
		
		//pts->SetPoint(i, x.begin() );
	};

	//CommonTools::SavePoints( pts, "points.vtk");

	//generate layers
	const double factor = 1.0/static_cast<double>(nLayers-1);
	for ( unsigned int i=0; i<endo_points.rows(); i++)
	{
		vnl_vector<double> pt1(3);
		vnl_vector<double> pt2(3);
		pt1 = endo_points.get_row(i);
		pt2 = located_epi_points.get_row(i);

		for(int j=0;j<nLayers; j++)
		{
			vnl_vector<double> pt = pt1+(pt2-pt1)*static_cast<double>(j)*factor;
			layers->SetPoint(i+j*nPointsPerLayer,pt.begin());
		}
	}

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


	std::cout<<"Dilating epi mask..."<<std::endl;
	vtkSmartPointer<vtkImageContinuousDilate3D> dilate = vtkSmartPointer<vtkImageContinuousDilate3D>::New();
	dilate->SetInputData(mask_epi);
	dilate->SetKernelSize(5,5,5);
	dilate->Update();

	std::cout<<"Computing Laplace field"<<std::endl;
	//combine masks
	vtkSmartPointer<vtkImageData> combinedMask = vtkSmartPointer<vtkImageData>::New();
	combinedMask->DeepCopy(mask_epi);

	char* mask_endo_data = static_cast<char*>(erode->GetOutput()->GetScalarPointer());
	char* mask_epi_data = static_cast<char*>(dilate->GetOutput()->GetScalarPointer());
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
	//inimage->SetScalarTypeToChar();
	//nimage->Update();
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
	filt->SetInputData( endo_closed );
	filt->Update();

	vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
	stencil->SetStencilData(filt->GetOutput());
	stencil->SetBackgroundValue(bg);
	stencil->SetInputData(inimage);
	stencil->SetOutput(res_image);
	stencil->Update();
}










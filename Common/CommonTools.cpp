/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <vtkCleanPolyData.h>
#include <vtkThreshold.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataSetWriter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageShrink3D.h>
#include <vtkImageData.h>
#include <vtkShortArray.h>
#include <vtkPolyData.h>
#include <vtkDataSet.h>
#include <vtkPoints.h>
#include <vtkStripper.h>
#include <vtkCutter.h>
#include <vtkPointLocator.h>
#include <vtkPlane.h>
#include <vtkStringArray.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkDataSetReader.h>
#include <vtkFeatureEdges.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPLYReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkIVWriter.h>
#include <vtkCellLocator.h>
#include <vtkUnstructuredGridWriter.h>

	#include "CommonTools.h"
	#include <stdio.h>
	#include <stdlib.h>
	#include <fstream>
	#include <stdexcept>
	#include <string>


#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>


	//---------------------------------------------------
	// OPERATIONS
	//---------------------------------------------------





void CommonTools::GenerateDecimationScript( const char* filename, int nfaces )
{
	char nfaces_str[100];
	sprintf(nfaces_str,"%d", nfaces);

	std::ofstream file(filename);
	file<<"<!DOCTYPE FilterScript>"<<std::endl;
	file<<"<FilterScript>"<<std::endl;
	file<<"<filter name=\"Quadric Edge Collapse Decimation\">"<<std::endl;
	file<<"<Param type=\"RichInt\" value=\""<<nfaces_str<<"\" name=\"TargetFaceNum\"/>"<<std::endl;
	file<<"<Param type=\"RichFloat\" value=\"0\" name=\"TargetPerc\"/>"<<std::endl;
	file<<"<Param type=\"RichFloat\" value=\"0.3\" name=\"QualityThr\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"false\" name=\"PreserveBoundary\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"false\" name=\"PreserveNormal\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"false\" name=\"PreserveTopology\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"true\" name=\"OptimalPlacement\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"false\" name=\"PlanarQuadric\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"false\" name=\"QualityWeight\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"true\" name=\"AutoClean\"/>"<<std::endl;
	file<<"<Param type=\"RichBool\" value=\"false\" name=\"Selected\"/>"<<std::endl;
	file<<"</filter>"<<std::endl;
	file<<"</FilterScript>"<<std::endl;
	file.close();
}


void CommonTools::SaveVtkShortArray( const char *filename, vtkShortArray* the_array )
{
	FILE* f = fopen(filename,"wb");

	if(!f)
	{
		std::cout << "Failed opening "<<filename<<std::endl;
		throw std::runtime_error( "Failed opening file" );
	}

	char id[] = "RD10";
	fwrite(id, 4, 1, f );

	int size = the_array->GetSize();
	fwrite(&size, sizeof(size), 1, f);

	std::string name;
	if ( the_array->GetName() )
	{
		name = the_array->GetName();
	}
	fwrite(name.c_str(), name.length()+1, 1, f); //+1 to include eol



	short *scalars = new short[size];

	for( int i=0; i<size; i++ ) 
		scalars[i] = the_array->GetValue(i);

	fwrite(scalars, 1, size*sizeof(scalars[0]), f);


	fclose(f);

	delete[] scalars;
}

void CommonTools::LoadVtkShortArray( const char *filename, vtkShortArray* the_array )
{
	if( strlen(filename) == 0 )
	{
		std::cout << "Failed opening empty file"<<std::endl;
		throw std::runtime_error( "Failed opening file" );
	}

	FILE* f = fopen(filename,"rb");
	if( !f )
	{
		std::cout <<"Failed opening "<<filename<<std::endl;
		throw std::runtime_error( "Failed opening file" );
	}

	char id[]="1234\x0";
	
	int waste;
	
	waste = fread(id, sizeof( id[0] ), 4, f ); //read four symbols
	if( strcmp(id, "RD10")!=0 )
	{
		fclose(f);
		std::cout <<"Region definitions has incorrect version number"<<std::endl;
		throw std::runtime_error( "Region definitions has incorrect version number" );
	}

	int size;
	waste = fread(&size, sizeof(size), 1, f);

	std::string name="";
	char symbol = ' ';
	while( symbol!='\x0' )
	{
		waste = fread(&symbol, sizeof(symbol), 1, f);
		name+=symbol;
	}

	the_array->SetNumberOfValues( size );
	the_array->SetName( name.c_str() );

	short *scalars = new short[size];
	int b = fread(scalars, 1,size*sizeof(scalars[0]),f);
	if( b!=size*static_cast<int>(sizeof(scalars[0])) )
	{
		fclose(f);
		std::cout  << "Number of items read from region definition file is incorrect. "
			<<b<<" out of "<<size<<". Aborting."<<std::endl;
		throw std::runtime_error( "Number of items read from region definition file is incorrect" );
	}

	for( int i=0; i<size; i++) the_array->SetValue(i, scalars[i]);

	fclose(f);

	delete[] scalars;
}


vtkPolyData* CommonTools::GetShapeSubSurface(vtkPolyData * inputShape, double tholdLower, double tholdUpper)
{
	//extract the subpart
	vtkSmartPointer<vtkThreshold> thold = vtkSmartPointer<vtkThreshold>::New();
	thold->SetInputData(inputShape);
	thold->ThresholdBetween(tholdLower, tholdUpper);
	thold->Update();

	//extract surface
	vtkSmartPointer<vtkDataSetSurfaceFilter> sur_filt = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
	sur_filt->SetInputData( (vtkDataObject*)thold->GetOutput() );
	sur_filt->Update();

	vtkSmartPointer<vtkCleanPolyData> pdCleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	pdCleaner -> SetInputData(sur_filt->GetOutput());
	pdCleaner -> SetTolerance (0.0);
	pdCleaner -> SetAbsoluteTolerance (1.0);
	//pdCleaner -> ToleranceIsAbsoluteOn ();
	pdCleaner -> ConvertStripsToPolysOn ();
	pdCleaner -> ConvertPolysToLinesOn ();
	pdCleaner -> PointMergingOn ();
	pdCleaner -> Update();

	///
	vtkPolyData * outputShape = vtkPolyData::New();
	outputShape->DeepCopy(pdCleaner->GetOutput());

	return outputShape;
}	



//------------------------------------------------------------------
vtkPolyData* CommonTools::GetShapeSubSurface(vtkPolyData * inputShape, unsigned int nSubPart)
//------------------------------------------------------------------
{
	return CommonTools::GetShapeSubSurface(inputShape, nSubPart-0.1, nSubPart+0.1);
}	











void CommonTools::SavePolydata( vtkPolyData* poly, const char* filename, bool binary )
{
	//vtkSTLWriter * writer = vtkSTLWriter::New();
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	if ( binary ) writer->SetFileTypeToBinary();
	writer -> SetInputData(poly);
	writer -> SetFileName (filename);
	writer -> Write();
}




void CommonTools::SaveImage( vtkDataSet* image, const char* filename )
{
	vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
	writer->SetFileTypeToBinary();
	writer->SetInputData( (vtkDataObject*)image );
	writer->SetFileName( filename );
	writer->Write();
}

void CommonTools::SaveUnstructuredGrid( vtkUnstructuredGrid* grid, const char* filename )
{
	vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetFileName(filename);
	writer->SetInputData(grid);
	writer->SetFileTypeToBinary();
	writer->Update();	
}




vtkImageData* CommonTools::LoadImage( const char* filename )
{
	vtkSmartPointer<vtkDataSetReader> imread = vtkSmartPointer<vtkDataSetReader>::New();
	imread->SetFileName( filename );
	imread->Update();
	
	vtkImageData *result = vtkImageData::New();
	result->DeepCopy( dynamic_cast<vtkImageData*>(imread->GetOutput()) );
	return result;
}






vtkPolyData* CommonTools::CloseSurface( vtkPolyData* shape )
{
	vtkSmartPointer<vtkFeatureEdges> fe = vtkSmartPointer<vtkFeatureEdges>::New();
	fe->SetInputData( shape );
	fe->BoundaryEdgesOn();
	fe->NonManifoldEdgesOff();
	fe->FeatureEdgesOff();
	fe->ManifoldEdgesOff();
	fe->Update();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> connect = 
				vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	connect->SetInputData(fe->GetOutput());
	connect->Update();

	const int ncontours = connect->GetNumberOfExtractedRegions();
		
	vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();
	append->AddInputData (shape);


	for( int i =0; i<ncontours ; i++ )
	{
		connect->AddSpecifiedRegion(i);
		connect->SetExtractionModeToSpecifiedRegions();
		connect->Update();
		vtkPolyData *edges = connect->GetOutput();
		vtkPolyData* cover = GenerateHoleCover(edges);

		//CommonTools::SavePolydata(edges,"edges.vtk");
		//CommonTools::SavePolydata(cover,"cover.vtk");

		append->AddInputData (cover);

		cover->Delete();
		connect->DeleteSpecifiedRegion(i);
	}

	append->Update();


	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner->SetInputData(append->GetOutput());
	cleaner->Update();

	vtkPolyData *result = vtkPolyData::New();
	result->DeepCopy(cleaner->GetOutput());

	return result;
}


vtkPolyData* CommonTools::GenerateHoleCover(vtkPolyData* edges)
{
	// We'll create the building blocks of polydata including data attributes.
	vtkPolyData* cover = vtkPolyData::New();
	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkCleanPolyData> sur_filt = vtkSmartPointer<vtkCleanPolyData>::New();
	sur_filt->SetInputData( edges );
	sur_filt->Update();

	points->DeepCopy(sur_filt->GetOutput()->GetPoints());

	//add centroid
	double centr[3] = {0,0,0};
	for(int i=0; i<points->GetNumberOfPoints(); i++)
	{
		for(int j=0; j<3; j++)
			centr[j] += points->GetPoint(i)[j];
	}
	for(int j=0; j<3; j++)
		centr[j] /= points->GetNumberOfPoints();

	vtkIdType cnt_pt = points->InsertNextPoint(centr);

	//add cells
	for (int i=0; i<sur_filt->GetOutput()->GetNumberOfCells(); i++)
	{
		vtkCell *cell = sur_filt->GetOutput()->GetCell(i);

		vtkIdType pts3[3];
		pts3[0] = cell->GetPointId(0);
		pts3[1] = cell->GetPointId(1);
		pts3[2] = cnt_pt;
		polys->InsertNextCell(3,pts3);
	}

	// We now assign the pieces to the vtkPolyData.
	cover->SetPoints(points);
	cover->SetPolys(polys);

	return cover;
}





void CommonTools::SavePoints( vtkPoints* pts, const char* filename )
{
	vtkPolyData* ply = CommonTools::Points2Polydata( pts );
	CommonTools::SavePolydata( ply, filename );
	ply->Delete();
}




bool CommonTools::FileExists( const char* filename, bool no_exception )
{
	bool result = true;

	FILE *fid = fopen(filename,"r");
	if( fid == NULL )
	{
		std::cerr<<"File not found: "<<filename;

		char error[500];
		sprintf(error,"File not found: %s",filename);

		if( !no_exception ) throw std::runtime_error( error );

		result = false;
	}
	else fclose(fid);

	return result;
}



//stores in "list"
void CommonTools::ReadFilelist( const char* file, std::vector<std::string>& list, bool check_existence )
{
	list.clear();

	FileExists( file );

	std::ifstream the_file( file );
	while( !the_file.eof() )
	{
		const int maxlen = 2000;
		char line[maxlen];
		the_file.getline( line, maxlen );
		
		if( strcmp(line,"")!=0 )
		{
			if( line[0]!=' ' )
			{
				if(check_existence) FileExists( line );
				list.push_back( line );
			}
		}
	}
}



vtkPolyData* CommonTools::Points2Polydata( vtkPoints* points, const double* scalars)
{
	vtkPolyData *polyData = vtkPolyData::New();

	vtkCellArray *candidateCells = vtkCellArray::New();
	vtkFloatArray *candidateScalars = vtkFloatArray::New();

	bool manual_scalars = (scalars==NULL);

	for (vtkIdType j = 0; j < points->GetNumberOfPoints(); j++) {
		candidateCells -> InsertNextCell(VTK_VERTEX,&j);
		const float candidateScalar = manual_scalars ? j : scalars[j];
		candidateScalars -> InsertNextTuple(&candidateScalar);
	}


	polyData -> SetPoints(points);
	polyData -> SetVerts(candidateCells);
	candidateCells->Delete();
	polyData -> GetPointData() -> SetScalars(candidateScalars);
	candidateScalars->Delete();
	polyData -> Modified();

	return  polyData;
}



vtkPolyData* CommonTools::Points2Polydata( vtkPoints* points, double scalar )
{
	vnl_vector<double> scalars(points->GetNumberOfPoints(),scalar);

	return  Points2Polydata( points, scalars.begin() );
}



void CommonTools::ExportPolyDataPoints( vtkPolyData* shape, vnl_matrix<double>& points )
{
	const int ndim = 3; //number of dimensions

	points.set_size( shape->GetNumberOfPoints(), ndim );
	for( int vert = 0; vert<shape->GetNumberOfPoints(); vert++ )
	{
		vnl_vector<double> point(ndim);
		shape->GetPoint(vert, point.begin());

		points.set_row(vert, point);
	}
}


void CommonTools::ImportPolyDataPoints( vtkPolyData* shape, vnl_matrix<double>& points )
{
	for(unsigned int vert = 0; vert<points.rows(); vert++ )
	{
		vnl_vector<double> point = points.get_row(vert);
		shape->GetPoints()->SetPoint(vert, point.begin());
	}

//	shape->Update();
}



void CommonTools::ExportPolyDataPoints( vtkPolyData* shape, vnl_vector<double>& points )
{
	const int ndim = 3; //number of dimensions

	points.set_size( shape->GetNumberOfPoints()*ndim );
	int c=0;
	for( int vert = 0; vert<shape->GetNumberOfPoints(); vert++ )
	{
		vnl_vector<double> point(ndim);
		shape->GetPoint(vert, point.begin());

		for( size_t k=0; k<point.size(); k++ ) 
			points[c++] = point[k];
	}
}


void CommonTools::ImportPolyDataPoints( vtkPolyData* shape, vnl_vector<double>& points )
{
	const int ndim = 3; //number of dimensions

	int c=0;
	for( int vert = 0; vert<shape->GetNumberOfPoints(); vert++ )
	{
		vnl_vector<double> point(ndim);

		for( size_t k=0; k<point.size(); k++ ) 
			point[k] = points[c++];

		shape->GetPoints()->SetPoint(vert, point.begin());
	}

//	shape->Update();
}


void CommonTools::ScaleShape(vtkPolyData *shapein, vtkPolyData *shapeout, float scale, 
														bool centerAfterScale /* = false*/)
{
	vtkSmartPointer<vtkTransform> shape_tform = vtkSmartPointer<vtkTransform>::New();
	//vtkSmartPointer<vtkTransformPolyDataFilter> shape_tformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	vtkSmartPointer<vtkTransformFilter> shape_tformer = vtkSmartPointer<vtkTransformFilter>::New();

	shape_tform->Scale( scale, scale, scale );

	if (centerAfterScale == true)
	{
		double* origBounds = shapein->GetBounds();	
		double origCenter[3];
		origCenter[0] = (origBounds[1]+origBounds[0])/2;
		origCenter[1] = (origBounds[3]+origBounds[2])/2;
		origCenter[2] = (origBounds[5]+origBounds[4])/2;
		
		shape_tform->Translate(origCenter[0]*(1-scale)/scale, 
								origCenter[1]*(1-scale)/scale, 
								origCenter[2]*(1-scale)/scale);
	}
	shape_tformer->SetInputData(shapein);
	shape_tformer->SetTransform(shape_tform);
	//shape_tformer->SetOutput( shapeout );
	shape_tformer->Update();
	shapeout->DeepCopy(shape_tformer->GetOutput());
	vtkSmartPointer<vtkDataArray> data = shapein->GetPointData()->GetVectors();
	if( data!= NULL ) 
	{
		shapeout->GetPointData()->RemoveArray(data->GetName());
		shapeout->GetPointData()->SetVectors(data);
	}
}

void CommonTools::ShrinkImage(vtkDataSet *imagein, vtkDataSet *imageout, int factor)
{
	vtkSmartPointer<vtkImageGaussianSmooth> smoother = vtkSmartPointer<vtkImageGaussianSmooth>::New();
	vtkSmartPointer<vtkImageShrink3D> resampler = vtkSmartPointer<vtkImageShrink3D>::New();

	smoother->SetInputData( imagein );
	smoother->SetStandardDeviation( 0.5*(float)factor ); //like in ITK
	smoother->Update();

	resampler->SetInputData( (vtkDataObject*)smoother->GetOutput() );
	resampler->SetShrinkFactors( factor, factor, factor );
	resampler->SetOutput(imageout);
	resampler->Update();

	//shrinker doubles the spacing to produce the image of the same physical size
	//this must be fixed

	double spacing[3];
	((vtkImageData*)imageout)->GetSpacing( spacing );
	for( int kk = 0; kk<3; kk++ ) spacing[kk] /= factor;
	((vtkImageData*)imageout)->SetSpacing( spacing );
//	imageout->Update();

	const int *dimensions = ((vtkImageData*)imageout)->GetDimensions();

	const double *old_origin = ((vtkImageData*)imagein)->GetOrigin();

	//this should change depending on modality. In this case it's needed for CT
	//this is just a quick solution.
	if( old_origin[0]!=0 || old_origin[1]!=0 || old_origin[2]!=0 )
	{
		std::cout<<"old_origin:"<<old_origin[0]<<" "<<old_origin[1]<<" "<<old_origin[2]<<std::endl;
		double origin[3];
		origin[0] = -(dimensions[0]*spacing[0]/2) + spacing[0]/2;
		origin[1] = -(dimensions[1]*spacing[1]/2) + spacing[1]/2;
		origin[2] = -(dimensions[2]*spacing[2]/2) + spacing[2]/2;
		((vtkImageData*)imageout)->SetOrigin( origin[0],origin[1],origin[2] );
		std::cout<<"origin:"<<origin[0]<<" "<<origin[1]<<" "<<origin[2]<<std::endl;
	}



}

void CommonTools::ScaleVolume(vtkUnstructuredGrid* volumein, vtkUnstructuredGrid* volumeout, float scale)
{
	vtkSmartPointer<vtkTransform> volume_tform = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkTransformFilter> volume_tformer = vtkSmartPointer<vtkTransformFilter>::New();
	volume_tform->Scale( scale, scale, scale );

	volume_tformer->SetInputData(volumein);
	volume_tformer->SetTransform(volume_tform);
	volume_tformer->Update();
	volumeout->DeepCopy(volume_tformer->GetOutput());
	vtkSmartPointer<vtkDataArray> data = volumein->GetPointData()->GetVectors();
	if( data!= NULL ) 
		{
			volumeout->GetPointData()->RemoveArray(data->GetName());
			volumeout->GetPointData()->SetVectors(data);
		}
}






//------------------------------------------------------------------
void CommonTools::GetP2S(vtkPolyData * manualPt, vtkPolyData *segmentedPt, 
						double& mean, double& std_dev, double& max, 
						double& last, bool b_array)
//------------------------------------------------------------------
{
	
	unsigned int n_manual_points = manualPt -> GetNumberOfPoints();

	vtkSmartPointer<vtkFloatArray> error_array;
	if (b_array) 
	{
		error_array = vtkSmartPointer<vtkFloatArray>::New();
		error_array -> SetName("p2s");
		error_array -> SetNumberOfValues(n_manual_points);
		manualPt -> GetPointData()-> AddArray(error_array);
	}

	// fill the CellLocator with the fitted shape
	vtkSmartPointer<vtkCellLocator> cell_locator = vtkSmartPointer<vtkCellLocator>::New();
	cell_locator->SetDataSet(segmentedPt);
	cell_locator->BuildLocator();
	
	vnl_vector<float> stdaux(n_manual_points);
	
	for (size_t i=0; i<n_manual_points; i++) 
	{
		double *p1 = manualPt -> GetPoint(i);
		double p2[] = {0.0, 0.0, 0.0};
		int subId;
		vtkIdType cellId;
		double dist;
		cell_locator->FindClosestPoint(p1, p2, cellId, subId, dist);
		stdaux[i] = sqrt(dist);
		error_array->SetValue(i,stdaux[i]);
	}

	max = stdaux.max_value();
	//min = stdaux.min_value();
	mean = stdaux.mean();
	
	double sum=0;
	double sum2=0;

	for (size_t i = 0; i < n_manual_points; i++)
	{
		sum2=sum2+(stdaux[i]-mean)*(stdaux[i]-mean);
		sum=sum+(stdaux[i]-mean);
	}

	std_dev = sqrt((sum2-sum/n_manual_points)/(n_manual_points-1));

	last = stdaux[n_manual_points-1];
}

//------------------------------------------------------------------
void CommonTools::GetP2P(vtkPolyData * manualPt, vtkPolyData *segmentedPt, 
						double& mean, double& std_dev, double& max, 
						double& last, bool b_array)
//------------------------------------------------------------------
{
	unsigned int n_manual_points = manualPt -> GetNumberOfPoints();

	vtkSmartPointer<vtkFloatArray> error_array;

	if (b_array) 
	{
		error_array = vtkSmartPointer<vtkFloatArray>::New();
		error_array -> SetName("p2p");
		error_array -> SetNumberOfValues(n_manual_points);
		manualPt -> GetPointData()-> AddArray(error_array);
	}
	
	vnl_vector<float> stdaux(n_manual_points);
	vnl_vector<double> p1(3);
	vnl_vector<double> p2(3);
	vnl_vector<double> dist(3);
	
	for (size_t i=0; i<n_manual_points; i++) 
	{
		manualPt -> GetPoint(i,p1.data_block());
		segmentedPt -> GetPoint(i,p2.data_block());
		dist = p1 - p2;
		stdaux[i] = dist.magnitude();
		error_array->SetValue(i,stdaux[i]);
	}

	max = stdaux.max_value();
	mean = stdaux.mean();
	
	double sum=0;
	double sum2=0;

	for (size_t i = 0; i < n_manual_points; i++)
	{
		sum2=sum2+(stdaux[i]-mean)*(stdaux[i]-mean);
		sum=sum+(stdaux[i]-mean);
	}

	std_dev = sqrt((sum2-sum/n_manual_points)/(n_manual_points-1));

	last = stdaux[n_manual_points-1];
}


vtkPolyData* CommonTools::GetP2S(vtkPolyData *shapePt1, vtkPolyData *shapePt2, std::vector< vnl_vector<double> >& distances)
{
	const double* scalar_range = shapePt1->GetScalarRange();


	vtkSmartPointer<vtkAppendPolyData> appender = vtkSmartPointer<vtkAppendPolyData>::New();

	for( int region = floor(scalar_range[0]); region<=floor(scalar_range[1]); region++ )
	{
		vnl_vector<double> dist;

		vtkSmartPointer<vtkPolyData> shape1 = vtkSmartPointer<vtkPolyData>::Take(
			GetShapeSubSurface(shapePt1,region) );
		vtkSmartPointer<vtkPolyData> shape2 = vtkSmartPointer<vtkPolyData>::Take(
			GetShapeSubSurface(shapePt2,region) );

		vtkSmartPointer<vtkPolyData> shape = vtkSmartPointer<vtkPolyData>::Take(
			GetP2S( shape1, shape2, dist ) );

		if( shape.GetPointer()!=NULL )
		{
			vnl_vector<double> regionID(1);
			regionID[0] = region;
			appender->AddInputData (shape);
			distances.push_back( regionID );
			distances.push_back( dist );
		}
	}

	appender->Update();

	vtkPolyData* result = vtkPolyData::New();
	result->DeepCopy( appender->GetOutput() );

	return result;
}


vtkPolyData* CommonTools::GetP2S(vtkPolyData *shapePt1, vtkPolyData *shapePt2, vnl_vector<double>& distances)
{
	if( shapePt1->GetNumberOfPoints()==0 || shapePt2->GetNumberOfPoints()==0 )
	{
		return NULL;
	}


	vtkPolyData* shape1 = vtkPolyData::New();
	shape1->DeepCopy( shapePt1 );



	vtkSmartPointer<vtkFloatArray> scalars1 = vtkSmartPointer<vtkFloatArray>::New();
	scalars1->SetName("P2SError");
	scalars1->SetNumberOfValues( shape1->GetNumberOfPoints() );

	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(shapePt2);
	cellLocator->BuildLocator();

	const unsigned int nPoints1 = shape1 -> GetNumberOfPoints();

	vtkSmartPointer<vtkPoints> intersectionPoints2 = vtkSmartPointer<vtkPoints>::New();

	distances.set_size(nPoints1);


	for (unsigned int i=0; i<nPoints1; i++) {
		double *p1 = shape1 -> GetPoint(i);
		double p2[] = {0.0, 0.0, 0.0};
		int subId;
		vtkIdType cellId;
		double dist;
		cellLocator->FindClosestPoint(p1, p2, cellId, subId, dist);
		dist = sqrt(dist);
		distances[i] = dist;
		scalars1 -> SetValue(i,dist);
		intersectionPoints2->InsertNextPoint(p2);
	}


	shape1->GetPointData()->AddArray( scalars1 );


	return shape1;
}


void CommonTools::GetStatistics( vnl_vector<double>& v, double& mean, double& std )
{
    mean = v.mean();
    std = 0;
    for(int i=0; i<v.size(); i++)
    {
        std = std + (v[i]-mean)*(v[i]-mean);
    }
    std /= (v.size()-1);
    std = sqrt(std);
}

void CommonTools::GetS2S(vtkPolyData *shapePt1, vtkPolyData *shapePt2, 
        vtkPolyData* outShape1, 
        vtkPolyData* outShape2,
         double& mean, double& stdev, double& maximum)
{
    vnl_vector<double> errors1;
    vnl_vector<double> errors2;
    
    vtkSmartPointer<vtkPolyData> outShape11 = vtkSmartPointer<vtkPolyData>::Take( GetP2S( shapePt1, shapePt2, errors1 ) );
    vtkSmartPointer<vtkPolyData> outShape21 = vtkSmartPointer<vtkPolyData>::Take( GetP2S( shapePt2, shapePt1, errors2 ) );

    
    outShape1->DeepCopy(outShape11.GetPointer()); 
    outShape2->DeepCopy(outShape21.GetPointer());

    vnl_vector<double> errors(errors1.size()+errors2.size());
    
    for(int i=0; i<errors1.size(); i++)
        errors[i]=errors1[i];
    
    for(int i=0; i<errors2.size(); i++)
        errors[i+errors1.size()]=errors2[i];
    
    GetStatistics(errors, mean, stdev);
    maximum = errors.max_value();
    
    /*
    double mean1, mean2, std1, std2;

    GetStatistics(errors1, mean1, std1);
    GetStatistics(errors2, mean2, std2);
            
    const double N1 = errors1.size();
    const double N2 = errors2.size();
    mean = (N1*mean1 + N2*mean2) / (N1+N2);
    
    const double var1 = std1*std1;
    const double var2 = std2*std2;
    
    stdev = sqrt( ((N1-1)*var1 + (N2-1)*var2) / (N1-1 + N2-1) );

    double max1 = errors1.max_value();
    double max2 = errors2.max_value();
    
    
    //take the maximum of the two maxima
    maximum = max1;
    if ( max2>maximum )
        maximum = max2;
     */
}


void CommonTools::GetS2S(vtkPolyData *shapePt1, vtkPolyData *shapePt2, std::vector< vnl_vector<double> >& distances )
{
	std::vector< vnl_vector<double> > errors;
	std::vector< vnl_vector<double> > errors1;
	vtkPolyData* deformation = GetP2S( shapePt1, shapePt2, errors );
	vtkPolyData* deformation1 = GetP2S( shapePt2, shapePt1, errors1 );
	deformation1->Delete();
	deformation->Delete();

	distances.resize(errors.size());

	for( unsigned int i=0; i<errors.size(); i+=2 )
	{
		distances[i] = errors[i];
		distances[i+1] = (errors[i+1] + errors1[i+1])/2.0; //get symmetric S2S error
	}

}









/* Input values (c) Ruben Cardenes + Constantine Butakoff
3 Outside domain
1 Exterior boundary //not used
0 Interior boundary //not used
2 Inside domain 
*/ 
//-------------------------------------------------------------------------------
int CommonTools::laplace3D_voxelsize( vtkImageData* inputImage, vtkImageData* outputImage, int iterations/*=100*/ )
//-------------------------------------------------------------------------------
{
	//int *dim = inputImage->GetDimensions();
	//const int max1=dim[0];
	//const int max2=dim[1];
	//const int max3=dim[2];
	double hx, hy, hz;
	inputImage->GetSpacing(hx, hy, hz);


	outputImage->DeepCopy(inputImage);
//	outputImage->SetOrigin(inputImage->GetOrigin());
//	outputImage->SetSpacing(inputImage->GetSpacing());
//	outputImage->SetDimensions(inputImage->GetDimensions());
//	outputImage->SetScalarType(VTK_FLOAT,NULL);
	//outputImage->SetNumberOfScalarComponents(1);
	//outputImage->Update();
	outputImage->AllocateScalars(VTK_FLOAT,1);


	char* input = static_cast<char*>(inputImage->GetScalarPointer());
	float* output = static_cast<float*>(outputImage->GetScalarPointer());

	int k,l;
//	int sum = 0;
	/* Initialize domain, inside=0, and boundaries values*/
	for (k=0;k<inputImage->GetNumberOfPoints();k++) 
	{
		if (input[k] == 2) {
			//output[k][i][j] = 0;
			output[k] = 0;  //calculation region
		} else if (input[k] == 3) { 
			//output[k][i][j] = -1;
			output[k] = -1; //exterior
		} else {
			//output[k][i][j] = input[sum];
			output[k] = input[k]; //interior, must be = 1... INterior and exterior can be interchanged
		}
	}

	/* Solve Laplacian */
	//compute increments for differentiation
	int ijk[3]={1,1,1};
	vtkIdType ptId= outputImage->ComputePointId(ijk);
	ijk[0]=2;
	vtkIdType ptId1= outputImage->ComputePointId(ijk);
	const register int incXp = ptId1-ptId;
	ijk[0]=1;
	ijk[1]=2;
	ptId1= outputImage->ComputePointId(ijk);
	const register int incYp = ptId1-ptId;
	ijk[1]=1;
	ijk[2]=2;
	ptId1= outputImage->ComputePointId(ijk);
	const register int incZp = ptId1-ptId;

	//precalculate constants
	const float hx2 = hx*hx;
	const float hy2 = hy*hy;
	const float hz2 = hz*hz;
	const float kkk = 0.5 *(hx2*hy2*hz2)/(hx2*hy2 + hy2*hz2 + hx2*hz2);
	
	//const float overhx2 = kkk/hx2;
	//const float overhy2 = kkk/hy2;
	const float overhz2 = kkk/hz2;

	for (l=0;l<iterations;l++) 
	{
		if(l%100==0)
		{
			printf("Iteration %05d/%05d\r",l,iterations);
			fflush(stdout);
		}


		register vtkIdType ptId=inputImage->GetNumberOfPoints()-1;

		while(ptId--) 
		{
			if (input[ptId] == 2) //make sure it is not on the edge 
			{
				//make use of uniform voxel size
				output[ptId] = ( 
					(output[ptId-incXp] + output[ptId+incXp]) + 
					(output[ptId-incYp] + output[ptId+incYp]) + 
					(output[ptId-incZp] + output[ptId+incZp])
					) *overhz2;

			}
		}
	}
	std::cout<<std::endl;

	return 0;
}


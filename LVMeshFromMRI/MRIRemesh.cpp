/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief From shortaxis contours in MRI + reference point creates a smoother mesh using splines longitudinally
    */
#include <vtkPolyData.h>

#include <vtkCellArray.h>

#include <vtkSmartPointer.h>

#include "VTKCommonTools.h"
#include "CommonTools.h"


#include <vtkPlane.h>

#include <vtkType.h>
#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkCellData.h>
#include <vtkCutter.h>
#include <vtkAppendPolyData.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

#include <vector>

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_cross.h>

//regionID
//%mark the uppermost contour (last ones)
//%epicardium: 5-apical contour, 6-intermediate, 7 - basal
//%endocardium: 0-apical contour, 1-intermediate, 2 - basal
//%reference point: 10

#define CONTOUR_ENDO_APICAL 0
#define CONTOUR_ENDO_BASAL 2
#define CONTOUR_ENDO 1

#define CONTOUR_EPI_APICAL 5
#define CONTOUR_EPI_BASAL 7
#define CONTOUR_EPI 6

#define REF_POINT 10


/* structure for qsort to sort an array and also have inices into the original array */ 
struct sorting_struct { 
    double value;
    int index;
};

/* qsort struct comparision function (product C-string field) */ 
int struct_cmp_by_value(const void *a, const void *b) 
{ 
    struct sorting_struct *ia = (struct sorting_struct *)a;
    struct sorting_struct *ib = (struct sorting_struct *)b;
	return (int)( ia->value - ib->value );    
} 
////////////////////////////////////////////////////////////////

//! compute centroid of polydata
void ComputeCentroid( vtkPolyData* shape, double* centroid );

/*! \brief Rotate a point p by angle theta around an arbitrary line segment p1-p2

   Return the rotated point.
   Positive angles are anticlockwise looking down the axis
   towards the origin.
   Assume right hand coordinate system. 
   */
vnl_vector<double> ArbitraryRotate(vnl_vector<double> p,double theta,vnl_vector<double> p1, vnl_vector<double> p2);

//! Something about point order for this specific application
vtkPolyData* OrderPoints( vtkPolyData* pts, int apical_id, int basal_id );

int main(int argc, char **argv)
{
	
	std::cout << "Version 0.99. From shortaxis contours in MRI + reference point creates a smoother mesh using splines longitudinally."<< std::endl;

	if( argc<2 ) 
	{
		std::cout << "Params: "<< std::endl;
		std::cout << "-shape <shape.vtk> \t\t- shape"<< std::endl;
		std::cout << "-i <int> \t\t- number of planes for intersections with contours"<< std::endl;
		std::cout << "-out <shape.vtk> \t\t- resulting shape"<< std::endl;
		
		return -1;
	}

	char *inshapefile=NULL;
	char *outshapefile=NULL;
	int nintersections = 8;
	
	for(int c=1; c<argc; c++)
	{
		if( strcmp(argv[c],"-shape")==0 )
		{
			inshapefile = argv[++c];
		}
		else if( strcmp(argv[c],"-i")==0 )
		{
			nintersections = atoi(argv[++c]);
		}
		else if( strcmp(argv[c],"-out")==0 )
		{
			outshapefile = argv[++c];
		}

	}


	std::cout<<"Shape: "<<inshapefile<<std::endl;
	std::cout<<"Output: "<<outshapefile<<std::endl;
	std::cout<<"Num. of intersections: "<<nintersections<<std::endl;
	
	vtkSmartPointer<vtkPolyData> inshape = vtkSmartPointer<vtkPolyData>::Take(
		CommonTools::LoadShapeFromFile(inshapefile));

	inshape->GetPointData()->SetActiveScalars("regionID");

	const double angle_increment = M_PI/nintersections;


	//////////////////////////////////////////////
	//
	// extract the centerline of the endocardium
	//
	vtkSmartPointer<vtkPolyData> endo_base = vtkSmartPointer<vtkPolyData>::Take( 
		CommonTools::GetShapeSubSurface(inshape, CONTOUR_ENDO_BASAL) );
	vtkSmartPointer<vtkPolyData> endo_apex = vtkSmartPointer<vtkPolyData>::Take( 
		CommonTools::GetShapeSubSurface(inshape, CONTOUR_ENDO_APICAL) );
	vtkSmartPointer<vtkPolyData> ref_point_pd = vtkSmartPointer<vtkPolyData>::Take( 
		CommonTools::GetShapeSubSurface(inshape, REF_POINT) );


	vnl_vector<double> endo_basal_cntr(3);
	vnl_vector<double> endo_api_cntr(3);
	vnl_vector<double> ref_point(3);

	//3 points that define initial plane
	ComputeCentroid( endo_base, endo_basal_cntr.begin() );
	ComputeCentroid( endo_apex, endo_api_cntr.begin() );
	ref_point_pd->GetPoint(0,ref_point.begin());

	vnl_vector<double> lax = endo_api_cntr - endo_basal_cntr;

	//this will be the cutting plane
	vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
	plane->SetOrigin(endo_basal_cntr.begin());

	//prepare the cutters
	vtkSmartPointer<vtkPolyData> endo = vtkSmartPointer<vtkPolyData>::Take( 
		CommonTools::GetShapeSubSurface(inshape, CONTOUR_ENDO_APICAL, CONTOUR_ENDO_BASAL ) );
	vtkSmartPointer<vtkPolyData> epi = vtkSmartPointer<vtkPolyData>::Take( 
		CommonTools::GetShapeSubSurface(inshape, CONTOUR_EPI_APICAL, CONTOUR_EPI_BASAL ) );

	vtkSmartPointer<vtkCutter> cutter_endo = vtkSmartPointer<vtkCutter>::New();
	vtkSmartPointer<vtkCutter> cutter_epi = vtkSmartPointer<vtkCutter>::New();
	cutter_endo->SetInputData(endo);
	cutter_epi->SetInputData(epi);

	//debug stuff - save the ntersection points for visualization
	vtkSmartPointer<vtkAppendPolyData> appender = vtkSmartPointer<vtkAppendPolyData>::New();

	for(int i=0; i<nintersections; i++)
	{
		double angle = i*angle_increment;
		
		//rotate the reference point about the long axis to get different cuts
		vnl_vector<double> rotated_pt = ArbitraryRotate( ref_point, angle, endo_basal_cntr, endo_api_cntr);
		//use endo_basal_cntr for origin and the other two points to compute the normal to define the cutting plane
		vnl_vector<double> sax = rotated_pt - endo_basal_cntr;
		vnl_vector<double> normal = vnl_cross_3d(lax, sax);
		normal.normalize();
		
		plane->SetNormal(normal.begin());
		cutter_endo->SetCutFunction( plane );
		cutter_epi->SetCutFunction( plane );
	
		cutter_endo->Update();
		cutter_epi->Update();
		
		vtkSmartPointer<vtkPolyData> endo_cut = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPolyData> epi_cut = vtkSmartPointer<vtkPolyData>::New();
		
		//order the points
		
		endo_cut->DeepCopy( cutter_endo->GetOutput() );
		epi_cut->DeepCopy( cutter_epi->GetOutput() );
		
		vtkSmartPointer<vtkPolyData> endo_cut_ordered = vtkSmartPointer<vtkPolyData>::Take(
			OrderPoints( endo_cut, CONTOUR_ENDO_APICAL, CONTOUR_ENDO_BASAL ) );
		
		vtkSmartPointer<vtkPolyData> epi_cut_ordered = vtkSmartPointer<vtkPolyData>::Take(
			OrderPoints( epi_cut, CONTOUR_EPI_APICAL, CONTOUR_EPI_BASAL ) );
			
		CommonTools::SavePolydata(endo_cut_ordered, outshapefile);
		
		appender->AddInputData( endo_cut_ordered );
		appender->AddInputData( epi_cut_ordered );
		
	}
	
	appender->Update();
	
	CommonTools::SavePolydata(appender->GetOutput(), outshapefile);

	return 0;

}


///////////////////////////////////////////////////////////////////////////////////////////
//
//   Functions
//


void ComputeCentroid( vtkPolyData* shape, double* centroid )
{
	centroid[0] = 0;
	centroid[1] = 0;
	centroid[2] = 0;
	
	const double k = 1.0/static_cast<double>(shape->GetNumberOfPoints());
	
	for(int i=0; i<shape->GetNumberOfPoints(); i++)
	{
		double pt[3];
		shape->GetPoint(i, pt);
		centroid[0] += pt[0]*k;
		centroid[1] += pt[1]*k;
		centroid[2] += pt[2]*k;
	}
}






/*
   Rotate a point p by angle theta around an arbitrary line segment p1-p2
   Return the rotated point.
   Positive angles are anticlockwise looking down the axis
   towards the origin.
   Assume right hand coordinate system.  
*/
vnl_vector<double> ArbitraryRotate(vnl_vector<double> p,double theta,vnl_vector<double> p1, vnl_vector<double> p2)
{
	vnl_vector<double> q(3,0);
	double costheta,sintheta;
	vnl_vector<double> r;

	r = p2-p1; 
	p = -p1;

	r.normalize();

	costheta = cos(theta);
	sintheta = sin(theta);

	q[0] += (costheta + (1 - costheta) * r[0] * r[0]) * p[0];
	q[0] += ((1 - costheta) * r[0] * r[1] - r[2] * sintheta) * p[1];
	q[0] += ((1 - costheta) * r[0] * r[2] + r[1] * sintheta) * p[2];

	q[1] += ((1 - costheta) * r[0] * r[1] + r[2] * sintheta) * p[0];
	q[1] += (costheta + (1 - costheta) * r[1] * r[1]) * p[1];
	q[1] += ((1 - costheta) * r[1] * r[2] - r[0] * sintheta) * p[2];

	q[2] += ((1 - costheta) * r[0] * r[2] - r[1] * sintheta) * p[0];
	q[2] += ((1 - costheta) * r[1] * r[2] + r[0] * sintheta) * p[1];
	q[2] += (costheta + (1 - costheta) * r[2] * r[2]) * p[2];

	q += p1;

   return(q);
}




vtkPolyData* OrderPoints( vtkPolyData* pts, int apical_id, int basal_id )
{
	vtkPolyData* result = vtkPolyData::New();
	vnl_matrix<double> points_in;
	
	CommonTools::ExportPolyDataPoints( pts, points_in );
	
	vtkShortArray* scalars = static_cast<vtkShortArray*>( pts->GetPointData()->GetArray("regionID") );
	
	//find the topmost points
	std::vector<int> basal;
	std::vector<int> apical; //there will be 2 points in the apical and basal contours
	std::vector<struct sorting_struct> pts4sorting;

	const int npoints = points_in.rows();

	for( int i=0; i<npoints; i++ )
	{
		if( scalars->GetValue(i)==basal_id ) basal.push_back(i);
		else if( scalars->GetValue(i)==apical_id ) apical.push_back(i);
		else
		{
			struct sorting_struct str;
			str.value = points_in[i][2];
			str.index = i;
			pts4sorting.push_back( str );
		}
	}
	
	vnl_matrix<double> points_out( points_in.rows(), points_in.cols() );
	
	//put first the points of the basal contour
	points_out.set_row( 0, points_in.get_row( basal[0] ) );
	points_out.set_row( npoints-1, points_in.get_row( basal[1] ) );

	//qsort other points based on Z coordinate
	qsort( &pts4sorting[0], pts4sorting.size(), sizeof(pts4sorting[0]), struct_cmp_by_value );
	
	for( int i=0; i<pts4sorting.size()/2; i++ )
	{
		points_out.set_row( i+1, points_in.get_row( pts4sorting[2*i].index ) );
		points_out.set_row( npoints-i-1, points_in.get_row( pts4sorting[2*i+1].index ) );
	}
	
	//Check the point pairs and flip when necessary
	for( int i=1; i<(npoints/2) -1; i++ )
	{
		vnl_vector<double> prev_point = points_out.get_row(i-1);
		vnl_vector<double> pt1 = points_out.get_row(i);
		vnl_vector<double> pt2 = points_out.get_row(npoints-i-1);
		double d1 = (prev_point-pt1).squared_magnitude();
		double d2 = (prev_point-pt2).squared_magnitude();
		if( d1>d2 )
		{
			points_out.set_row(i, pt2);
			points_out.set_row(npoints-i-1, pt1);
		}
	}

	result->DeepCopy( pts );
	CommonTools::ImportPolyDataPoints( result, points_out );
	
	//debug scalars
	vtkSmartPointer<vtkShortArray> ds = vtkSmartPointer<vtkShortArray>::New();
	ds->SetName("pointID");
	ds->SetNumberOfValues( npoints );
	for( int i=0; i<npoints; i++ ) ds->SetValue(i, i);
	result->GetPointData()->AddArray( ds );
	//:~ debug scalars
	
	//result->Update();
	
	return result;
}



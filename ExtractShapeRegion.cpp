/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Extract a part of polydata
*/
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkCleanPolyData.h>

#include <vtkPointData.h>
#include <vtkShortArray.h>
#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkObject.h>

#include "CommonTools.h"

//#include <meHoleFiller.h>

//#include <vtkstd/exception>

#include <vtkFeatureEdges.h>
#include <vtkCleanPolyData.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_svd_economy.h>

#include <vtkXMLPolyDataWriter.h>



#include <vtkFeatureEdges.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkAppendPolyData.h>

//------------------------------------------------------------------
int main( int argc, char *argv[] )
{	

	vtkObject::GlobalWarningDisplayOff();

	if( argc < 4 )
	{
		std::cout << "Subpart extractor v1.0 " << std::endl;
		std::cout << "Usage: " << std::endl;
		std::cout << "extract_region <input shape> <arrayname> <lower threshold> <upper threshold> <output shape> [closed]" << std::endl;
		std::cout << "Shapes must be stored in vtk or XML" << std::endl;
		std::cout << "output shape type is determined by its extension, make sure it's vtp(xml) or vtk(binary). ALL the letters of the extension should have the same case." << std::endl;
		return EXIT_FAILURE;
	}


	char type[4];
	int cc = 1;
	const char* filename = argv[cc++];
	const char* arrayname = argv[cc++];
	const char* low_thold = argv[cc++];
	const char* up_thold = argv[cc++];
	const char* out_shape = argv[cc++];

	const char* closed_flag = cc<argc?argv[cc++]:"";


	strcpy( type, &filename[strlen(filename)-3] );

	cout << " opening input shape" << filename << endl;

	vtkPolyData* shape = vtkPolyData::New();
	if( strcmp( type, "vtp" )==0 )
	{
		vtkXMLPolyDataReader *reader = vtkXMLPolyDataReader::New();
		reader -> SetFileName(filename);
		reader -> Update();
		shape->DeepCopy( reader->GetOutput() );
		reader->Delete();
	}
	else
	{
		vtkPolyDataReader *reader = vtkPolyDataReader::New();
		reader -> SetFileName(filename);
		reader -> Update();
		shape->DeepCopy( reader->GetOutput() );
		reader->Delete();
	}

	//cout << " opening region defs" << argv[2] << endl;

	//vtkShortArray* regions = vtkShortArray::New();
	//blVTKHelperTools::LoadVtkShortArray( argv[2], regions );
	//
	//if( shape->GetPointData()->GetArray( regions->GetName() ) == NULL )
	//{
	// shape->GetPointData()->AddArray( regions );
	//}

	//shape->GetPointData()->SetActiveScalars( regions->GetName() );
	shape->GetPointData()->SetActiveScalars( arrayname );

	float region_thold_lower = atof(low_thold);
	float region_thold_upper = atof(up_thold);
	cout << "extracting region [" << region_thold_lower<<", "<<region_thold_upper<<"]"<< endl;

	// vtkPolyData* subpart = blVTKHelperTools::GetShapeSubSurface(shape, region_thold_lower, region_thold_upper);
	vtkPolyData* subpart;

	if( strcmp(closed_flag,"closed") == 0 )
	{
		subpart = CommonTools::GetShapeSubSurface( shape, region_thold_lower, region_thold_upper);

		vtkPolyData *closed_subpart = CommonTools::CloseSurface( subpart );
		subpart->DeepCopy( closed_subpart );
		closed_subpart->Delete();

	}
	else
	{
		subpart = CommonTools::GetShapeSubSurface(shape, region_thold_lower, region_thold_upper);
	}



	char extension[4];
	strcpy( extension, &out_shape[strlen(out_shape)-3] );

	if( strcmp(extension,"vtp")==0 || strcmp(extension,"VTP")==0 )
	{
		vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
		writer->SetInputData( subpart );
		writer->SetFileName( out_shape );
		writer->SetDataModeToAscii();
		writer->Write();

		writer->Delete();
	}
	else if( strcmp(extension,"vtk")==0 || strcmp(extension,"VTK")==0 )
	{
		vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
		writer->SetInputData( subpart );
		writer->SetFileName( out_shape );
		writer->SetFileTypeToBinary();
		writer->Write();

		writer->Delete();
	}
	else
	{
		std::cout<<"Unable to determine output file format. Check the extension."<<std::endl;
	}

	shape->Delete();
	subpart->Delete();


	return EXIT_SUCCESS;
}


//vtkPolyData* FindEdge( vtkPolyData* shape )
//{
//	vtkIdList* shape_edge_ids = vtkIdList::New();
//
//	shape->BuildLinks();
//
//	vtkIdList* 	neighbour_cells = vtkIdList::New();	 
//
//	//iterate over all the cells and edges
//	for( int cell_no=0; cell_no<shape->GetNumberOfCells(); cell_no++)
//	{
//		vtkCell* cell = shape->GetCell( cell_no );
//
//		const int no_edges = cell->GetNumberOfEdges();
//		for( int edge_no=0; edge_no<no_edges; edge_no++ )
//		{
//			vtkCell* edge = cell->GetEdge( edge_no );
//			vtkIdList* pt_ids = edge->GetPointIds();
//
//			shape->GetCellEdgeNeighbors( cell_no, pt_ids->GetId(0), pt_ids->GetId(1), neighbour_cells );
//
//			if( neighbour_cells->GetNumberOfIds()==0 ) //a shape edge has been found
//			{
//				for( int i=0; i<2; i++ )
//					shape_edge_ids->InsertUniqueId( pt_ids->GetId(i) );
//			}
//
//		}
//
//	}
//	neighbour_cells->Delete();
//
//	vtkPoints* pts = vtkPoints::New();
//	for( int i=0; i<shape_edge_ids->GetNumberOfIds(); i++ )
//		pts->InsertNextPoint( shape->GetPoint( shape_edge_ids->GetId(i) ) );
//
//
//	shape_edge_ids->Delete();
//
//	return blVTKHelperTools::Points2Polydata( pts );
//	
//}



/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file CommonTools.h
    \brief Some common functions
*/
#ifndef __CommonTools_h
#define __CommonTools_h


	//---------------------------------------------------
	// HEADERS
	//---------------------------------------------------

	#include <vector>
	#include <vnl/vnl_vector.h>
	
	#include <vtkImageData.h>


	//---------------------------------------------------
	// FORWARD DECLARATIONS
	//---------------------------------------------------

	class vtkShortArray;
	class vtkPolyData;
	class vtkUnstructuredGrid;
	class vtkPoints;
	class vtkDataSet;
	class vtkStringArray;


	//---------------------------------------------------
	// CLASS DEFINITION
	//---------------------------------------------------


	namespace CommonTools
	{

    	//! Valid volume formats io functions can handle 
		enum VTKSurfaceMeshFormats { UnknownType, VTKPolyDataType, VTKXMLPolyDataType, STLType, PLYType };
		
    	//! Valid volume formats io functions can handle 
		enum VTKVolumeMeshFormats { UnknownVolumeType, VTKUnstructuredGridType };
		
		/**
		\brief Save a vtk short array using a specific format 
		\note If there's an error, an std::exception is thrown
		*/
		void SaveVtkShortArray( const char *filename, vtkShortArray* the_array );

		/**
		\brief Load a vtk short array using a specific format 
		\note If there's an error, an std::exception is thrown
		 */
		void LoadVtkShortArray( const char *filename, vtkShortArray* the_array );


		//! Call to GetShapeSubSurface( ) with nSubPart-0.1, nSubPart+0.1
		vtkPolyData* GetShapeSubSurface(vtkPolyData * inputShape, unsigned int nSubPart);

		/**
		\brief Apply threshold to extract the subpart, apply 
		vtkDataSetSurfaceFilter and vtkCleanPolyData
		\note The caller to this function should delete the output shape
		*/
		vtkPolyData* GetShapeSubSurface(vtkPolyData * inputShape, double tholdLower, double tholdUpper);

		//! Closes only 1 hole, make sure there are no more. 
		vtkPolyData* CloseSurface(vtkPolyData* shape);

		//! Generate a cover for a small hole. Uses centroid. 
		vtkPolyData* GenerateHoleCover(vtkPolyData* edge);


		//! Save polydata to a file. 
		void SavePolydata( vtkPolyData* poly, const char* filename, bool binary = false );
		
		//! Save image to a file. 
		void SaveImage( vtkDataSet* image, const char* filename );

		//! Load image from a file. 
		vtkImageData* LoadImage( const char* filename );

		//! Save unstructured grid to a file. 
		void SaveUnstructuredGrid( vtkUnstructuredGrid* grid, const char* filename );


		//! Save vtkPoints to a file. 
		void SavePoints( vtkPoints* pts, const char* filename );

		//! Check if file exists and throw an exception if needed 
		bool FileExists( const char* filename, bool no_exception = false );

		//! Generate a filelist
		void ReadFilelist( const char* file, std::vector<std::string>& list, bool check_existence = false );


		/**
		\brief Saves points such that they can e visualized in paraview
		also saves a scalar corresponding to a position so it is easy to see the ordering of points 
		(if sampling is correct)
		if scalars != NULL, the corresponding scalar values will be assigned to the points
		\note The caller to this function should call points->Delete( )
		*/
		vtkPolyData* Points2Polydata( vtkPoints* points, double scalar );

		/**
		\brief Saves points such that they can e visualized in paraview
		also saves a scalar corresponding to a position so it is easy to see the ordering of points 

		(if sampling is correct)
		if scalars != NULL, the corresponding scalar values will be assigned to the points
		\note The caller to this function should call points->Delete( )
		*/
		vtkPolyData* Points2Polydata( vtkPoints* points, const double* scalars=NULL );

		//! extract points from polydata
		void ExportPolyDataPoints( vtkPolyData* shape, vnl_vector<double>& points );
		//! extract points from polydata
		void ExportPolyDataPoints( vtkPolyData* shape, vnl_matrix<double>& points );
		//!copy points to polydata
		void ImportPolyDataPoints( vtkPolyData* shape, vnl_vector<double>& points );
		//!copy points to polydata
		void ImportPolyDataPoints( vtkPolyData* shape, vnl_matrix<double>& points );


		//! generate a Meshlab script for mesh decimation, specify number of faces
		void GenerateDecimationScript( const char* filename, int nfaces );



		//! Rescale polydata
		void ScaleShape(vtkPolyData* shapein, vtkPolyData* shapeout, float scale, bool centerAfterScale = false); 
		//! Resize image
		void ShrinkImage(vtkDataSet* imagein, vtkDataSet* imageout, int factor);
		//! Resize unstructured grid
		void ScaleVolume(vtkUnstructuredGrid* volumein, vtkUnstructuredGrid* volumeout, float scale); 


		//! Load polydata from file
		vtkPolyData* LoadShapeFromFile(const char *shapeFileName);
		//! Load unstructured grid from file
		vtkUnstructuredGrid* LoadVolumeFromFile(const char *volumeFileName);
		//! Save unstructured grid to file
		void SaveVolumeToFile( vtkUnstructuredGrid *volumePt, const char *volumeFileName, const char *header);
		//! Save polydata to file
		void SaveShapeToFile(
								vtkPolyData *shapePt, 
								const char *shapeFileName,
								const char *header=NULL );
								
		//! Calculate point-to-surface distance
		void GetP2S(vtkPolyData * manualPt, vtkPolyData *segmentedPt, 
								double& mean, double& std_dev, double& max, 
								double& last, bool b_array);

		//! Calculate point-to-point distance
		void GetP2P(vtkPolyData * manualPt, vtkPolyData *segmentedPt, 
								double& mean, double& std_dev, double& max, 
								double& last, bool b_array);
		//! Calculate point-to-surface distance
		vtkPolyData* GetP2S(vtkPolyData *shapePt1, vtkPolyData *shapePt2, std::vector< vnl_vector<double> >& distances);
		//! Calculate point-to-surface distance
		vtkPolyData* GetP2S(vtkPolyData *shapePt1, vtkPolyData *shapePt2, vnl_vector<double>& distances);
		//! Calculate surface-to-surface distance
		void GetS2S(vtkPolyData *shapePt1, vtkPolyData *shapePt2, std::vector< vnl_vector<double> >& distances );

        //! check if the file has valid extension for saving. To remove one day.
		bool CheckSaveFileExtension(const char *shapeFileName);
        //! identify volume data type
		VTKVolumeMeshFormats GetTypeOfVTKVolumeData(const char *volumeFileName);
        //! identify VTK data type
		VTKSurfaceMeshFormats GetTypeOfVTKData(const char *shapeFileName);
		
		
		/*! \brief Explicit solution to Laplace Eq. (c) Ruben Cardenes + Constantine Butakoff
		
		 Explicit solution to Laplace Eq. (c) Ruben Cardenes + Constantine Butakoff
		3 Outside domain
		1 Exterior boundary 
		0 Interior boundary
		2 Inside domain 
		*/ 
		//-------------------------------------------------------------------------------
		int laplace3D_voxelsize( vtkImageData* inputImage, vtkImageData* outputImage, int iterations/*=100*/ );
		//-------------------------------------------------------------------------------

	}

#endif

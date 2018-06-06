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
#ifndef __VTKCommonTools_h
#define __VTKCommonTools_h


//---------------------------------------------------
// HEADERS
//---------------------------------------------------

#include <vector>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkAlgorithm.h>


//#define LINEBREAK "\x0D\x0A"
#define LINEBREAK "\x0A"

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


namespace CommonTools {

    //! Valid volume formats io functions can handle 

    enum VTKSurfaceMeshFormats {
        UnknownType, VTKPolyDataType, VTKXMLPolyDataType, STLType, PLYType
    };

    //! Valid volume formats io functions can handle 

    enum VTKVolumeMeshFormats {
        UnknownVolumeType, VTKUnstructuredGridType
    };

    void CommandLineProgressIndicator(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData); 
    vtkSmartPointer<vtkCallbackCommand> AssociateProgressFunction(vtkAlgorithm* algorithm, vtkCallbackCommand* existing_function = NULL);
    
    
    //! Load polydata from file
    vtkPolyData* LoadShapeFromFile(const char *shapeFileName);
    //! Load unstructured grid from file
    vtkUnstructuredGrid* LoadVolumeFromFile(const char *volumeFileName);
    //! Save unstructured grid to file
    void SaveVolumeToFile(vtkUnstructuredGrid *volumePt, const char *volumeFileName, const char *header);
    //! Save polydata to file
    void SaveShapeToFile(
            vtkPolyData *shapePt,
            const char *shapeFileName,
            const char *header = NULL);

    //! check if the file has valid extension for saving. To remove one day.
    bool CheckSaveFileExtension(const char *shapeFileName);
    //! identify volume data type
    VTKVolumeMeshFormats GetTypeOfVTKVolumeData(const char *volumeFileName);
    //! identify VTK data type
    VTKSurfaceMeshFormats GetTypeOfVTKData(const char *shapeFileName);

    void SaveVolMeshBSC(const char* infile, const char* outfile_prefix, float scale, bool correct_orientation=true);
    void SaveVolMeshBSC(vtkDataSet* volmesh, const char* outfile_prefix, float scale, bool correct_orientation=true);

}

#endif

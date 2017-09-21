/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#ifndef VTKCOMMANDLINEPROGRESS_H
#define VTKCOMMANDLINEPROGRESS_H

#endif /* VTKCOMMANDLINEPROGRESS_H */






//
// GenerateCubesFromLabels
//   Usage: GenerateCubesFromLabels InputVolume Startlabel Endlabel
//          where
//          InputVolume is a meta file containing a 3 volume of
//            discrete labels.
//          StartLabel is the first label to be processed
//          EndLabel is the last label to be processed
//          NOTE: There can be gaps in the labeling. If a label does
//          not exist in the volume, it will be skipped.
//
//
#include <vtkMetaImageReader.h>
#include <vtkImageAccumulate.h>
#include <vtkImageWrapPad.h>
#include <vtkMaskFields.h>
#include <vtkThreshold.h>
#include <vtkTransformFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>

#include <vtkTransform.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <sstream>

#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkImageData.h>
#include <vtkCallbackCommand.h>


#include "vtkPolyDataAlgorithm.h"

/* To use, call  AssociateProgressFunction for the objects that you need
 */






class vtkTestProgressReportFilter : public vtkPolyDataAlgorithm {
public:
    static vtkTestProgressReportFilter *New();
    vtkTypeMacro(vtkTestProgressReportFilter, vtkAlgorithm);

protected:

    vtkTestProgressReportFilter() {
    }

    ~vtkTestProgressReportFilter() VTK_OVERRIDE {
    }

    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *) VTK_OVERRIDE;

private:
    vtkTestProgressReportFilter(const vtkTestProgressReportFilter&) VTK_DELETE_FUNCTION;
    void operator=(const vtkTestProgressReportFilter&) VTK_DELETE_FUNCTION;

};




/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file CloseBVMesh.cpp  
    \brief Closes biventricular mesh by connecting endocardial edge to epicardial at the base. Does not close ventricles.

  It was made for Rafa's biventricular model. The mesh must have epicardium, rv endo and lv endo separable. No scalars are necessary. 
  The scalars must be vtkShortArray!!! (type short)
 */
#include "CommonTools.h"

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkShortArray.h>
#include <vtkCellData.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolygon.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDataWriter.h>



//------------------------------------------------------------------

int main(int argc, char *argv[]) {
    const char* mesh_file = argv[1];
    const char* output_file = argv[2];

    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(mesh_file);
    rdr->Update();
    
    vtkPolyData* pd = rdr->GetOutput();
   

    //feature edges, all off , leave boudary only
    vtkSmartPointer<vtkFeatureEdges> fe = vtkSmartPointer<vtkFeatureEdges>::New();
    fe->SetInputData(pd);
    fe->BoundaryEdgesOn();
    fe->FeatureEdgesOff();
    fe->NonManifoldEdgesOff();
    fe->ManifoldEdgesOff();
    fe->Update();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> split = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    split->SetInputData(fe->GetOutput());
    split->SetExtractionModeToAllRegions();
    split->ColorRegionsOn();
    split->Update();

    vtkSmartPointer<vtkAppendPolyData> final = vtkSmartPointer<vtkAppendPolyData>::New();
    //final->AddInputData (pd);

    for (int i = 0; i < split->GetNumberOfExtractedRegions(); i++) {
        vtkSmartPointer<vtkPolyData> sub = vtkSmartPointer<vtkPolyData>::Take(
                CommonTools::GetShapeSubSurface(split->GetOutput(), i, i));

        //		char filename[100];
        //		sprintf(filename,"sub %03d.vtk",i);
        //		CommonTools::SaveShapeToFile( sub, filename, NULL);

        std::cout << "Hole id :" << i << " pts: " << sub->GetNumberOfPoints() << std::endl;

        /*
                        vtkSmartPointer<vtkStripper> strip = vtkSmartPointer<vtkStripper>::New();
                        strip->SetInputData(sub);
                        strip->Update();
		
                        vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
                        cleaner->SetInputData(strip->GetOutput());
                        cleaner->Update();
         */
        vtkSmartPointer<vtkPolyData> cover = vtkSmartPointer<vtkPolyData>::Take(
                CommonTools::GenerateHoleCover(sub));


        final->AddInputData(cover);


        //		sprintf(filename,"delaunay_%03d.vtk",i);
        //		CommonTools::SaveShapeToFile( cover, filename, NULL);

    }

    final->Update();

    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(final->GetOutput());
    cleaner->Update();

    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(output_file);
    wr->SetInputData(cleaner->GetOutput());
    wr->Write();

    return 0;
}


/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/



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
#include <vtkImageAccumulate.h>
#include <vtkMaskFields.h>
#include <vtkThreshold.h>
#include <vtkSmartPointer.h>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <sstream>

#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkImageData.h>

#include <VTKCommonTools.h>
#include <vtkCallbackCommand.h>

#include <vtkEnSightWriter.h>
#include <vtkShortArray.h>

void SaveUnstructuredGrid(vtkUnstructuredGrid *data, const char* filename) {

    bool FoundCellData = false;

    if (data->GetCellData()) {
        if (data->GetCellData()->GetArray("BlockId")) {
            FoundCellData = true;
        }
    }
    if (FoundCellData == false) {
        //Looks like we have to create it.. crap
        vtkSmartPointer<vtkShortArray> bids = vtkSmartPointer<vtkShortArray>::New();
        bids->SetName("BlockId");
        for (vtkIdType i = 0; i < data->GetNumberOfCells(); i++)
            bids->InsertNextValue(1);
        data->GetCellData()->SetScalars(bids);
    }

    int blockids[2];
    blockids[0] = 1;
    blockids[1] = 0;
    
    vtkSmartPointer<vtkEnSightWriter> ptr_ens = vtkSmartPointer<vtkEnSightWriter>::New();
    CommonTools::AssociateProgressFunction(ptr_ens);
    ptr_ens->SetNumberOfBlocks(1);
    ptr_ens->SetBlockIDs(blockids);
    ptr_ens->SetTimeStep(0);
    ptr_ens->SetInputData(data);
    ptr_ens->SetFileName(filename);

    ptr_ens->Write();
    ptr_ens->WriteCaseFile(1);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        cout << "Usage: " << argv[0] << " InputVolume StartLabel EndLabel outprefix tetra/hexa" << endl;
        return EXIT_FAILURE;
    }



    // Create all of the classes we will need
    vtkSmartPointer<vtkDataSetReader> reader =
            vtkSmartPointer<vtkDataSetReader>::New();
    vtkSmartPointer<vtkImageAccumulate> histogram =
            vtkSmartPointer<vtkImageAccumulate>::New();
    vtkSmartPointer<vtkMaskFields> scalarsOff =
            vtkSmartPointer<vtkMaskFields>::New();
    vtkSmartPointer<vtkThreshold> selector =
            vtkSmartPointer<vtkThreshold>::New();
    vtkSmartPointer<vtkDataSetTriangleFilter> triangulator =
            vtkSmartPointer<vtkDataSetTriangleFilter>::New();



    vtkSmartPointer<vtkCallbackCommand> progressCallback =
            CommonTools::AssociateProgressFunction(reader);
    CommonTools::AssociateProgressFunction(selector, progressCallback);
    CommonTools::AssociateProgressFunction(triangulator, progressCallback);
    CommonTools::AssociateProgressFunction(histogram, progressCallback);



    // Define all of the variables
    unsigned int startLabel = atoi(argv[2]);
    if (startLabel > VTK_SHORT_MAX) {
        std::cout << "ERROR: startLabel is larger than " << VTK_SHORT_MAX << std::endl;
        return EXIT_FAILURE;
    }
    unsigned int endLabel = atoi(argv[3]);
    if (endLabel > VTK_SHORT_MAX) {
        std::cout << "ERROR: endLabel is larger than " << VTK_SHORT_MAX << std::endl;
        return EXIT_FAILURE;
    }
    std::string filePrefix = "Cubes";



    bool make_tetras = true;
    if (strcmp(argv[5], "hexa") == 0) make_tetras = false;


    // Generate cubes from labels
    // 1) Read the meta file
    // 2) Generate a histogram of the labels
    // 3) Convert point data to cell data
    // 4) Output each cube model into a separate file

    reader->SetFileName(argv[1]);
    reader->Update();

    histogram->SetInputConnection(reader->GetOutputPort());
    histogram->SetComponentExtent(0, endLabel, 0, 0, 0, 0);
    histogram->SetComponentOrigin(0, 0, 0);
    histogram->SetComponentSpacing(1, 1, 1);
    histogram->Update();


    // Copy the scalar point data of the volume into the scalar cell data

    selector->SetInputData(reader->GetOutput());
    selector->SetInputArrayToProcess(0, 0, 0,
            vtkDataObject::FIELD_ASSOCIATION_POINTS,
            vtkDataSetAttributes::SCALARS);






    for (unsigned int i = startLabel; i <= endLabel; i++) {
        // see if the label exists, if not skip it
        double frequency =
                histogram->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i);
        if (frequency == 0.0) {
            continue;
        }

        // select the cells for a given label
        selector->ThresholdBetween(i, i);
        selector->Update();

        // Strip the scalars from the output
        scalarsOff->SetInputData(selector->GetOutput());
        scalarsOff->CopyAttributeOff(vtkMaskFields::POINT_DATA,
                vtkDataSetAttributes::SCALARS);
        scalarsOff->CopyAttributeOff(vtkMaskFields::CELL_DATA,
                vtkDataSetAttributes::SCALARS);
        scalarsOff->Update();




        // output the polydata
        std::stringstream ss;
        ss << filePrefix << i; // << ".vtk";
        cout << argv[0] << " writing " << ss.str() << endl;


        if (make_tetras) {
            triangulator->SetInputConnection(scalarsOff->GetOutputPort());
            triangulator->Update();

            SaveUnstructuredGrid( triangulator->GetOutput(), ss.str().c_str() ) ;
        } else {
            SaveUnstructuredGrid( vtkUnstructuredGrid::SafeDownCast(scalarsOff->GetOutput()), ss.str().c_str() ) ;
        }

        std::cout << "Writing finished" << std::endl;
    }
    return EXIT_SUCCESS;
}

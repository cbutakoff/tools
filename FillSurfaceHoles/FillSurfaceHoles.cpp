/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
//Remesh-like hole filling
//based on Eurographics Symposium on Geometry Processing(2003), Filling Holes in Meshes, Peter Liepa
//and Filling Gaps in the Boundary of a Polyhedron, Gill Barequet, Micha Sharir
#include "SurfaceHoleFiller.h"

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>



//======================================================================
//
//
//      Main function
//
//
//
//======================================================

int main(int argc, char **argv) {
    std::cout << "Hole filling based on Eurographics Symposium on Geometry Processing(2003), Filling Holes in Meshes, Peter Liepa" << std::endl;
    std::cout << "usage: FillSurfaceHoles -i mesh.vtk -o mesh.vtk -smooth [none|cot|edgelen]" << std::endl;
    std::cout << "Default: smoothing on, cotangent weights" << std::endl;
    std::cout << "Smoothing: cotangent - standard cotangent weights, angle preserving, edgelen - inverse edge length, god knows why it's good, but authors say that in some cases it's nicer" << std::endl;
    std::cout << "version 1.0" << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;

    if (argc < 3) {
        return -1;
    }

    char *in_filename = NULL;
    char *out_filename = NULL;

    bool do_smoothing = true;
    int smoothing_type = EDGE_WEIGHT_COTANGENT;
    
    
    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-i") == 0) {
            in_filename = argv[++c];
        } else if (strcmp(argv[c], "-o") == 0) {
            out_filename = argv[++c];
        } else if (strcmp(argv[c], "-smooth") == 0) {
            const char* methodtype = argv[++c];
            if (strcmp(methodtype, "none") == 0)
                do_smoothing = false;
            if (strcmp(methodtype, "cot") == 0)
                smoothing_type = EDGE_WEIGHT_COTANGENT;
            if (strcmp(methodtype, "edgelen") == 0)
                smoothing_type = EDGE_WEIGHT_InvEdgeLength;
            else
                std::cout << "Unknown filling algorithm, using default." << std::endl;
        }
    }

    std::cout << "Input file: " << in_filename << std::endl;
    std::cout << "Output file: " << out_filename << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;

    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(in_filename);

    if (!rdr->IsFilePolyData()) {
        std::cout << "Input file is not plydata. Verify your file. Aborting" << std::endl;
        return -1;
    }

    rdr->Update();

    SurfaceHoleFiller hole_filler;
    hole_filler.SetInput(rdr->GetOutput());
    hole_filler.SetSmoothing(do_smoothing);
    hole_filler.SetEdgeWeightingType(smoothing_type);
    hole_filler.Update();
    

    //save the result
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(out_filename);
    wr->SetInputData(hole_filler.GetOutput());
    wr->Write();

    return 0;
}




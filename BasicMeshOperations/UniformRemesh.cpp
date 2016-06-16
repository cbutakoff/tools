/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Generate layers of Left atrium. 
 */
#include <vtkCleanPolyData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include "CommonTools.h"
#include "vtkSmartPointer.h"
#include <vtkvmtkPolyDataSurfaceRemeshing.h>


//! Uniform remeshing. Requires admesh in path for mesh fixing
void UniformRemesh(vtkPolyData* mesh, float targetarea = 1.0);

void usage(char *exe) {
    std::cout << "Generate layers of the atrium. " << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "-i <mesh.vtk> \t - input mesh" << std::endl;
    std::cout << "-o <mesh.vtk> \t - output mesh" << std::endl;
    std::cout << "-ta <float> \t - target cell area" << std::endl;

    exit(0);
}

int main(int argc, char **argv) {
    std::cout << "Version 1.0" << std::endl;
    if (argc < 4) usage(argv[0]);

    char* ifilename;
    char* ofilename;
    float targetarea = 1.0;


    std::cout << "Parsing parameters" << std::endl;
    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-i") == 0) {
            ifilename = argv[++c];
        } else if (strcmp(argv[c], "-o") == 0) {
            ofilename = argv[++c];
        } else if (strcmp(argv[c], "-ta") == 0) {
            targetarea = atof(argv[++c]);
        }
    }


    std::cout << "Parameters parsed" << std::endl;

    std::cout << "Checking for file existence" << std::endl;
    CommonTools::FileExists(ifilename);

    std::cout << "Loading " << ifilename << std::endl;
    vtkSmartPointer<vtkPolyData> imesh = vtkSmartPointer<vtkPolyData>::Take(CommonTools::LoadShapeFromFile(ifilename));



    std::cout << "Cleaning the mesh..." << std::endl;

    vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputData(imesh);
    clean->Update();

    imesh->DeepCopy(clean->GetOutput());

    //generate normals
    if (targetarea > 0) {
        UniformRemesh(imesh, targetarea);
    }

    clean->SetInputData(imesh);
    clean->Update();

   
    CommonTools::SaveShapeToFile(clean->GetOutput(),ofilename);

}

void UniformRemesh(vtkPolyData* mesh, float targetarea) {
    char filename[100];
    //post process: fix and uniformly remesh
    sprintf(filename, "tempshape3274899816.stl");
    CommonTools::SaveShapeToFile(mesh, filename);

    char fixcmd[100];
    std::cout << "Running admesh.." << std::endl;
    sprintf(fixcmd, "admesh -n -u -f -d %s -b%s -t0.1 -i50", filename, filename);
    int rubbish = system(fixcmd);

    vtkSmartPointer<vtkPolyData> temp_mesh = vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(filename));

    remove(filename);

    //remesh
    vtkSmartPointer<vtkvmtkPolyDataSurfaceRemeshing> remesh =
            vtkSmartPointer<vtkvmtkPolyDataSurfaceRemeshing>::New();
    remesh->SetInputData(temp_mesh);
    remesh->SetNumberOfIterations(10);
    remesh->SetElementSizeModeToTargetArea();
    remesh->SetTargetArea(targetarea);
    remesh->Update();

    mesh->DeepCopy(remesh->GetOutput());
}

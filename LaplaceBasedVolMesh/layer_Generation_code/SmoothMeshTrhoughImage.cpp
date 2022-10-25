/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Smooth mesh using image mask
 */
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

//#define USE_VMTK

#include "vtkPolyData.h"
#include "vtkCleanPolyData.h"
#include "VTKCommonTools.h"
#include "CommonTools.h"
// here go your #includes
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkImageMarchingCubes.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkDataSetWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataNormals.h>

//#include "itkPluginUtilities.h"


#include <iostream>
#include <fstream>
#include <stdio.h>
#include "vtkDataSetSurfaceFilter.h"

#ifdef USE_VMTK
#include <vtkvmtkPolyDataSurfaceRemeshing.h>
#endif


// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace {
} // end of anonymous namespace

int main(int argc, char * argv[]) {
    std::cout << "Version 1.0" << std::endl;

    if (argc < 2) {
        std::cout << "Params: " << std::endl;
        std::cout << "-i <shape.vtk> \t\t- shape" << std::endl;
        std::cout << "-o <shape.vtk> \t\t- resulting shape" << std::endl;
        std::cout << "-v <float> \t\t- voxel size (def=2.5)" << std::endl;
        std::cout << "-f <int> \t\t- Target number of faces (def=4000)" << std::endl;
        std::cout << "-s <float> \t\t- smoothing sigma (def=3)" << std::endl;
        std::cout << "-m <float> \t\t- threshold for marching cubes in [0,1] (def=0.5)" << std::endl;
        std::cout << "-ta <float> \t\t - target area for vmtk remesh (only for layers, default 1, <0 - disable, requires admesh)" << std::endl;

        return -1;
    }

    char *inputGeometry = NULL;
    char *outputGeometry = NULL;
    float VoxelSize = 0.5;
    float GaussSigma = 3.0;
    int nFaces = 4000;
    float mcthold = 0.5;
    float targetarea = 1.0;


    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-i") == 0) {
            inputGeometry = argv[++c];
        } else if (strcmp(argv[c], "-o") == 0) {
            outputGeometry = argv[++c];
        } else if (strcmp(argv[c], "-v") == 0) {
            VoxelSize = atof(argv[++c]);
        } else if (strcmp(argv[c], "-s") == 0) {
            GaussSigma = atof(argv[++c]);
        } else if (strcmp(argv[c], "-f") == 0) {
            nFaces = atoi(argv[++c]);
        } else if (strcmp(argv[c], "-s") == 0) {
            GaussSigma = atof(argv[++c]);
        } else if (strcmp(argv[c], "-m") == 0) {
            mcthold = atof(argv[++c]);
        } else if (strcmp(argv[c], "-ta") == 0) {
            targetarea = atof(argv[++c]);
        }
    }


    //VoxelSize
    //GaussSigma

    vtkSmartPointer< vtkPolyData > inshape = vtkSmartPointer< vtkPolyData > ::Take(CommonTools::LoadShapeFromFile(inputGeometry));

    std::cout << "Cleaning the mesh" << std::endl;

    vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputData(inshape);
    clean->Update();

    vtkSmartPointer< vtkPolyData > polyOut = vtkSmartPointer< vtkPolyData > ::New();

    polyOut->SetPoints(clean->GetOutput()->GetPoints());
    polyOut->SetPolys(clean->GetOutput()->GetPolys());

    /*
            std::cout<<"Filling the holes"<<std::endl;

            meCloseHoles::Pointer meshcloser = meCloseHoles::New();
            meshcloser->SetAlgorithm(meCloseHoles::LINEAR_TO_CENTER);
            meshcloser->SetInputData(polyOut);
            meshcloser->Update();
     */

    double origin[3];
    double spacing[3];
    int wextent[6];

    const int padding = 10; //space to add to shape extremes

    int fg = 1;
    int bg = 0;

    //std::cout<<"Foreground: "<<fg<<std::endl;
    //std::cout<<"Background: "<<bg<<std::endl;
    std::cout << "Voxel: " << VoxelSize << std::endl << std::endl;


    vtkSmartPointer<vtkImageData> inimage = vtkSmartPointer<vtkImageData>::New();

    std::cout << "Generating an image" << std::endl;

    spacing[0] = VoxelSize;
    spacing[1] = VoxelSize;
    spacing[2] = VoxelSize;

    double bounds[6];
    inshape->GetBounds(bounds);
    double t[3];
    t[0] = -bounds[0] + padding * spacing[0];
    t[1] = -bounds[2] + padding * spacing[1];
    t[2] = -bounds[4] + padding * spacing[2];


    origin[0] = -t[0];
    origin[1] = -t[1];
    origin[2] = -t[2];
    int dims[3];
    inshape->GetBounds(bounds);
    dims[0] = floor((bounds[1] - bounds[0]) / spacing[0]) + padding * 2;
    dims[1] = floor((bounds[3] - bounds[2]) / spacing[1]) + padding * 2;
    dims[2] = floor((bounds[5] - bounds[4]) / spacing[2]) + padding * 2;

//    inimage->SetScalarTypeToUnsignedChar();
    inimage->SetOrigin(origin);
    inimage->SetSpacing(spacing);
    inimage->SetDimensions(dims);
//    inimage->SetScalarTypeToFloat();
//    inimage->Update();
//    inimage->AllocateScalars();
    inimage->AllocateScalars(VTK_FLOAT, 1);
    inimage->GetExtent(wextent);

    float* voxels = static_cast<float*> (inimage->GetScalarPointer());
    for (int i = 0; i < dims[0] * dims[1] * dims[2]; i++) {
        voxels[i] = fg;
    }


    vtkSmartPointer<vtkPolyDataToImageStencil> filt =
            vtkSmartPointer<vtkPolyDataToImageStencil>::New();


    filt->SetOutputOrigin(origin);
    filt->SetOutputSpacing(spacing);
    filt->SetOutputWholeExtent(wextent);
    //	filt->SetInputData( meshcloser->GetOutput() );
    filt->SetInputData(polyOut);
    filt->Update();

    vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
    stencil->SetStencilData(filt->GetOutput());
    stencil->SetBackgroundValue(bg);
    stencil->SetInputData(inimage);

    std::cout << "Smoothing the image" << std::endl;

    vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    smooth->SetStandardDeviation(GaussSigma);
    smooth->SetInputData(stencil->GetOutput());
    smooth->Update();

    //vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
    //writer->SetInputData( smooth->GetOutput() );
    //writer->SetFileTypeToBinary();
    //writer->SetFileName( "image.vtk");
    //writer->Update();

    std::cout << "Running marching cubes" << std::endl;

    vtkSmartPointer<vtkImageMarchingCubes> marchingcubes = vtkSmartPointer<vtkImageMarchingCubes>::New();
    marchingcubes->SetInputData(smooth->GetOutput());
    marchingcubes->ComputeNormalsOn();
    marchingcubes->SetValue(0, mcthold);
    marchingcubes->Update();



    //reduce the number of faces
    vtkSmartPointer<vtkPolyData> tempPD;
    if (nFaces > 0) {
        ////cleanup
        //std::cout<<"Cleaning the mesh"<<std::endl;

        //clean->SetInputData(marchingcubes->GetOutput());
        //clean->Update();
        //polyOut->SetPoints(clean->GetOutput()->GetPoints());
        //polyOut->SetPolys(clean->GetOutput()->GetPolys());

        char temp_filename[200];
        strcpy(temp_filename, outputGeometry);
        strcat(temp_filename, ".ply");


        CommonTools::SaveShapeToFile(marchingcubes->GetOutput(), temp_filename);
        CommonTools::GenerateDecimationScript("decimation_script_196239462.mlx", nFaces);
        char cmdline[1000];
        sprintf(cmdline, "meshlabserver -i %s -o %s -s decimation_script_196239462.mlx", temp_filename, temp_filename);
        int rubbish = system(cmdline);
        tempPD = vtkSmartPointer< vtkPolyData > ::Take(CommonTools::LoadShapeFromFile(temp_filename));
        remove("decimation_script_196239462.mlx");
        remove(temp_filename);
    } else
        tempPD = marchingcubes->GetOutput();



    //cleanup
    std::cout << "Cleaning the mesh" << std::endl;

    clean->SetInputData(tempPD);
    clean->Update();
    polyOut->SetPoints(clean->GetOutput()->GetPoints());
    polyOut->SetPolys(clean->GetOutput()->GetPolys());


    //generate normals
    vtkSmartPointer<vtkPolyDataNormals> normalgen = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalgen->SetInputData(polyOut);
    normalgen->SplittingOff();
    normalgen->Update();

    vtkSmartPointer<vtkDataSetSurfaceFilter> extrsurf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    extrsurf->SetInputData(normalgen->GetOutput());
    extrsurf->Update();

    vtkSmartPointer<vtkPolyData> mesh = extrsurf->GetOutput();
    
    if (targetarea > 0) {
        //post process: fix and uniformly remesh
        char filename[100];
        sprintf(filename, "mymesh_196239462.stl");
        CommonTools::SaveShapeToFile(extrsurf->GetOutput(), filename);

        char fixcmd[500];
        std::cout << "Running admesh.." << std::endl;
        sprintf(fixcmd, "admesh -n -u -f -d %s -b%s -t0.1 -i50", filename, filename);
        int rubbish = system(fixcmd);

        mesh = vtkSmartPointer<vtkPolyData>::Take(
                CommonTools::LoadShapeFromFile(filename));

        remove(filename);

        //remesh
#ifdef USE_VMTK
        vtkSmartPointer<vtkvmtkPolyDataSurfaceRemeshing> remesh =
                vtkSmartPointer<vtkvmtkPolyDataSurfaceRemeshing>::New();
        remesh->SetInputData(mesh);
        remesh->SetNumberOfIterations(10);
        remesh->SetElementSizeModeToTargetArea();
        remesh->SetTargetArea(targetarea);
        remesh->Update();

        mesh->DeepCopy(remesh->GetOutput());
#endif
    }


    //generate normals
    normalgen->SetInputData(mesh);
    normalgen->SplittingOff();
    normalgen->Update();

    CommonTools::SaveShapeToFile(normalgen->GetOutput(), outputGeometry);

    return EXIT_SUCCESS;
}



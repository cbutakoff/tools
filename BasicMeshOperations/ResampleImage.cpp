/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Resample image, smooth, reconstruct shape.
 */
#include <vtkImageData.h>

#include <vtkSmartPointer.h>
#include <vtkImageResample.h>

#include "VTKCommonTools.h"
#include "CommonTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

int main(int argc, char **argv) {
    std::cout << "Resample image. Version 1.0" << std::endl;

    if (argc < 2) {
        std::cout << "Params: " << std::endl;
        std::cout << "-i <image.vtk> \t\t- image (mask)" << std::endl;
        std::cout << "-fx <float> \t\t- resampling factor x (def=1)" << std::endl;
        std::cout << "-fy <float> \t\t- resampling factor y (def=1)" << std::endl;
        std::cout << "-fz <float> \t\t- resampling factor z (def=1)" << std::endl;
        std::cout << "-vx <float> \t\t- spacing for x (def=0.2)" << std::endl;
        std::cout << "-vy <float> \t\t- spacing for y (def=0.2)" << std::endl;
        std::cout << "-vz <float> \t\t- spacing for z (def=0.2)" << std::endl;
        std::cout << "-o <image.vtk> \t\t- optional, to save resampled image" << std::endl;
        return -1;
    }

    char *inimagefile = NULL;
    char *resampledimagefile = NULL;

    float factorx = 1.0;
    float factory = 1.0;
    float factorz = 1.0;
    float vx = 0.2;
    float vy = 0.2;
    float vz = 0.2;


    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-i") == 0) {
            inimagefile = argv[++c];
        } else if (strcmp(argv[c], "-o") == 0) {
            resampledimagefile = argv[++c];
        } else if (strcmp(argv[c], "-fx") == 0) {
            factorx = atof(argv[++c]);
        } else if (strcmp(argv[c], "-fy") == 0) {
            factory = atof(argv[++c]);
        } else if (strcmp(argv[c], "-fz") == 0) {
            factorz = atof(argv[++c]);
        } else if (strcmp(argv[c], "-vx") == 0) {
            vx = atof(argv[++c]);
        } else if (strcmp(argv[c], "-vy") == 0) {
            vy = atof(argv[++c]);
        } else if (strcmp(argv[c], "-vz") == 0) {
            vz = atof(argv[++c]);
        }


    }



    vtkSmartPointer<vtkImageData> imagein = vtkSmartPointer<vtkImageData>::Take(
            CommonTools::LoadImage(inimagefile));


    if (resampledimagefile != NULL) {
        //just resample
        std::cout << "Resampling" << std::endl;

        vtkSmartPointer<vtkImageResample> resamp = vtkSmartPointer<vtkImageResample>::New();
        resamp->SetInputData(imagein);
        resamp->SetAxisMagnificationFactor(0, factorx);
        resamp->SetAxisMagnificationFactor(1, factory);
        resamp->SetAxisMagnificationFactor(2, factorz);
        resamp->SetAxisOutputSpacing(0, vx);
        resamp->SetAxisOutputSpacing(1, vy);
        resamp->SetAxisOutputSpacing(2, vz);
        resamp->SetInterpolationModeToCubic();
        resamp->Update();

        CommonTools::SaveImage(resamp->GetOutput(), resampledimagefile);
    }



    return 0;
}


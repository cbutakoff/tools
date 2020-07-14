/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file 
    \brief Create mask for a given shape with user defined dimensions
 */
#include <vtkPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>

#include <vtkImageData.h>
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkMetaImageWriter.h>
#include <vtkMetaImageReader.h>

#include <vtkSmartPointer.h>

#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "VTKCommonTools.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char **argv) {
    std::cout << "Version 1.1" << std::endl;

    if (argc < 2) {
        std::cout << "Params: " << std::endl;
        std::cout << "-shape <shape.vtk> \t\t- shape" << std::endl;
        std::cout << "-image <image.mhd> \t\t- image, voxels outside the shape are replaced with bg value" << std::endl;
        std::cout << "-out_image <image.vtk|image.mhd> \t\t- resulting image" << std::endl;
        std::cout << "-bg <int> \t\t- value to replace background voxels" << std::endl;
        std::cout << "-fg <int> \t\t- value to replace foreground voxels [optional]. If not specified - foreground is not changed." << std::endl;        
        std::cout << std::endl;
        std::cout << "If -image is not specified, it will be generated around the shape. The following must be specified:" << std::endl;
        std::cout << "-fg <int> \t\t- value for foreground voxels" << std::endl;
        std::cout << "-vx <float> \t\t- voxel dimension x" << std::endl;
        std::cout << "-vy <float> \t\t- voxel dimension y" << std::endl;
        std::cout << "-vz <float> \t\t- voxel dimension z" << std::endl;
        std::cout << "-pad <int> \t\t- number of pixels to pad to each dimension (10 default)" << std::endl;

        return -1;
    }

    char *inshapefile = NULL;
    char *inimagefile = NULL;
    char *outimagefile = NULL;
    //char *outshapefile=NULL;
    unsigned char bg = 255;
    unsigned char fg = 0;
    double vx = 1, vy = 1, vz = 1;
    int padding = 10; //space to add to shape extremes

    bool requested_fg = false;

    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-shape") == 0) {
            inshapefile = argv[++c];
        } else if (strcmp(argv[c], "-image") == 0) {
            inimagefile = argv[++c];
        } else if (strcmp(argv[c], "-out_image") == 0) {
            outimagefile = argv[++c];
        } else if (strcmp(argv[c], "-out_shape") == 0) {
            //outshapefile = argv[++c];
        } else if (strcmp(argv[c], "-bg") == 0) {
            bg = atoi(argv[++c]);
        } else if (strcmp(argv[c], "-fg") == 0) {
            fg = atoi(argv[++c]);
            requested_fg = true;
        } else if (strcmp(argv[c], "-vx") == 0) {
            vx = atof(argv[++c]);
        } else if (strcmp(argv[c], "-vy") == 0) {
            vy = atof(argv[++c]);
        } else if (strcmp(argv[c], "-vz") == 0) {
            vz = atof(argv[++c]);
        } else if (strcmp(argv[c], "-pad") == 0) {
            padding = atoi(argv[++c]);
        }

    }

    vtkSmartPointer<vtkPolyData> inshape = vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(inshapefile));

    double origin[3];
    double spacing[3];
    int wextent[6];


    if (inimagefile != NULL)
        std::cout << "Image: " << inimagefile << std::endl;
    else
        std::cout << "Image: will be generated" << std::endl;

    std::cout << "Shape: " << inshapefile << std::endl;
    std::cout << "Output image: " << outimagefile << std::endl;
    //cout<<"Output shape: "<<outshapefile<<endl;
    std::cout << "Foreground: " << int(fg) << std::endl;
    std::cout << "Background: " << int(bg) << std::endl;
    std::cout << "Voxel: " << vx << " x " << vy << " x " << vz << std::endl;
    std::cout << "Padding (pix): " << padding<< std::endl << std::endl;


    vtkSmartPointer<vtkImageData> inimage = vtkSmartPointer<vtkImageData>::New();
    if (inimagefile != NULL) {
        std::cout << "Creating mask from a supplied image" << std::endl;
        vtkSmartPointer<vtkMetaImageReader> imread = vtkSmartPointer<vtkMetaImageReader>::New();
        imread->SetFileName(inimagefile);
        imread->Update();
        inimage->DeepCopy(dynamic_cast<vtkImageData*> (imread->GetOutput()));

        inimage->GetOrigin(origin);
        inimage->GetSpacing(spacing);
        inimage->GetExtent(wextent);

        if (requested_fg) {
            int* dims = inimage->GetDimensions();

            // Fill every entry of the image data with "2.0"
            for (int z = 0; z < dims[2]; z++) {
                for (int y = 0; y < dims[1]; y++) {
                    for (int x = 0; x < dims[0]; x++) {
                        inimage->SetScalarComponentFromDouble(x, y, z, 0, fg);
                    }
                }
            }
        }


    } else {
        std::cout << "Creating mask and generating an image" << std::endl;

        spacing[0] = vx;
        spacing[1] = vy;
        spacing[2] = vz;

        double bounds[6];
        inshape->GetBounds(bounds);
        double t[3];
        t[0] = -bounds[0] + padding * spacing[0];
        t[1] = -bounds[2] + padding * spacing[1];
        t[2] = -bounds[4] + padding * spacing[2];

        //vtkSmartPointer<vtkTransform> tform = vtkSmartPointer<vtkTransform>::New();
        //tform->Translate(t);

        //vtkSmartPointer<vtkTransformPolyDataFilter> tf = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        //tf->SetInputData(inshape);
        //tf->SetTransform(tform);
        //tf->Update();

        //inshape->DeepCopy(tf->GetOutput());

        //if(outshapefile!=NULL) 
        //	blShapeUtils::ShapeUtils::SaveShapeToFile(inshape,outshapefile);

        origin[0] = -t[0];
        origin[1] = -t[1];
        origin[2] = -t[2];
        int dims[3];
        inshape->GetBounds(bounds);
        dims[0] = floor((bounds[1] - bounds[0]) / spacing[0]) + padding * 2;
        dims[1] = floor((bounds[3] - bounds[2]) / spacing[1]) + padding * 2;
        dims[2] = floor((bounds[5] - bounds[4]) / spacing[2]) + padding * 2;

        //inimage->SetScalarTypeToUnsignedChar();
        inimage->SetOrigin(origin);
        inimage->SetSpacing(spacing);
        inimage->SetDimensions(dims);
        //inimage->Update();
        inimage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
        inimage->GetExtent(wextent);
        memset(inimage->GetScalarPointer(), fg, dims[0] * dims[1] * dims[2] * sizeof (unsigned char));
        std::cout << "Dimensions: " << dims[0] << ", " << dims[1] << ", " << dims[2] << std::endl;
    }

    vtkSmartPointer<vtkPolyDataToImageStencil> filt =
            vtkSmartPointer<vtkPolyDataToImageStencil>::New();


    filt->SetOutputOrigin(origin);
    filt->SetOutputSpacing(spacing);
    filt->SetOutputWholeExtent(wextent);
    filt->SetInputData(inshape);
    filt->Update();

    vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
    stencil->SetStencilData(filt->GetOutput());
    stencil->SetBackgroundValue(bg);
    stencil->SetInputData(inimage);
    stencil->Update();


    if (outimagefile[ strlen(outimagefile)-1 ] == 'k')  //vtk
    {
        vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
        writer->SetInputData(stencil->GetOutput());
        writer->SetFileTypeToBinary();
        writer->SetFileName(outimagefile);
        writer->Update();
    }
    else
    {
        std::string filePath = string(outimagefile)+string(".mhd");
        std::string filePathRaw =  string(outimagefile)+string(".raw");
        vtkSmartPointer<vtkMetaImageWriter> writer =
          vtkSmartPointer<vtkMetaImageWriter>::New();
        writer->SetInputData(stencil->GetOutput());
        writer->SetFileName(filePath.c_str());
        writer->SetRAWFileName(filePathRaw.c_str());
        writer->Write();
    }

    return 0;
}

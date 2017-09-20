/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageCast.h>
#include <vtkPolyDataWriter.h>
#include <vtkMarchingCubes.h>

int main(int argc, char** argv)
{
    std::cout << "Version 1.0"<< std::endl;

    if( argc<2 ) 
    {
            std::cout << "Params: "<< std::endl;
            std::cout << "-i <image.vtk> \t\t- input"<< std::endl;
            std::cout << "-o <mesh.vtk> \t\t- output"<< std::endl;
            std::cout << "-s <float> \t\t- smoothing sigma"<< std::endl;
            std::cout << "-l <float> \t\t- isolevel"<< std::endl;
            return EXIT_FAILURE;
    }

    char *inputimagefile=NULL;
    char *outputshapefile=NULL;
    float sigma = 1;
    float level = 0.5;


    for(int c=1; c<argc; c++)
    {
            if( strcmp(argv[c],"-i")==0 )
            {
                    inputimagefile = argv[++c];
            }
            else if( strcmp(argv[c],"-o")==0 )
            {
                    outputshapefile = argv[++c];
            }
            else if( strcmp(argv[c],"-s")==0 ) //open value
            {
                    sigma = atof(argv[++c]);
            }
            else if( strcmp(argv[c],"-l")==0 ) //open value
            {
                    level = atof(argv[++c]);
            }
    }

    
    std::cout<<"Input image: "<<inputimagefile<<std::endl;
    std::cout<<"Output shape: "<<outputshapefile<<std::endl;
    std::cout<<"Gaussian sigma: "<<sigma<<std::endl;
    std::cout<<"Isosurface level: "<<level<<std::endl;
    
    
    vtkSmartPointer<vtkDataSetReader> rdr = vtkSmartPointer<vtkDataSetReader>::New();
    rdr->SetFileName(inputimagefile);
    rdr->Update();
    
    std::cout<<"Casting image to float"<<std::endl;
    vtkSmartPointer<vtkImageCast> caster = vtkSmartPointer<vtkImageCast>::New();
    caster->SetInputData(rdr->GetOutput());
    caster->SetOutputScalarTypeToFloat();
    caster->Update();
    
    std::cout<<"Smoothing the image"<<std::endl;
    vtkSmartPointer<vtkImageGaussianSmooth> smooth = vtkSmartPointer<vtkImageGaussianSmooth>::New();
    smooth->SetInputData((vtkDataObject *)caster->GetOutput());
    smooth->SetStandardDeviation(sigma);
    smooth->Update();
    
    std::cout<<"Running marching cubes"<<std::endl;
    vtkSmartPointer<vtkMarchingCubes> mc = vtkSmartPointer<vtkMarchingCubes>::New();
    mc->SetInputData((vtkDataObject *)smooth->GetOutput());
    mc->SetNumberOfContours(1);
    mc->SetValue(0,level);
    mc->Update();
    
    std::cout<<"Saving"<<std::endl;
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetFileName(outputshapefile);
    wr->SetFileTypeToBinary();
    wr->SetInputData(mc->GetOutput());
    wr->Write();
           
}


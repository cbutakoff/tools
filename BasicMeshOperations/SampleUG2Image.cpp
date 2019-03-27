/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/
// Given a mesh (unstructured grid) and and images , 
// create a pointdata array by sampling the image with image values

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellLocator.h>
#include <vtkCell.h>
#include <vtkCellData.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkPoint.h>

#include <iostream>
#include <limits>


int main(int argc, char** argv)
{
    if(argc<4)
    {
        std::cout<<"Usage: SampleUG2Image mesh.vtk image.nrrd out_image.nrrd array_name component_number"<<std::endl;
        std::cout<<"pixeltype is uchar"<<std::endl;
        std::cout<<"The image is sampled at vertex positions"<<std::endl;
        std::cout<<"!! Does lazy interpolation: assigns value from the first point of the cell to the pixel"<<std::endl;
        
        exit(-1);
    }

    
    int c=1;
    const char* meshfilename = argv[c++];
    const char* imagefilename = argv[c++];
    const char* outfilename = argv[c++];
    const char * arrayname = argv[c++];
    const int component = atoi(argv[c++]);
    const int Dimension = 3; //image dimensionality
    
    std::cout<<"Mesh: "<<meshfilename<<std::endl;
    std::cout<<"Output: "<<outfilename<<std::endl;
    std::cout<<"Image: "<<imagefilename<<std::endl;
    std::cout<<"Array: "<<arrayname<<std::endl;
    std::cout<<"Component: "<<component<<std::endl;
    
    
    
    
    vtkSmartPointer<vtkDataSetReader> mesh_reader = vtkSmartPointer<vtkDataSetReader>::New();
    std::cout<<"Reading mesh: "<<meshfilename<<std::endl;
    mesh_reader->SetFileName(meshfilename);
    mesh_reader->Update();
    
    vtkUnstructuredGrid* mesh = vtkUnstructuredGrid::SafeDownCast(mesh_reader->GetOutput());

    auto array = mesh->GetPointData()->GetArray(arrayname);
    if( array==nullptr )
    {
        std::cout<<"Pointdata array not present: "<<arrayname<<std::endl;
        std::cout<<"Available point arrays:"<<std::endl;
        for(int i=0; i<mesh->GetPointData()->GetNumberOfArrays(); i++)
        {
            std::cout<<"["<<mesh->GetPointData()->GetArray(i)->GetName()<<"]"<<std::endl;            
        }

        std::cout<<"Available cell arrays:"<<std::endl;
        for(int i=0; i<mesh->GetCellData()->GetNumberOfArrays(); i++)
        {
            std::cout<<"["<<mesh->GetCellData()->GetArray(i)->GetName()<<"]"<<std::endl;            
        }
        
        throw;
    }
    
    double * bounds = mesh->GetBounds(); //(xmin,xmax, ymin,ymax, zmin,zmax)
    
    
    //main program
    std::cout<<"Reading image: "<<imagefilename<<std::endl;
    
    
    typedef unsigned char InputPixelType;
    typedef float OutputPixelType;
    typedef itk::Image< InputPixelType, Dimension > InputImageType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    reader->SetFileName(imagefilename);
    writer->SetFileName(outfilename);
    reader->Update();
    InputImageType::Pointer image = reader->GetOutput();
    

    InputImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    const long int npixels = size[0]*size[1]*size[2];
    
    OutputImageType::Pointer output_image = OutputImageType::New();
    output_image->SetSpacing( image->GetSpacing() );
    output_image->SetOrigin( image->GetOrigin() );
    output_image->SetRegions( image->GetLargestPossibleRegion() );
    output_image->Allocate();
    
    //build locator
    vtkSmartPointer<vtkCellLocator> cl =     vtkSmartPointer<vtkCellLocator>::New();
    cl->SetDataSet(mesh);
    cl->BuildLocator();
    
    
    itk::ImageRegionIteratorWithIndex<InputImageType> imageIterator(image, image->GetLargestPossibleRegion());    
    itk::ImageRegionIterator<OutputImageType> imageIteratorOut(output_image, output_image->GetLargestPossibleRegion());    
    
    long int cc=0;
    while(!imageIterator.IsAtEnd())
    {
        if (cc%20000 == 0)
            std::cout<<"Sampling progress: "<<100*cc/npixels<<"%\r";
            
        itk::Point<double, Dimension> pt;
        image->TransformIndexToPhysicalPoint( imageIterator.GetIndex(), pt);
    
        if( (pt[0]>=bounds[0] && pt[0]<=bounds[1]) &&
           (pt[1]>=bounds[2] && pt[1]<=bounds[3]) &&
           (pt[2]>=bounds[4] && pt[2]<=bounds[5]) )
        {
            long int cellid = cl->FindCell(pt.Begin());
            
            //std::cout<<"cellid "<<cellid<<std::endl;
            if (cellid<0) imageIteratorOut.Set( std::numeric_limits<double>::max() );
            else
            {
                long int  pointid = mesh->GetCell(cellid)->GetPointId(0);
                double v = array->GetComponent(pointid, component);
                imageIteratorOut.Set(v);
//                std::cout<<"cell: "<<cellid<<", v="<<v<<std::endl;
            }
        }
        else
            imageIteratorOut.Set( std::numeric_limits<double>::max() );
            
        ++imageIteratorOut;
        ++imageIterator;
        ++cc;
    }
    std::cout<<std::endl;
    
    std::cout<<"Saving image: "<<outfilename<<std::endl;
    writer->SetInput(output_image);
    writer->Update();
}



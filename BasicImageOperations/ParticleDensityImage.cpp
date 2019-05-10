/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>

typedef float    float32_t;
typedef double float64_t;



using namespace std;

int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Fill  image voxels with the particle density (0-100 percent) (3d histogram). Points must be in vtp format" << std::endl;
        std::cerr << "Both input and output images have uint8 (unsigned char) voxel type" << std::endl;
        std::cerr << "Usage: ParticleDensityImage inputImageFile outputImageFile points.vtp" << std::endl;
        return EXIT_FAILURE;
    }

    //poarsing parameters
    int c = 1;
    const char* imagefilename = argv[c++];
    const char* outputfilename = argv[c++];
    const char* pointsfile = argv[c++];


    const unsigned int Dimension = 3;

    //main program
    typedef uint8_t InputPixelType;
    typedef float32_t OutputPixelType;
    typedef itk::Image< InputPixelType, Dimension > InputImageType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    typedef itk::Point<double, Dimension> PointType;


    OutputImageType::Pointer outputImage = OutputImageType::New();
    {
        cout<<"Reading image"<<imagefilename<<endl;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(imagefilename);
        reader->Update();
        InputImageType * inputImage = reader->GetOutput();

        cout<<"Creating output image"<<endl;

        outputImage->SetOrigin( inputImage->GetOrigin() );
        outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
        outputImage->SetDirection( inputImage->GetDirection() );
        outputImage->SetSpacing( inputImage->GetSpacing() );
        outputImage->Allocate();

        typedef itk::ImageRegionIterator<OutputImageType> ImageIteratorType;
        ImageIteratorType it(outputImage, outputImage->GetLargestPossibleRegion());
        for( it = it.Begin(); !it.IsAtEnd(); ++it)
        {
            it.Set(0); //initialize image to 0
        }

    }

    cout<<"Reading points: "<<pointsfile<<endl;
    vtkSmartPointer<vtkXMLPolyDataReader> ptsrdr =     vtkSmartPointer<vtkXMLPolyDataReader>::New();
    ptsrdr->SetFileName(pointsfile);
    ptsrdr->Update();

    uint64_t npoints = ptsrdr->GetOutput()->GetNumberOfPoints();



    for( uint64_t i=0; i<npoints; i++ )
    {
        PointType pt;
        ptsrdr->GetOutput()->GetPoint(i, pt.GetDataPointer());

        OutputImageType::IndexType index;
        outputImage->TransformPhysicalPointToIndex(pt, index);

        outputImage->SetPixel( index, outputImage->GetPixel(index)+double(1.0)/double(npoints) );
    }



    typedef itk::ImageRegionIterator<OutputImageType> ImageIteratorType;
    ImageIteratorType it(outputImage, outputImage->GetLargestPossibleRegion());
    for( it = it.Begin(); !it.IsAtEnd(); ++it)
    {
        if(it.Get()>0)
            it.Set( -log(it.Get()) ); //initialize image to 0
    }
    

    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputfilename);
    writer->SetInput(outputImage);

    //  The output of the resampling filter is connected to a writer and the
    //  execution of the pipeline is triggered by a writer update.
    cout<<"Writing to "<<outputfilename<<endl;
    try {
        writer->Update();
    } catch (itk::ExceptionObject & excep) {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
    }

    

    
    return EXIT_SUCCESS;
}


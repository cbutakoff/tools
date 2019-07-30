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
#include "ImageRegionIterator.h"

#include <set>

using namespace std;

int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "From an existing labelmap, cteate a new labelmap with only specified labels. 
        Expected datatype: uint8. Pixels with the indicated labes are set to 255, everythign else to 0." << std::endl;
        std::cerr << "Usage: ExtractLabels inputImageFile outputImageFile [list of labels]" << std::endl;
        std::cerr << "Example: ExtractLabels inputImageFile.nrrd outputImageFile.nrrd 1 2 5 8" << std::endl;
        return EXIT_FAILURE;
    }

    //poarsing parameters
    int c = 1;
    const char* imagefilename = argv[c++];
    const char* outputfilename = argv[c++];

    set<int> output_labels;
    for(;c<argc; c++)
    {
        output_labels.insert( atoi(argc[c] ) );
    }


    const unsigned int Dimension = 3;
    typedef itk::AffineTransform< double, Dimension > TransformType;
    TransformType::OutputVectorType axis;
    axis[0] = atof(argv[c++]);
    axis[1] = atof(argv[c++]);
    axis[2] = atof(argv[c++]);

    const double angle = atof(argv[c++]);
    const int padding = atoi(argv[c++]);

    std::cout << "Input: " << imagefilename << std::endl;
    std::cout << "Output: " << outputfilename << std::endl;
    std::cout << "Labels: ";
    for( lbl:output_labels )
        cout<<lbl<<", ";
    cout<<endl;
    


    //main program
    typedef uint8_t InputPixelType;
    typedef uint8_t OutputPixelType;
    typedef itk::Image< InputPixelType, Dimension > InputImageType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    reader->SetFileName(imagefilename);
    writer->SetFileName(outputfilename);
    reader->Update();


    typedef itk::ImageRegionIterator<InputImageType> InIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutIteratorType;

    

    const InputImageType * inputImage = reader->GetOutput();


    OutputImageType::Pointer outImage = OutputImageType::New();
    outImage->SetSpacing( inputImage->GetSpacing() );
    outImage->SetRegions( inputImage->GetLargestPossibleRegion() );
    outImage->SetOrigin( inputImage->GetOrigin() );
    outImage->Allocate();


    InputImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
    uint64_t npixels = 1;
    for( int i=0; i<size.GetSizeDimension(); i++ ) npixels *= size[i];


    InImageIterator  itin( inputImage, inputImage->GetLargestPossibleRegion() );
    InImageIterator  itout( outImage, outImage->GetLargestPossibleRegion() );
    itin.GoToBegin();
    itout.GoToBegin();

    uint64_t c=0;
    while( !itin.IsAtEnd() )
    {
        if( c % 100000 == 0 )
            cout<<"Progress "<<c*100/npixels<<"\r";


        uint8_t v = itin.Get();

        

        auto search = output_labels.find( v );
        if (search != output_labels.end()) {
            itout.Set(255);
        } else {
            itout.Set(0);
        }        

        ++itin;
        ++itout;
        ++c;
    }
    cout<<endl;



    try {
        writer->SetInput( outImage );
        writer->Update();
    } catch (itk::ExceptionObject & excep) {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
    }

    

    
    return EXIT_SUCCESS;
}


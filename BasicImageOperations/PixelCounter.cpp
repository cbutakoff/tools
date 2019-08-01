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
#include "itkImageRegionIterator.h"


using namespace std;

int main(int argc, char** argv) {

    if (argc < 2) {
        std::cerr << "Count nonzero pixels. Image should be uint8." << std::endl;
        std::cerr << "Usage: PixelCounter inputImageFile" << std::endl;
        return EXIT_FAILURE;
    }

    //poarsing parameters
    int c = 1;
    const char* imagefilename = argv[c++];



    const unsigned int Dimension = 3;

    std::cout << "Input: " << imagefilename << std::endl;
    


    //main program
    typedef uint8_t InputPixelType;
    typedef itk::Image< InputPixelType, Dimension > InputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(imagefilename);
    reader->Update();


    typedef itk::ImageRegionIterator<InputImageType> InIteratorType;

    

    InputImageType * inputImage = reader->GetOutput();



    InputImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
    uint64_t npixels = 1;
    for( int i=0; i<size.GetSizeDimension(); i++ ) npixels *= size[i];


    InIteratorType  itin( inputImage, inputImage->GetLargestPossibleRegion() );
    itin.GoToBegin();

    uint64_t count=0;
    uint64_t pixel_counter = 0;
    while( !itin.IsAtEnd() )
    {
        if( count % 100000 == 0 )
            cout<<"Progress "<<count *100/npixels<<"\r";


        if ( itin.Get() != 0 )  ++pixel_counter;

        
        ++itin;
        ++count;
    }
    cout<<endl;

    cout<<"Number of nonzero pixels: "<<pixel_counter<<endl;


    
    return EXIT_SUCCESS;
}


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
#include "itkImageRegionIterator.h"

#include <set>
#include <string>

using namespace std;


int main(int argc, char** argv) {

    if (argc < 3) {
        std::cerr << "Threshold input ikmages and combine labels"<< std::endl;
        std::cerr << "Usage: CombineImageLabels -i image1 threshold1 label1 -i image2 threshold2 label2 ... -o output_image" << std::endl;
        return EXIT_FAILURE;
    }

    //poarsing parameters
    typedef struct __image_info {
        string imagefilename;
        double threshold;
        int label;
    } ImageInfo;

    const char* outputfilename = nullptr;

    vector<ImageInfo> input_images;

    for(int c=1; c<argc;)    {
        if( strcmp(argv[c],"-i")==0 ){
            ImageInfo info;
            info.imagefilename = string(argv[++c]);
            info.threshold = atof(argv[++c]);
            info.label = atoi(argv[++c]);
            input_images.push_back(info);
        }
        else if( strcmp(argv[c],"-o")==0 ){
            outputfilename = argv[++c];
        }
        c++;
    }



    const unsigned int Dimension = 3;

    for( auto& info: input_images) {
        std::cout << "Input: " << info.imagefilename <<", t: "<<info.threshold<<", l: "<<info.label<< std::endl;
    }

    std::cout << "Output: " << outputfilename << std::endl;
    


    //main program
    typedef int16_t InputPixelType;
    typedef uint8_t OutputPixelType;
    typedef itk::Image< InputPixelType, Dimension > InputImageType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    reader->SetFileName(input_images[0].imagefilename);
    writer->SetFileName(outputfilename);

    cout<<"Reading image"<<endl;
    reader->Update();


    typedef itk::ImageRegionIterator<InputImageType> InIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType> OutIteratorType;

    

    InputImageType * inputImage = reader->GetOutput();


    OutputImageType::Pointer outImage = OutputImageType::New();
    outImage->SetSpacing( inputImage->GetSpacing() );
    outImage->SetRegions( inputImage->GetLargestPossibleRegion() );
    outImage->SetOrigin( inputImage->GetOrigin() );
    outImage->Allocate();


    InputImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
    uint64_t npixels = 1;
    for( int i=0; i<size.GetSizeDimension(); i++ ) npixels *= size[i];

    int image_no=0;
    for( auto& image_info: input_images ){
        cout<<"Processing image "<<image_info.imagefilename<<endl;
        reader->SetFileName( image_info.imagefilename );
        reader->Update();
        inputImage = reader->GetOutput();

        InIteratorType  itin( inputImage, inputImage->GetLargestPossibleRegion() );
        OutIteratorType  itout( outImage, outImage->GetLargestPossibleRegion() );
        itin.GoToBegin();
        itout.GoToBegin();

        uint64_t count=0;
        while( !itin.IsAtEnd() )
        {
            if( count % 1000000 == 0 )
                cout<<"Progress "<<count *100/npixels<<"% \r";


            auto v = itin.Get();
            if( v>image_info.threshold )
                itout.Set(image_info.label);
            else
                if(image_no==0)
                    itout.Set(0);

    
            ++itin;
            ++itout;
            ++count;
        }
        cout<<endl;
        image_no++;
    }

    cout<<"saving the result"<<endl;
    try {
        writer->SetInput( outImage );
        writer->Update();
    } catch (itk::ExceptionObject & excep) {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
    }

    

    
    return EXIT_SUCCESS;
}


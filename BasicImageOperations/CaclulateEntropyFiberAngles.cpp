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
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkImageRegionIterator.h>

#include <Eigen/Dense>


int main(int argc, char** argv) {

    if (argc < 5) {
        std::cerr << "From the 3 images representing fiber orientation estimate angles and calculate angle homogeneity in the neighborhoods." << std::endl;
        std::cerr << "Usage: CalculateEntropyFiberAngles coordx_image coordy_image coordz_image out_image neigborhood_radius" << std::endl;
        return EXIT_FAILURE;
    }

    //parsing parameters
    int c = 1;
    const char* imagex_filename = argv[c++];
    const char* imagey_filename = argv[c++];
    const char* imagez_filename = argv[c++];
    const char* output_filename = argv[c++];
    int neighborhood_radius = atoi(argv[c++]);


    const unsigned int Dimension = 3;

    std::cout << "Input x: " << imagex_filename << std::endl;
    std::cout << "Input y: " << imagey_filename << std::endl;
    std::cout << "Input z: " << imagez_filename << std::endl;
    std::cout << "Output: " << output_filename << std::endl;
    std::cout << "Neighborhood radius: " << neighborhood_radius << std::endl;



    //main program
    typedef float InputPixelType;
    typedef float OutputPixelType;
    typedef itk::Image< InputPixelType, Dimension > InputImageType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    ReaderType::Pointer readerx = ReaderType::New();
    ReaderType::Pointer readery = ReaderType::New();
    ReaderType::Pointer readerz = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    readerx->SetFileName(imagex_filename);
    readery->SetFileName(imagex_filename);
    readerz->SetFileName(imagex_filename);
    writer->SetFileName(output_filename);



    std::cout<<"Reading X"<<std::endl;
    readerx->Update();
    std::cout<<"Reading Y"<<std::endl;
    readery->Update();
    std::cout<<"Reading Z"<<std::endl;
    readerz->Update();

    
    const InputImageType * inputx = readerx->GetOutput();
    const InputImageType * inputy = readery->GetOutput();
    const InputImageType * inputz = readerz->GetOutput();


    //create the output image
    OutputImageType::Pointer image = OutputImageType::New();
    image->SetRegions(inputx->GetLargestPossibleRegion());
    image->Allocate();
    image->FillBuffer(0);

    
    InputImageType::SizeType radius;
    radius[0] = neighborhood_radius;
    radius[1] = neighborhood_radius;
    radius[2] = neighborhood_radius;
    
    itk::ConstNeighborhoodIterator<InputImageType> 
            iteratorx(radius, inputx, image->GetRequestedRegion());
    itk::ConstNeighborhoodIterator<InputImageType> 
            iteratory(radius, inputy, image->GetRequestedRegion());
    itk::ConstNeighborhoodIterator<InputImageType> 
            iteratorz(radius, inputz, image->GetRequestedRegion());

    itk::ImageRegionIterator<OutputImageType>  it( image, image->GetRequestedRegion() );
    
    iteratorx.GoToBegin();
    iteratory.GoToBegin();
    iteratorz.GoToBegin();
    it.GoToBegin();
  
    const unsigned int nneighbors = (neighborhood_radius*2+1)*(neighborhood_radius*2+1)*(neighborhood_radius*2+1);
    Eigen::MatrixXf neighbors( nneighbors-1, 3 ); //to store the vectors
    Eigen::Vector3f central_vector;
    
    int c_id = (int) (nneighbors / 2); // get offset of center pixel
    
    while ( ! iteratorx.IsAtEnd() )
    {
        unsigned int j = 0;
        for (unsigned int i = 0; i < nneighbors; ++i)
        {    
            if(i!=c_id) //exclude the center
            {
                neighbors(j, 0) = iteratorx.GetPixel(i);  
                neighbors(j, 1) = iteratory.GetPixel(i);  
                neighbors(j, 2) = iteratorz.GetPixel(i);  
                j++;
            }
            else
            {
                central_vector(0) = iteratorx.GetPixel(i);  
                central_vector(1) = iteratory.GetPixel(i);  
                central_vector(2) = iteratorz.GetPixel(i);                  
            }
        }
        
        // calculate the angles 
        Eigen::VectorXf dp = neighbors*central_vector;
        
        const float variance = ((dp.array()-dp.mean()).pow(2)).mean();
        it.Set(variance);        
        
        ++iteratorx;
        ++iteratory;
        ++iteratorz;
        ++it;
    }    
    
    //  The output of the resampling filter is connected to a writer and the
    //  execution of the pipeline is triggered by a writer update.
    try {
        writer->SetInput(image);
        writer->Update();
    } catch (itk::ExceptionObject & excep) {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
    }

    

    
    return EXIT_SUCCESS;
}


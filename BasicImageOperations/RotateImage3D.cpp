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
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkConstantPadImageFilter.h"

int main(int argc, char** argv) {

    if (argc < 8) {
        std::cerr << "Rotate image about its center (uses linear interpolation)" << std::endl;
        std::cerr << "Usage: RotateImage3D inputImageFile outputImageFile axis_X axis_Y axis_Z angle_radians padding_pix" << std::endl;
        return EXIT_FAILURE;
    }

    //poarsing parameters
    int c = 1;
    const char* imagefilename = argv[c++];
    const char* outputfilename = argv[c++];


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
    std::cout << "Axis: " << axis << std::endl;
    std::cout << "Angle(rad): " << angle << std::endl;
    std::cout << "Padding(pix): " << padding << std::endl;



    //main program
    typedef unsigned char InputPixelType;
    typedef unsigned char OutputPixelType;
    typedef itk::Image< InputPixelType, Dimension > InputImageType;
    typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    ReaderType::Pointer reader = ReaderType::New();
    WriterType::Pointer writer = WriterType::New();
    reader->SetFileName(imagefilename);
    writer->SetFileName(outputfilename);


    typedef itk::ResampleImageFilter<InputImageType, OutputImageType > FilterType;
    FilterType::Pointer filter = FilterType::New();


    reader->Update();

    

    const InputImageType * inputImage = reader->GetOutput();


    // Software Guide : EndCodeSnippet
    typedef itk::LinearInterpolateImageFunction<
            InputImageType, double > InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    filter->SetInterpolator(interpolator);
    filter->SetDefaultPixelValue(0);



    const InputImageType::SpacingType & spacing = inputImage->GetSpacing();
    InputImageType::PointType origin = inputImage->GetOrigin();
    InputImageType::SizeType size =
            inputImage->GetLargestPossibleRegion().GetSize();


    InputImageType::PointType origin_new;
    InputImageType::SizeType size_new;
    
    origin_new[0] = origin[0]-padding*spacing[0];
    origin_new[1] = origin[1]-padding*spacing[1];
    origin_new[2] = origin[2]-padding*spacing[2];
    size_new[0] = size[0]+padding*2;
    size_new[1] = size[1]+padding*2;
    size_new[2] = size[2]+padding*2;
    
    std::cout << "Origin: " << origin << std::endl;
    std::cout << "Size: " << size << std::endl;

   std::cout << "Output Origin: " << origin_new << std::endl;
    std::cout << "Output Size: " << size_new << std::endl;    
    
    filter->SetOutputOrigin(origin_new);
    filter->SetOutputSpacing(spacing);
    filter->SetOutputDirection(inputImage->GetDirection());
    filter->SetSize(size_new);

    filter->SetInput(inputImage);
    writer->SetInput(filter->GetOutput());

    //translate to the center
    TransformType::Pointer transform = TransformType::New();

    TransformType::OutputVectorType translation1;
    const double imageCenterX = origin[0] + spacing[0] * size[0] / 2.0;
    const double imageCenterY = origin[1] + spacing[1] * size[1] / 2.0;
    const double imageCenterZ = origin[2] + spacing[2] * size[2] / 2.0;
    translation1[0] = -imageCenterX;
    translation1[1] = -imageCenterY;
    translation1[2] = -imageCenterZ;
    transform->Translate(translation1);
    // Software Guide : EndCodeSnippet
    std::cout << "imageCenterX = " << imageCenterX << std::endl;
    std::cout << "imageCenterY = " << imageCenterY << std::endl;
    std::cout << "imageCenterZ = " << imageCenterZ << std::endl;

    //  Rotate


    transform->Rotate3D(axis, angle, false);

    //  The third and final step requires translating the image origin back to
    //  its previous location. This is be done by applying a translation equal
    //  to the origin values.
    TransformType::OutputVectorType translation2;
    translation2[0] = imageCenterX;
    translation2[1] = imageCenterY;
    translation2[2] = imageCenterZ;
    transform->Translate(translation2, false);
    filter->SetTransform(transform);

    //  The output of the resampling filter is connected to a writer and the
    //  execution of the pipeline is triggered by a writer update.
    try {
        writer->Update();
    } catch (itk::ExceptionObject & excep) {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
    }

    

    
    return EXIT_SUCCESS;
}


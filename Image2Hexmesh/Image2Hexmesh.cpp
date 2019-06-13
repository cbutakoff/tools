#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageIOBase.h"


template <unsigned int VDimension>
int ReadScalarImage(const char *inputFileName,
                    const itk::ImageIOBase::IOComponentType componentType);

template <class TImage>
int ReadImage(const char *fileName,
              typename TImage::Pointer image);


template <class TImage>
int CreateHexMeshFromImage(const char *outputFileName,
              typename TImage::Pointer image);


int main(int argc, char *argv[])
{
    int nthreads, tid;

    if (argc < 2)
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0];
        std::cerr << " <InputImage.mhd> <OutputMesh.vtu>"<<std::endl;
        std::cerr << "Value 0 in the image is reserved for the background"<<std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }

    int c = 1;
    const char *inputFileName = argv[c++];
    const char *outputFileName = argv[c++];

    //read image
    std::cout<<"Reading image "<<inputFileName<<std::endl;
    itk::ImageIOBase::Pointer imageIO =
        itk::ImageIOFactory::CreateImageIO(
            inputFileName,
            itk::ImageIOFactory::ReadMode);

    imageIO->SetFileName(inputFileName);
    imageIO->ReadImageInformation();

    using IOPixelType = itk::ImageIOBase::IOPixelType;
    const IOPixelType pixelType = imageIO->GetPixelType();

    std::cout << "Pixel Type is "
              << itk::ImageIOBase::GetPixelTypeAsString(pixelType)
              << std::endl;

    using IOComponentType = itk::ImageIOBase::IOComponentType;
    const IOComponentType componentType = imageIO->GetComponentType();

    std::cout << "Component Type is "
              << imageIO->GetComponentTypeAsString(componentType)
              << std::endl;

    const unsigned int imageDimension = imageIO->GetNumberOfDimensions();

    std::cout << "Image Dimension is " << imageDimension << std::endl;

    switch (pixelType)
    {
    case itk::ImageIOBase::SCALAR:
    {
        if (imageDimension == 3)
        {
            return ReadScalarImage<3>(inputFileName, componentType);
        }
        else 
        {
            std::cout<<"Only 3D images are supported"<<std::endl;
            return EXIT_FAILURE;
        }
    }

    default:
        std::cerr << "Nonscalar pixel types are not implemented yet!" << std::endl;
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;

}

template <class TImage>
int ReadImage(const char *fileName,
              typename TImage::Pointer image)
{
    using ImageType = TImage;
    using ImageReaderType = itk::ImageFileReader<ImageType>;

    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    try
    {
        reader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    image->Graft(reader->GetOutput());

    return EXIT_SUCCESS;
}

template <unsigned int VDimension>
int ReadScalarImage(const char *inputFileName, const char *outputFileName
                    const itk::ImageIOBase::IOComponentType componentType)
{
    switch (componentType)
    {
    default:
    case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        std::cerr << "Unknown and unsupported component type!" << std::endl;
        return EXIT_FAILURE;

    case itk::ImageIOBase::UCHAR:
    {
        using PixelType = unsigned char;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::CHAR:
    {
        using PixelType = char;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::USHORT:
    {
        using PixelType = unsigned short;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::SHORT:
    {
        using PixelType = short;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::UINT:
    {
        using PixelType = unsigned int;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::INT:
    {
        using PixelType = int;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::ULONG:
    {
        using PixelType = unsigned long;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::LONG:
    {
        using PixelType = long;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::FLOAT:
    {
        using PixelType = float;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }

    case itk::ImageIOBase::DOUBLE:
    {
        using PixelType = double;
        using ImageType = itk::Image<PixelType, VDimension>;

        typename ImageType::Pointer image = ImageType::New();

        if (ReadImage<ImageType>(inputFileName, image) == EXIT_FAILURE)
        {
            return EXIT_FAILURE;
        }

        std::cout << image << std::endl;
        CreateHexMeshFromImage(outputFileName, image);
        break;
    }
    }


    return EXIT_SUCCESS;
}




template <class TImage>
int CreateHexMeshFromImage(const char *outputFileName,
              typename TImage::Pointer image)
{
    using IteratorType = itk::ImageRegionConstIteratorWithIndex<TImage>;
    using PointType = TImage::PointType;

    using idtype = int64_t;

    auto spacing = image->GetSpacing();
    auto size = image->GetLargestPossibleRegion().GetSize();

    const idtype sx = size[0];
    const idtype sy = size[1];
    const idtype sz = size[2];


    std::cout<<"Creating elements"<<std::endl;

    vtkSmartPointer<vtkCellArray> hexas =     vtkSmartPointer<vtkCellArray> ::New();

    //create a mask to know which elements exist
    const idtype nvertices_max = sx*sy*(sz+1) + sy*(sy+1) + sx+1 +1;
    char* mask = new char[ nvertices_max ];
    memset(mask, 0, nvertices_max);

    IteratorType it(image, image->GetLargestPossibleRegion());
    it.Begin();

    while (!it.IsAtEnd())
    {
        auto label = it.Get();

        if (label!=0) 
        {
            auto index = it.GetIndex();

            //the order in which the vetices of the mesh will be generated : cols(x), then rows(y), then slices (z)
            const idtype x = index[0];
            const idtype y = index[1];
            const idtype z = index[2];

            hexas->InsertNextCell(8);
            hexas->InsertCellPoint( sx*sy*z + sy*(y+1) + x );
            hexas->InsertCellPoint( sx*sy*z + sy*(y+1) + (x+1) );
            hexas->InsertCellPoint( sx*sy*(z+1) + sy*(y+1) + (x+1) );
            hexas->InsertCellPoint( sx*sy*(z+1) + sy*(y+1) + x );

            hexas->InsertCellPoint( sx*sy*z + sy*y + x );
            hexas->InsertCellPoint( sx*sy*z + sy*y + (x+1) );
            hexas->InsertCellPoint( sx*sy*(z+1) + sy*y + (x+1) );
            hexas->InsertCellPoint( sx*sy*(z+1) + sy*y + x );

            mask[ sx*sy*z + sy*(y+1) + x ] = 1;
            mask[ sx*sy*z + sy*(y+1) + (x+1) ] = 1;
            mask[ sx*sy*(z+1) + sy*(y+1) + (x+1) ] = 1;
            mask[ sx*sy*(z+1) + sy*(y+1) + x ] = 1;
            mask[ sx*sy*z + sy*y + x ] = 1;
            mask[ sx*sy*z + sy*y + (x+1) ] = 1;
            mask[ sx*sy*(z+1) + sy*y + (x+1) ] = 1;
            mask[ sx*sy*(z+1) + sy*y + x ] = 1;
        }
        
        ++it;
    }

    //see which vertices were added and generate coordintes
    for(idtype i=0; i<nvertices_max; i++)
    {
        if(mask>0)
        {
            
        }
    }



    delete[] mask;

    return EXIT_SUCCESS;
}

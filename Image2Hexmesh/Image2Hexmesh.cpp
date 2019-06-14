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
int ReadScalarImage(const char *inputFileName, const char *outputFileName,
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
            return ReadScalarImage<3>(inputFileName, outputFileName, componentType);
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
int ReadScalarImage(const char *inputFileName, const char *outputFileName,
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
        CreateHexMeshFromImage<ImageType>(outputFileName, image);
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
    using PointType = typename TImage::PointType;

    using IdType = int64_t;
    using CoordType = float_t;

    auto spacing = image->GetSpacing();
    auto size = image->GetLargestPossibleRegion().GetSize();

    const IdType sx = size[0]+1;
    const IdType sy = size[1]+1;
    const IdType sz = size[2]+1;


    std::cout<<"Creating elements"<<std::endl;


    //create a mask to know which elements exist
    //alternatively: sx sy (sz-1) + sx (sy-1) + sx-1 + 1, 
    //note that with these sx,sy,sz, the last image coordinate is sx-2, sy-2, sz-2
    const IdType nvertices_max = sx*sy*sz; 
    const IdType ncells_max = (sx-1)*(sy-1)*(sz-1); 

    std::vector<char>  mask( nvertices_max, 0 );

    IteratorType it(image, image->GetLargestPossibleRegion());
    it.GoToBegin();

    std::vector<std::vector<IdType> > hexas;


    vtkSmartPointer<vtkFloatArray> labels =     vtkSmartPointer<vtkFloatArray> ::New();
    labels->SetName( "Label" );
    labels->SetNumberOfComponents( 1 );
    labels->Allocate( ncells_max );

    while (!it.IsAtEnd())
    {
        auto label = it.Get();

        //if (label!=0) 
        {
            auto index = it.GetIndex();


            //the order in which the vetices of the mesh will be generated : cols(x), then rows(y), then slices (z)
            const IdType x = index[0];
            const IdType y = index[1];
            const IdType z = index[2];

            IdType ids[8];
            ids[0] = sx*sy*z + sx*(y+1) + (x+1) ;
            ids[1] = sx*sy*(z+1) + sx*(y+1) + (x+1) ;
            ids[2] = sx*sy*(z+1) + sx*(y+1) + x ;
            ids[3] = sx*sy*z + sx*y + x ;
            ids[4] = sx*sy*z + sx*y + (x+1) ;
            ids[5] = sx*sy*(z+1) + sx*y + (x+1) ;
            ids[6] = sx*sy*(z+1) + sx*y + x ;
            ids[7] = sx*sy*z + sx*(y+1) + x ;

            std::cout<<"Pixel "<<x<<","<<y<<","<<z<<". Ids: ";
            for (int kk=0; kk<8; kk++)
                std::cout<<ids[kk]<<", ";
            std::cout<<std::endl;

            std::cout<<"inserting cell "<<std::endl;

            std::vector<IdType> cell(8);
            for(int kk=0; kk<8; kk++)
            {
                cell[kk] = ids[kk];
                mask[ ids[kk] ] = 1;
            }

            hexas.push_back( cell );

            labels->InsertNextTuple1( label );


            if( ids[3]%50000 ==0 )
                std::cout<<"Progress: "<<ids[3]*100/nvertices_max<<"\r";
        }
        
        ++it;
    }

    std::cout<<std::endl;


    std::cout<<"Mask: ";
    for(int kk=0; kk<nvertices_max; kk++)
    {
        std::cout<<(int) mask[kk]<<", ";
    }
    std::cout<<std::endl;


    //go over the mask and create the points
    //the vertices can be interpreted as pixel centers padded by one row/col/slice at the end and shifted 
    //half voxel 

    std::cout<<"Generating vertex coordinates"<<std::endl;


    vtkSmartPointer<vtkPoints> coord =     vtkSmartPointer<vtkPoints> ::New();
    coord->Allocate(nvertices_max);

    //this permutation is for renumbering the points. To skip nonexisting points. Initialize to -1
    std::vector<IdType> pt_permutation(nvertices_max, -1);

    IdType c=0;
    for(IdType x=0; x<sx; x++)
    {
        for(IdType y=0; y<sy; y++)
        {
            for(IdType z=0; z<sz; z++)
            {
                typename TImage::PointType point;
                typename TImage::IndexType idx;

                idx[0] = x;
                idx[1] = y;
                idx[2] = z;

                const IdType linear_idx = sx*sy*z + sx*y + x;
                if( mask[linear_idx]>0 )
                {
                    image->TransformIndexToPhysicalPoint( idx, point );

                    coord->InsertNextPoint(point[0] - spacing[0]/2.0, point[1] - spacing[1]/2.0, point[2] - spacing[2]/2.0);
                    pt_permutation[linear_idx] = c++;

                    if( linear_idx%50000 ==0 )
                        std::cout<<"Progress: "<<linear_idx*100/nvertices_max<<"\r";
                }
            }
        }
    }

    std::cout<<std::endl;


    std::cout<<"Permutations: ";
    for(int kk=0; kk<nvertices_max; kk++)
    {
        std::cout<<pt_permutation[kk]<<", ";
    }
    std::cout<<std::endl;



    std::cout<<"Renumbering elements and creating mesh"<<std::endl;

    vtkSmartPointer<vtkCellArray> cells =     vtkSmartPointer<vtkCellArray> ::New();

    for(IdType elno = 0; elno<hexas.size(); elno++)
    {
        if( elno%50000 ==0 )
            std::cout<<"Progress: "<<elno*100/hexas.size()<<"\r";

        auto oldids = hexas[elno];
        std::cout<<"Element "<<elno<<" old ids: ";
        for(int kk=0; kk<8; kk++)
            std::cout<<oldids[kk]<<", " ; 
        std::cout<<std::endl;


        assert(oldids->GetNumberOfIds()==8);

        cells->InsertNextCell(8);
        for(int kk=0; kk<8; kk++)
        {
            cells->InsertCellPoint( pt_permutation[ oldids[kk] ] );

            if ( pt_permutation[ oldids[kk] ]<0 )
                std::cout<<"Bad old id: "<<oldids[kk]<<std::endl; 
        }

    }
    std::cout<<std::endl;

    std::cout<<"Saving the mesh "<<outputFileName<<std::endl;
    vtkSmartPointer<vtkUnstructuredGrid> ug =     vtkSmartPointer<vtkUnstructuredGrid> ::New();
    ug->SetPoints(coord);
    ug->SetCells(VTK_HEXAHEDRON,cells);
    ug->GetCellData()->AddArray(labels);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> wr =     vtkSmartPointer<vtkXMLUnstructuredGridWriter> ::New();
    wr->SetFileName(outputFileName);
    wr->SetInputData(ug);
    wr->Write();



    return EXIT_SUCCESS;
}

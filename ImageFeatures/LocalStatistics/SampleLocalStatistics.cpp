//sample voxel statistics in a box around each mesh point
#include <itkImageFileReader.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <stdlib.h>
#include <algorithm>
#include <itkMath.h>
#include <math.h>

#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>

#include <armadillo>

//#define DEBUG_MESSAGES_HOG 


template <typename TImage>
void SampleStatistics(double *pt, double radius, typename TImage::Pointer image, arma::vec& features);


//edrivativeFilter->GetOutput()->SetRequestedRegion(smallRegion);
//derivativeFilter->Update();

int main(int argc, char **argv) {
    
    std::cout<<"Usage: SampleLocalStatistics image.vtk mesh.vtp label_array_name output_features.csv output_label.csv"<<std::endl;
    
//    const char* inputImageFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D-LGE-320.vtk";
//    const char* inputMeshFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D_320_autolabels.vtp";
//    const char* outputHogFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D_320_hog.csv";
//    const char* outputLabelsFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D_320_lbl.csv";
//    const char* labelArrayName = "autolabels";

    
    

    if(argc<6) return EXIT_FAILURE;
    
    int c=1;
    const char* inputImageFileName = argv[c++];
    const char* inputMeshFileName =  argv[c++];
    const char* labelArrayName =  argv[c++];
    const char* outputFeaturesFileName =  argv[c++];
    const char* outputLabelsFileName =  argv[c++];
    //float radius = atof(argv[c++]);
    
    

    
   
    
    //define some region of interest around the point
    arma::vec radii; //does not support more than 1 value for now
    radii << 0.3 << 0.7 <<1.0 <<1.6<<3.5<<5.0; //mm


    //----------------------------------------------------
    //
    //
    //  read the image
    //
    //
    typedef unsigned int PixelType;
    typedef itk::Image< PixelType, 3 > ImageType;
    typedef itk::ImageFileReader< ImageType > ImageReaderType;
    typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
    imageReader->SetFileName(inputImageFileName);

    try {
        imageReader->Update();
    } catch (itk::ExceptionObject& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }


    //
    //
    //----------------------------------------------------
    
    
    //----------------------------------------------------
    //
    //   read the mesh
    //
    //
    //    
    vtkSmartPointer<vtkXMLPolyDataReader> rd = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    rd->SetFileName(inputMeshFileName);
    rd->Update();
    vtkPolyData* mesh = rd->GetOutput();
    
    
    //
    //
    //----------------------------------------------------
        



    //----------------------------------------------------
    //
    //   save the labels
    //
    //
    //    
    std::cout<<"Saving the labels"<<std::endl;
    vtkDataArray* mesh_labels = mesh->GetPointData()->GetArray(labelArrayName);
    arma::uvec labels(mesh_labels->GetNumberOfTuples());
    for(int i=0; i<mesh_labels->GetNumberOfTuples(); i++)
    {
        labels(i) = mesh_labels->GetTuple1(i);
    }
    
    labels.save(outputLabelsFileName, arma::raw_ascii);
    

    //
    //
    //----------------------------------------------------
 
    

    //----------------------------------------------------
    //
    //   for each point of the mesh, take normal and
    //   for each scale (radius+sigma) get the HOG
    //
    //
    //
    
    //create the matrix to store everything
    //every row - HOGS of one point
    
    
    
    std::cout<<"Calculating Features"<<std::endl;
    
    //sample something to know the size of the neighborhood
    double pt[3];
    mesh->GetPoint(0, pt);    
    arma::vec features;
    
    
    const int nfeatures = 4; //4 statistical features
    arma::mat feature_matrix(mesh->GetNumberOfPoints(), nfeatures*radii.size()); 
    //
    
    
 

   
    for(int j=0; j<radii.size(); j++)
    {

        std::cout<<"Radius "<<radii[j]<<std::endl;
    
        for(int i=0; i<mesh->GetNumberOfPoints(); i++)
        {
            if(i%1000==0)
                std::cout<<"Point "<<i<<"/"<<mesh->GetNumberOfPoints()<<"\r"<<std::flush;


        
            double pt[3];
            mesh->GetPoint(i, pt);



            SampleStatistics<ImageType>(pt, radii[j], imageReader->GetOutput(), features);

            if(nfeatures!=features.size())
                std::cout<<std::endl<<"Nonconstant number of samples (vertex "<<i<<"): "<<nfeatures<<" vs "<<features.size()<<std::endl<<std::flush;
            
            //copy the histogram into the matrix
            for(int k=0; k<features.size(); k++)
            {
                feature_matrix(i,k+j*nfeatures) = features(k);
            }
        }
        
        std::cout<<std::endl;
    }
    //
    //
    //----------------------------------------------------
    
    std::cout<<std::endl<<"Saving"<<std::endl;
             
    //save the matrix
    feature_matrix.save(outputFeaturesFileName, arma::arma_binary);
    
    //
    
    return EXIT_SUCCESS;
}



template <typename TImage>
void SampleStatistics(double *pt, double radius, typename TImage::Pointer image, arma::vec& features)
{
    features.set_size(4); //just four features for now
    features.zeros();
    
    
    
    typedef TImage ImageType;

    int radius_pix[3]; //radius in voxels, to be calculated


    radius_pix[0] = std::ceil(radius / image->GetSpacing()[0]);
    radius_pix[1] = std::ceil(radius / image->GetSpacing()[1]);
    radius_pix[2] = std::ceil(radius / image->GetSpacing()[2]);
    
#ifdef DEBUG_MESSAGES_HOG
    std::cout << "Requested point: " << pt[0] << ", " << pt[1] << ", " << pt[2] << std::endl;
    std::cout << "Neighborhood size in mm: " << radius << std::endl;
    std::cout << "Neighborhood size in pix: " << radius_pix[0] << ", " << radius_pix[1] << ", " << radius_pix[2] << std::endl;
#endif
    
    if (radius_pix[0] == 0 || radius_pix[1] == 0 || radius_pix[2] == 0) {
        std::cout << "One of the neighborhood dimensions is zero. Please correct the radius. Aborting." << std::endl;
        return;
    }


    //transform the point from physical s1pace
    typename ImageType::PointType point;
    point[0] = pt[0];
    point[1] = pt[1];
    point[2] = pt[2];

    typename ImageType::IndexType pt_pix;
    image->TransformPhysicalPointToIndex(point, pt_pix);


    //define the region around the point of interest
    typename ImageType::IndexType rgn_idx = {
        {pt_pix[0] - radius_pix[0], pt_pix[1] - radius_pix[1], pt_pix[2] - radius_pix[2]}};
    typename ImageType::SizeType rgn_size = {
        {2 * radius_pix[0] + 1, 2 * radius_pix[1] + 1, 2 * radius_pix[2] + 1}};

    //crop the region so that it is inside
    rgn_idx[0] = std::max(rgn_idx[0], image->GetLargestPossibleRegion().GetIndex(0));
    rgn_idx[1] = std::max(rgn_idx[1], image->GetLargestPossibleRegion().GetIndex(1));
    rgn_idx[2] = std::max(rgn_idx[2], image->GetLargestPossibleRegion().GetIndex(2));

    //set it first as a corner for comparison and then undo that operation
    rgn_size[0] = std::min(rgn_size[0]+rgn_idx[0], image->GetLargestPossibleRegion().GetSize(0))-rgn_idx[0];
    rgn_size[1] = std::min(rgn_size[1]+rgn_idx[1], image->GetLargestPossibleRegion().GetSize(1))-rgn_idx[1];
    rgn_size[2] = std::min(rgn_size[2]+rgn_idx[2], image->GetLargestPossibleRegion().GetSize(2))-rgn_idx[2];

    typename ImageType::RegionType window(rgn_idx, rgn_size);

#ifdef DEBUG_MESSAGES_HOG
    std::cout << "Region: " << window << std::endl;
#endif

    
    itk::ImageRegionConstIterator<ImageType> imageIterator(image, window);
    imageIterator.GoToBegin();

    const int nvoxels = rgn_size[0]*rgn_size[1]*rgn_size[2];
    arma::vec voxels;
    voxels.set_size(nvoxels);
    voxels.zeros();
    
    

    long int i=0;
    while (!imageIterator.IsAtEnd()) {
        // Get the value of the current pixel
        typename ImageType::PixelType val = imageIterator.Get();
    
        //std::cout << val << std::endl;
            

        voxels[i++] = val;
        ++imageIterator;
    }

    //do the statistics

    const double mean = arma::mean(voxels);
    const double std = arma::stddev(voxels);
    
    double E3 = 0; //Accumulate expectation of (X-u)^3
    double E4 = 0; //Accumulate expectation of (X-u)^4
    
    for(int i=0; i<voxels.size(); i++)
    {
        const double diff = voxels[i] - mean;
        const double k = diff*diff*diff/nvoxels;
        E3 += k; //cubed
        E4 += diff*k; //4th power       
    }
       
    features[0] = mean;
    features[1] = std;
    features[2] = E3/(std*std*std); //skewness
    features[3] = E3/(std*std*std*std); //kurtosis
}

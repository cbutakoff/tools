//calculate gradient orientation with respect to the surface normal
#include <itkImageFileReader.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <stdlib.h>
#include <algorithm>
#include <itkMath.h>
#include <math.h>
#include <armadillo>

#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>

//#define DEBUG_MESSAGES_HOG 

template <typename TImage>
void SampleGradHistogram(double *pt, double* normal, double radius, double ngrad_bins, typename TImage::Pointer gradient, arma::uvec& hist);


//edrivativeFilter->GetOutput()->SetRequestedRegion(smallRegion);
//derivativeFilter->Update();

int main(int argc, char **argv) {
    
    std::cout<<"Usage: GetHOG3D image.vtk mesh.vtp label_array_name output_hog.csv output_label.csv"<<std::endl;
    
//    const char* inputImageFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D-LGE-320.vtk";
//    const char* inputMeshFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D_320_autolabels.vtp";
//    const char* outputHogFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D_320_hog.csv";
//    const char* outputLabelsFileName = "/home/costa/Dropbox/code/HOG3D/bin/01D_320_lbl.csv";
//    const char* labelArrayName = "autolabels";

    
    

    if(argc<5) return EXIT_FAILURE;
    
    int c=1;
    const char* inputImageFileName = argv[c++];
    const char* inputMeshFileName =  argv[c++];
    const char* labelArrayName =  argv[c++];
    const char* outputHogFileName =  argv[c++];
    const char* outputLabelsFileName =  argv[c++];
    
    
    

    
   
    const int ngrad_bins = 8; //number of gradient orientation bins
    
    //define some region of interest around the point
    arma::vec radii;
    radii << 2.5 << 5 << 7.5 << 10; //mm
    arma::vec sigmas = radii/3;


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
    
    //calculate normals
    vtkSmartPointer<vtkPolyDataNormals> ng = vtkSmartPointer<vtkPolyDataNormals>::New();
    ng->SetInputData(mesh);
    ng->SplittingOff();
    ng->Update();
    
    vtkDataArray* pnormals = ng->GetOutput()->GetPointData()->GetNormals();
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
    
    std::cout<<"Calculating HOGs"<<std::endl;
    arma::umat hogs(mesh->GetNumberOfPoints(), ngrad_bins*radii.size());
    
    //
    
    
 
    for(int j=0; j<radii.size(); j++)
    {

        std::cout<<"Radius "<<radii[j]<<std::endl;

        typedef itk::GradientRecursiveGaussianImageFilter<ImageType> GradientFilterType;
        typename GradientFilterType::Pointer gradient_filter = GradientFilterType::New();
        gradient_filter->SetInput(imageReader->GetOutput());
        gradient_filter->SetSigma(sigmas[j]);
        gradient_filter->Update();


        GradientFilterType::OutputImageType::Pointer gradient = gradient_filter->GetOutput();           
    
    
        for(int i=0; i<mesh->GetNumberOfPoints(); i++)
         {
            if(i%1000==0)
                std::cout<<"Point "<<i<<"/"<<mesh->GetNumberOfPoints()<<"\r"<<std::flush;
        
            double pt[3];
            mesh->GetPoint(i, pt);


            double n[3];
            pnormals->GetTuple(0, n);

            arma::uvec hist;

            SampleGradHistogram<GradientFilterType::OutputImageType>(pt, n, radii[j], ngrad_bins, gradient, hist);

            //copy the histogram into the matrix
            for(int k=0; k<ngrad_bins; k++)
            {
                hogs(i,k+j*ngrad_bins) = hist(k);
            }
        }
        
        std::cout<<std::endl;
    }
    //
    //
    //----------------------------------------------------
    
    std::cout<<std::endl<<"Saving"<<std::endl;
             
    //save the matrix
    hogs.save(outputHogFileName, arma::csv_ascii);
    
    //
    
    return EXIT_SUCCESS;
}



template <typename TImage>
void SampleGradHistogram(double *pt, double* normal, double radius, double ngrad_bins, typename TImage::Pointer gradient, arma::uvec& hist)
{
    typedef TImage ImageType;

    itk::CovariantVector<double> n;
    n[0] = normal[0];
    n[1] = normal[1];
    n[2] = normal[2];


    const double rad_per_bin = itk::Math::pi/(ngrad_bins-1);

    int radius_pix[3]; //radius in voxels, to be calculated


    radius_pix[0] = std::ceil(radius / gradient->GetSpacing()[0]);
    radius_pix[1] = std::ceil(radius / gradient->GetSpacing()[1]);
    radius_pix[2] = std::ceil(radius / gradient->GetSpacing()[2]);
    
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
    gradient->TransformPhysicalPointToIndex(point, pt_pix);


    //define the region around the point of interest
    typename ImageType::IndexType rgn_idx = {
        {pt_pix[0] - radius_pix[0], pt_pix[1] - radius_pix[1], pt_pix[2] - radius_pix[2]}};
    typename ImageType::SizeType rgn_size = {
        {2 * radius_pix[0] + 1, 2 * radius_pix[1] + 1, 2 * radius_pix[2] + 1}};

    //crop the region so that it is inside
    rgn_idx[0] = std::max(rgn_idx[0], gradient->GetLargestPossibleRegion().GetIndex(0));
    rgn_idx[1] = std::max(rgn_idx[1], gradient->GetLargestPossibleRegion().GetIndex(1));
    rgn_idx[2] = std::max(rgn_idx[2], gradient->GetLargestPossibleRegion().GetIndex(2));

    //set it first as a corner for comparison and then undo that operation
    rgn_size[0] = std::min(rgn_size[0]+rgn_idx[0], gradient->GetLargestPossibleRegion().GetSize(0))-rgn_idx[0];
    rgn_size[1] = std::min(rgn_size[1]+rgn_idx[1], gradient->GetLargestPossibleRegion().GetSize(1))-rgn_idx[1];
    rgn_size[2] = std::min(rgn_size[2]+rgn_idx[2], gradient->GetLargestPossibleRegion().GetSize(2))-rgn_idx[2];

    typename ImageType::RegionType window(rgn_idx, rgn_size);

#ifdef DEBUG_MESSAGES_HOG
    std::cout << "Region: " << window << std::endl;
#endif

    //apply gradient filter to the window
//    typedef itk::GradientRecursiveGaussianImageFilter<ImageType> GradientFilterType;
//    typename GradientFilterType::Pointer gradient_filter = GradientFilterType::New();
//    gradient_filter->SetInput(image);
//    gradient_filter->GetOutput()->SetRequestedRegion(window);
//    gradient_filter->SetSigma(Gauss_sigma);
//    gradient_filter->Update();

    //typedef itk::ImageFileWriter<GradientFilterType::OutputImageType> GradientImageWriterType;
    //GradientImageWriterType::Pointer wr = GradientImageWriterType::New();
    //wr->SetInput(gradient_filter->GetOutput());
    //wr->SetFileName("gradient.vtk");
    //wr->Update();

    hist.set_size(ngrad_bins);
    hist.zeros();
    
//    itk::ImageRegionConstIterator<typename GradientFilterType::OutputImageType> imageIterator(gradient_filter->GetOutput(), window);
    itk::ImageRegionConstIterator<ImageType> imageIterator(gradient, window);
    imageIterator.GoToBegin();
    while (!imageIterator.IsAtEnd()) {
        // Get the value of the current pixel
        typename ImageType::PixelType val = imageIterator.Get();
    
        //std::cout << val << std::endl;
            
        const double angle_cos = val*n/val.GetNorm();
        const double angle = std::acos(angle_cos);
        const int nbin = round(angle/rad_per_bin); 
        //std::cout << "Angle with the normal " << angle*180/itk::Math::pi <<" deg, bin: "<<nbin<<std::endl;
        hist[nbin]++;    
        
        ++imageIterator;
    }

    arma::vec cntr(ngrad_bins);
    for(int i=0; i<ngrad_bins; i++) cntr[i]=i*rad_per_bin;
    
    
#ifdef DEBUG_MESSAGES_HOG
    std::cout << "Center angles:" << cntr*180/itk::Math::pi <<std::endl;
    std::cout << "Histogram    :" << hist <<std::endl;    
#endif
    
}

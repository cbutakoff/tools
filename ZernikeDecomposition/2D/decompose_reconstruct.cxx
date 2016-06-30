// COMPUTE THE POLYNOMIALS Vmn and save it as *.mat to visualize the Zernike polynomials in MATLAB
// DECOMPOSE IMAGE using the computed polynomials
// RECOMPOSE the image
// USING BOOST

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <complex>
#include <vector>
#include <iostream>
#include <string>
#include <armadillo>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/complex.hpp>

#include "../defs.h"

// GLOBAL TYPE DEFINITIONS

//typedef boost::multiprecision::cpp_dec_float_50 mp_float; // slower



// Generate the vectors of M and N 
void generateMN(std::vector<int>& MVect, std::vector<int>& NVect, int order);



//__________________MAIN__________________

int main(int argc, char *argv[]) {

    // Check command line arguments
    if (argc < 3) {
        std::cerr << "Usage: " << std::endl;
        std::cerr << "---> " << argv[0] << " input_image.format" << " order" << std::endl;

        return EXIT_FAILURE;
    }

    // Reading input image
    const unsigned int Dimension = 2;
    typedef unsigned short PixelType;

    typedef itk::Image <PixelType, Dimension> ImageType;
    typedef itk::ImageFileReader <ImageType> ReaderType;

    ReaderType::Pointer reader = ReaderType::New();
    const char *filename = argv[1];
    reader->SetFileName(filename);
    reader->Update();


    // Image Iterator definitions
    ImageType::Pointer image = reader->GetOutput();
    typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ConstIteratorType;
    ImageType::RegionType region = image->GetLargestPossibleRegion();
    ConstIteratorType imageIt(image, region);


    // Take Image Height and Width to define Vmn
    ImageType::SizeType size = region.GetSize();

    std::cout << "\nInput image size: " << size << std::endl;
    int imageHeight = region.GetSize()[0];
    int imageWidth = region.GetSize()[1];


    // Generate The M and N array used for Computing the polynomials, given order
    int order = atoi(argv[2]); 
    std::vector<int> MVect; 
    std::vector<int> NVect;

    generateMN(MVect, NVect, order); // Generate vectors N and M  
    int Nsize = NVect.size(); 
    int Msize = MVect.size(); 
    
     
    // Zmn matrix definition
    VectorType zVect(Nsize); // ** Each slide contains a Zmn  
    
    // Decomposition
    char filename1[1000];
    MatrixType vMatrix(imageHeight, imageWidth);
    for (int i = 0; i < Nsize; ++i) {
        // imageIt is defined in the start of main
        imageIt.GoToBegin();
        int height_cont = 0;
        int width_cont = 0;
        std::complex<mp_float> sumZ(0,0);
        
        //load the matrix
        sprintf(filename1, "poly/poly_o%05d_r%05d",MVect[i], NVect[i]);
      
	{
	  std::ifstream file1(filename1, std::ios::binary);
        
            boost::archive::binary_iarchive ia(file1);
            ia & vMatrix;
        }        
  
                
        
        while(!imageIt.IsAtEnd()){ 
                if(width_cont == imageWidth){
                    width_cont = 0;
                    ++height_cont;
                }
		
                mp_float tempImageData = (imageIt.Get()); // Read image pixel 
               // mp_float tempMatrixData = vMatrix[i];
                sumZ += std::conj(vMatrix(height_cont, width_cont)) * tempImageData;
                
                //std::cout << "X: " << width_cont << ", Y: " << height_cont << ", Value: " << imageIt.Get() <<"\r"<< std::flush;
                ++imageIt;
                ++width_cont;
        }
        
        std::cout << "Decomposing image: " << i + 1 << " / " << Nsize << "\r" << std::flush;        
        //zVect[i] = std::complex<mp_float>((mp_float(MVect[i]+1)/( boost::math::constants::pi<mp_float>())),0) * sumZ; // (m+1)/pi and Complex dot product
        zVect[i] = std::complex<mp_float>((mp_float(MVect[i]+1)/((imageWidth - 1)*(imageWidth - 1))),0) * sumZ; // (m+1)/pi and Complex dot product

    }
    std::cout << std::endl; // Aesthetic purpose


    // <><><><><><><><><>  IMAGE RECONSTRUCTION  <><><><><><><><><>

    MatrixType reconstrMatrix(imageHeight, imageWidth); 

    int leftOut = 0;
    for (int i = 0; i < Nsize-leftOut; ++i) {  // Starting from 4th polynomial ?--> Search Why
        //load the matrix
        sprintf(filename1, "poly/poly_o%05d_r%05d",MVect[i], NVect[i]);
        std::ifstream file1(filename1, std::ios::binary);
        {
            boost::archive::binary_iarchive ia(file1);
            ia & vMatrix;
        }        
        
    
        reconstrMatrix += (zVect[i + leftOut]) * (vMatrix); 
        std::cout << "Reconstructing image: " << i + 1 << " / " << Nsize - leftOut << "\r" << std::flush;
    }
    std::cout << std::endl; 


    // Convert to Armadillo type and double
    arma::Mat<double> ArmaReconstrMatrix(imageHeight, imageWidth);
    for (int i = 0; i < imageHeight; ++i){
        for (int j = 0; j < imageWidth; ++j){
            ArmaReconstrMatrix(i,j) = (real(reconstrMatrix(i,j)).convert_to<double>());
//            std::cout <<"matr("<<i<< ", "<<j<<"): "<<ArmaReconstrMatrix(i,j)<<std::endl;        
        }
    }
    

    char outname[100];
    sprintf(outname,"reconstr_%d.mat",order);

    // SAVE RECONSTRUCTION --> How ?
    std::cout << "Saving " << " [Reconstructed Image Size = " << ArmaReconstrMatrix.n_cols << " x " << ArmaReconstrMatrix.n_rows  << "] ... ..." << std::endl;
    ArmaReconstrMatrix.save(outname, arma::raw_binary); 

    return EXIT_SUCCESS;
}





//Generation of the vector of M and N

void generateMN(std::vector<int>& MVect, std::vector<int>& NVect, int order) {
    //std::cout << M.size() << std::endl;

    for (int m = 0; m <= order; ++m) {
        for (int n = -m; n <= m; ++n) {
            if (fmodf(m - abs(n), 2) == 0) {
                MVect.push_back(m);
                NVect.push_back(n);
                //std::cout << MVect[i] << " "<< NVect[i] << std::endl;
            }
        }
    }
    //    MVect.resize(i); // To have the actual "used" size
    //    NVect.resize(i);
}




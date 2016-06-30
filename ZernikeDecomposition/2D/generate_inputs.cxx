// COMPUTE THE POLYNOMIALS Vmn and save it as *.mat to visualize the Zernike polynomials in MATLAB
// DECOMPOSE IMAGE using the computed polynomials
// RECOMPOSE the image
// USING BOOST


#include <complex>
#include <vector>
#include <iostream>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

// GLOBAL TYPE DEFINITIONS



//__________________PROTOTYES__________________

// Generate the vectors of M and N 
void generateMN(std::vector<int>& MVect, std::vector<int>& NVect, int order);


//__________________MAIN__________________

int main(int argc, char *argv[]) {

    // Check command line arguments
    if (argc < 3) {
        std::cerr << "Usage: " << std::endl;
        std::cerr << "---> " << argv[0] << " height width max_order" << std::endl;

        return EXIT_FAILURE;
    }

    // Reading input image

    int imageHeight = atoi(argv[1]);
    int imageWidth = atoi(argv[2]);
    int maxorder = atoi(argv[3]);
    std::cout << "Height: " << imageHeight << std::endl;
    std::cout << "imageWidth: " << imageWidth << std::endl;    
    std::cout << "order: " << maxorder << std::endl;

    
    
    // Generate The M and N array used for Computing the polynomials, given order
    std::vector<int> MVect; 
    std::vector<int> NVect;

    generateMN(MVect, NVect, maxorder); // Generate vectors N and M  
    int Nsize = NVect.size(); 

    system("mkdir input");
    for (int i = 0; i < Nsize; ++i) {
        char filename[100];
        sprintf(filename,"input/in-%d.log",i+1);
        std::ofstream f(filename);
        f<< imageHeight << std::endl;
        f<< imageWidth << std::endl;
        f<< MVect[i] << std::endl;
        f<< NVect[i] << std::endl;
    } 
    
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




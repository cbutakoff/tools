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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/format.hpp>

#include "defs.h"

// GLOBAL TYPE DEFINITIONS



//__________________PROTOTYES__________________

// Generate the vectors of M and N 
void generateMN(std::vector<int>& MVect, std::vector<int>& NVect, int order);

// Computes Rmn using gmp for large numbers
mp_float CalcRmn(int m, int n, mp_float r);

// Computes Vmn. It uses CalcRmn.
void ComputePolynomial(int M, int N, double rc, double cc, double scale, MatrixType &vMatrix, int nrows, int ncols);
//void ComputePolynomial(int M, int N, double rc, double cc, double scale, int nrows, int ncols);

typedef std::vector<mp_int> FactorialArrayType;
void CalcFactorials(FactorialArrayType& F, int n);

//__________________MAIN__________________

int main(int argc, char *argv[]) {

    // Check command line arguments
    if (argc < 2) {
        std::cerr << "Usage: " << std::endl;
        std::cerr << "---> " << argv[0] << " input_file" << std::endl;

        return EXIT_FAILURE;
    }

    // Reading input file
    // Has height, width, order, repeat
    int imageHeight;
    int imageWidth;
    int order;
    int repeat;

    char *infile = argv[1];
    std::cout<<"reading "<<infile;
    std::string str;
    std::ifstream ff(infile);
    getline (ff,str);
    imageHeight = atoi(str.c_str());
    getline (ff,str);
    imageWidth = atoi(str.c_str());
    getline (ff,str);
    order = atoi(str.c_str());
    getline (ff,str);
    repeat = atoi(str.c_str());
    
        
    std::cout << "Height: " << imageHeight << std::endl;
    std::cout << "imageWidth: " << imageWidth << std::endl;
    std::cout << "order: " << order << std::endl;
    std::cout << "repeat: " << repeat << std::endl;
    
    
    // Vmn matrix definition
    MatrixType vMatrix(imageHeight,imageWidth);

    // Parameters for scaling definition
    double rc = ((imageHeight + 1) / 2);
    double cc = ((imageWidth + 1) / 2);
    double scale = (imageHeight / 2);

  
    // COMPUTING THE POLYNOMIALS [

    ComputePolynomial(order, repeat, rc, cc, scale, vMatrix, imageHeight, imageWidth);
   

    char filename1[1000];
    system("mkdir poly");
    sprintf(filename1, "poly/poly_o%05d_r%05d",order, repeat);
    std::ofstream file(filename1, std::ios::binary);
    {
        boost::archive::binary_oarchive oa(file);
        oa & vMatrix;
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




//factorial array must contain 0!, 1!, ..., n! (n-first argument, order)
mp_float CalcRmn(int n, int m, mp_float r, FactorialArrayType& F) {
    //n - order
    //m - repetition

    if((m-n)%2!=0) return 0;
    
    //create storage for the factorials
//    mpz_t *F = new mpz_t[n+1];
/*    
    std::vector<mp_int> F(n+1);
    
    //precalculate the factorials including 0! and n!
    F[0] = 1; //0!
    for(int i=1; i<=n; i++)
    {
        F[i] = F[i-1]*i;
    }
*/    
    
    mp_int k;
    mp_int k2;
    
    const int nm2p = (n+abs(m))/2;
    const int nm2m = (n-abs(m))/2;
    
    mp_float ff, rpow, result;
    result = 0;
    
    
    for(int l=0; l <= nm2m; l++)
    {   

        //numerator
        k  = (l%2?-1:1)*F[n-l];

        //denominator
        k2 = F[l]*F[nm2p-l]*F[nm2m-l];
        
        //fraction
        //calc factorial fraction (-1)^l (n-l)!/[ l! ((n+m)/2-s)! ((n-m)/2-s)!]
        k /= k2;
                                        
        ff.assign(k); //convert the factorial fraction to floating point
        rpow = pow(r,n-2*l); //calculate r^(n-2l)
        ff = ff*rpow; //calculate (-1)^l (n-l)!/[ l! ((n+m)/2-s)! ((n-m)/2-s)!] *r^(n-2l)
        
        //accumulate the results
        result += ff;
    }

    //gmp_printf("GMP Rmn = %Ff\n", result);

    return result;
}


void ComputePolynomial(int M, int N, double rc, double cc, double scale, MatrixType &vMatrix, int nrows, int ncols){
//void ComputePolynomial(int M, int N, double rc, double cc, double scale, int nrows, int ncols) {
    // As it was done in MATLAB

    FactorialArrayType F;
    CalcFactorials(F, M);
    
    for (int r = 0; r < nrows; ++r) {
        //std::cout<<"Row "<<r<<"/"<<nrows<<"\r"<<std::flush;
        mp_float rn = (r - rc) / scale;

        for (int c = 0; c < ncols; ++c) {
            mp_float cn = (c - cc) / scale;
            mp_float rho = sqrt(rn * rn + cn * cn);
            mp_float theta = atan2(rn, cn);

            if (theta < 0) {
                theta = theta + 2 * boost::math::constants::pi<mp_float>();
            }
            if (rho > 1) {
                continue;
            }

            mp_float Rmn = CalcRmn(M, N, rho, F);
            //vMatrix(r,c) = Rmn * exp((-1)) * N * theta;
            std::complex<mp_float> Vmn = std::complex<mp_float>(Rmn,0) * std::complex<mp_float>(cos(N * theta), sin(N * theta));
            vMatrix(r,c) = Vmn;
            mp_float abs_val = abs(Vmn);
            if(abs_val>1)
            {
                std::cout<<std::endl<<boost::format("|Vmn|: %.f") % abs_val<<std::endl;
            }
            // std::cout<<"Rmn: "<<Rmn<<std::endl;
            // std::cout<<"Rho: "<<rho<<std::endl;
            // std::cout<<"Theta: "<<theta<<std::endl;
        }
    }
    //return V;
}

//calculate factorials 0!, 1!, ..., n!
void CalcFactorials(FactorialArrayType& F, int n)
{
    F.resize(n+1);
    
    //precalculate the factorials including 0! and n!
    F[0] = 1; //0!
    for(int i=1; i<=n; i++)
    {
        F[i] = F[i-1]*i;
    }
}

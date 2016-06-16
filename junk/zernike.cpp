/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

//#include <math.h>
//#include <cmath>


#include<armadillo>


#include<vtkSmartPointer.h>
#include<vtkPolyData.h>
#include<vtkCell.h>
#include<vtkIdList.h>
#include<vtkTriangle.h>
#include<vtkPolyDataReader.h>
#include<vtkPolyDataWriter.h>
#include<vtkFloatArray.h>
#include<vtkPointData.h>

#include <math.h>
#include <complex>

#include "imported/pgm2itkvol/mcutil.h"

class Zernike
{
public: 
    typedef arma::Cube<long int> LICube;
    typedef std::complex<double> Complex;
    typedef arma::Cube< Complex > CCube;
    
    void SetMesh(vtkPolyData* mesh) {m_Surface = mesh; };
    vtkPolyData* GetMesh() {return m_Surface; };

    int GetOrder() {return m_Order; };
    void SetOrder(int order) {m_Order = order;};
    
    void SaveGeometricMoments(const char* filename, bool binary=true) 
        { m_G.save(filename, binary?arma::arma_binary:arma::arma_ascii); };
    void LoadGeometricMoments(const char* filename) { m_G.load(filename); m_Order=m_G.n_cols-1; };
    double GetGeometricMoment(int i, int j, int k) {return m_G(i,j,k);};


    void SaveZernikeMoments(const char* filename, bool binary=true) 
        { m_Z.save(filename, binary?arma::arma_binary:arma::arma_ascii); };
    void LoadZernikeMoments(const char* filename) { m_Z.load(filename); m_Order=m_Z.n_cols-1; };
    Complex GetZernikeMoment(int i, int j, int k) {return m_Z(i,j,k);};

    
    void RunTests();
    

    void CalculateSurfaceGeometricMoments();
    void CalculateZernikeMoments();
    
    vtkPolyData* Reconstruct(double x_min, double x_max, 
                            double y_min, double y_max,
                            double z_min, double z_max, double step) const;

    Zernike():m_Order(0) {};

protected:

    void ComputeTrinomialFactors(bool disable_order_check=false);
    long int BinomialFactor(int, int) const;
    long int TrinomialFactor(int, int, int) const;
    long int Factorial(int) const;
    
private:
    
    LICube m_TrinomialFactors;
    int m_Order; //requested order
    arma::cube m_G; //geometric moments
    CCube m_Z; //Zernike moments
    
    vtkPolyData* m_Surface;
};


long int Zernike::TrinomialFactor(int a, int b, int c) const
{
    return Factorial(a+b+c)/(Factorial(a)*Factorial(b)*Factorial(c));
}


long int Zernike::BinomialFactor(int p, int a) const
{
    return Factorial(p)/(Factorial(a)*Factorial(p-a));
}



void Zernike::RunTests()
{
    //Test the Factorial calculation
    std::cout<<"Testing factorial calculation: ";
    bool result= Factorial(0) == 1 && Factorial(3)==6;
    if(result)
        std::cout<<"passed"<<std::endl;
    else
        std::cout<<"failed"<<std::endl;
    
    //Test the trinomial calculation
    std::cout<<"Testing the trinomial calculation: ";
    m_Order = 2;
    ComputeTrinomialFactors(true);
    LICube::iterator a = m_TrinomialFactors.begin();
    LICube::iterator b = m_TrinomialFactors.end();

    //m_TrinomialFactors.save("trinomial.dat",arma::arma_ascii);


    long int correct_cube[] = { 1, 1, 1, 1, 2, 3, 1, 3, 6, 1, 2, 3, 2, 6, 12, 3,
        12, 30, 1, 3, 6, 3, 12, 30, 6, 30, 90};
    
    bool passed = true;
    int j=0;
    for(LICube::iterator i=a; i!=b; ++i)
    {
        passed = (*i) == correct_cube[j++];
    }
    if(passed)
        std::cout<<"passed"<<std::endl;
    else
        std::cout<<"failed"<<std::endl;


    std::cout<<"Testing binomial calculation: ";
    result= BinomialFactor(0,0) == 1 && BinomialFactor(5,3)==10;
    if(result)
        std::cout<<"passed"<<std::endl;
    else
        std::cout<<"failed"<<std::endl;
    

    std::cout<<"Testing trinomial calculation: ";
    result= TrinomialFactor(2,3,4) == 1260;
    if(result)
        std::cout<<"passed"<<std::endl;
    else
        std::cout<<"failed"<<std::endl;

    
}

long int Zernike::Factorial(int x) const
{
    long int result = x;
    
    if(x!=0)
        for(int i=x-1; i>0; i--)
        {
            result *= i; 
        }
    else
        result = 1;
    
    return result;
}



void Zernike::ComputeTrinomialFactors(bool disable_order_check)
{
    const int N = m_Order+1;
    m_TrinomialFactors.set_size(N,N,N);
    m_TrinomialFactors.fill(-1);
    
    
    for(int i=0; i<N; i++) 
    {
        const long int i_factor = Factorial(i);
        for(int j=0; j<N; j++)
        {
            const long int ij_factor = Factorial(j)*i_factor;

            for(int k=0; k<N; k++)
            {    
                //std::cout<<int(!disable_order_check)<<std::endl;
                //std::cout<<int(i+j+k<=m_Order)<<std::endl;
                //std::cout<<int((!disable_order_check) && (i+j+k<=m_Order))<<std::endl;
                if( disable_order_check || (i+j+k<=m_Order) )
                {
                    //std::cout<<"in";
                    const long int ijk_factor = Factorial(k)*ij_factor;
                    m_TrinomialFactors(i,j,k) = Factorial(i+j+k)/ijk_factor;
                    //std::cout<<"m_TrinomialFactors("<<i<<","<<j<<","<<k<<")="<<m_TrinomialFactors(i,j,k)<<std::endl;
                }
                
            }

        }

    }
}




void Zernike::CalculateSurfaceGeometricMoments()
{
    ComputeTrinomialFactors();

    const int N = m_Order+1;
    m_G.set_size(N, N, N);
    m_G.fill(0);

//    double total_area = 0;
    for(int i=0; i<m_Surface->GetNumberOfCells(); i++)
    {
        //std::cout<<"Processing triangle "<<i<<"/"<<m_Surface->GetNumberOfCells()<<std::endl;
        vtkCell* cell = m_Surface->GetCell(i);
        vtkIdList* ids = cell->GetPointIds();
        double p1[3];
        double p2[3];
        double p3[3];
        
        assert(ids->GetNumberOfIds()==3);
        m_Surface->GetPoint(ids->GetId(0),p1);
        m_Surface->GetPoint(ids->GetId(1),p2);
        m_Surface->GetPoint(ids->GetId(2),p3);
        
        //calculating C^(a)_{ijk}
        arma::cube C1(N, N, N);
        arma::cube C2(N, N, N);
        arma::cube C3(N, N, N);
        C1.fill(-1);
        C2.fill(-1);
        C3.fill(-1);

        for(int i=0; i<N; i++)
        {
            const double pp10 = std::pow(p1[0],i);
            const double pp20 = std::pow(p2[0],i);
            const double pp30 = std::pow(p3[0],i);
            
            for(int j=0; j<N; j++)
            {

                const double pp11 = std::pow(p1[1],j);
                const double pp21 = std::pow(p2[1],j);
                const double pp31 = std::pow(p3[1],j);

                for(int k=0; k<N; k++)
                {
                    if(i+j+k<=m_Order)
                    {
                        const double pp12 = std::pow(p1[2],k);
                        const double pp22 = std::pow(p2[2],k);
                        const double pp32 = std::pow(p3[2],k);

                        assert(m_TrinomialFactors(i,j,k)!=-1);
                        
                        C1(i,j,k) = m_TrinomialFactors(i,j,k)*pp10*pp11*pp12;
                        C2(i,j,k) = m_TrinomialFactors(i,j,k)*pp20*pp21*pp22;
                        C3(i,j,k) = m_TrinomialFactors(i,j,k)*pp30*pp31*pp32;
                    }
                    
                }
            }
        }
        
        //calculate D
        arma::cube D(N, N, N);
        D.fill(-1);

        for(int a=0; a<N; a++)
            for(int b=0; b<N; b++)            
                for(int c=0; c<N; c++)
                    if(a+b+c<=m_Order)
                    {
                        D(a,b,c)=0;
                        for(int i2=0; i2<=a; i2++)
                            for(int j2=0; j2<=b; j2++)
                                for(int k2=0; k2<=c; k2++)
                                {
                                    assert(C2(i2,j2,k2)!=-1);
                                    assert(C3(a-i2,b-j2,c-k2)!=-1);
                                    D(a,b,c) += C2(i2,j2,k2)*C3(a-i2,b-j2,c-k2);
                                }
                    }
        
        
        //Calculate Sijk and increment the geometric moment
        const double area = dynamic_cast<vtkTriangle*>(cell)->ComputeArea();

        arma::cube S(N, N, N);
        S.fill(-1);
        
        for(int i=0; i<N; i++)
            for(int j=0; j<N; j++)            
                for(int k=0; k<N; k++)
                    if(i+j+k<=m_Order)
                    {
                        S(i,j,k)=0;
                        for(int i1=0; i1<=i; i1++)
                            for(int j1=0; j1<=j; j1++)
                                for(int k1=0; k1<=k; k1++)
                                {
                                    assert(C1(i1,j1,k1)!=-1);
                                    assert(D(i-i1,j-j1,k-k1)!=-1);
                                    S(i,j,k) += C1(i1,j1,k1)*D(i-i1,j-j1,k-k1);
                                }
                        S(i,j,k) = S(i,j,k)*Factorial(i)*Factorial(j)*Factorial(k)/Factorial(m_Order+2);

                        m_G(i,j,k) += S(i,j,k)*2.0*area;                        
                    }
//        total_area += area;
    }
//    std::cout<<"total area: "<<total_area<<std::endl;
//    std::cout<<"m000: "<<m_G(0,0,0)<<std::endl;
}


void Zernike::CalculateZernikeMoments()
{
    const int N = m_Order+1;
    CCube V(N,N,N);
    CCube W(N,N,N);
    CCube X(N,N,N);
    CCube Yhat(N,N,N);
    V.fill(-1);
    W.fill(-1);
    X.fill(-1);
    Yhat.fill(-1);


    Complex im = 1.0i;
    
    //std::cout<<"Calculating V"<<std::endl;
    //Calculate V
    for(int a=0; a<=m_Order; a++)
        for(int b=0; b<=m_Order; b++)
            for(int c=0; c<=m_Order; c++)
                if(2*a+b+c<=m_Order)
                {
                    V(a,b,c)=0;
                    for(int alpha=0; alpha<=a+c; alpha++)
                    {
                        V(a,b,c) += BinomialFactor(a+c,alpha)*
                                m_G(2*a+c-alpha,alpha,b)*
                                std::pow(im,alpha);
                    }
                }
    //V.save("v.dat",arma::arma_ascii);
    
    //std::cout<<"Calculating W"<<std::endl;
    //Calculate W
    for(int a=0; a<=m_Order; a++)
        for(int b=0; b<=m_Order; b++)
            for(int c=0; c<=m_Order; c++)
                if(2*a+b+c<=m_Order)
                {
                    W(a,b,c)=0;
                    for(int alpha=0; alpha<=a; alpha++)
                    {
                        //assert(V(a-alpha, b, c+2*alpha)!=-1);
                        
                        W(a,b,c) += std::pow(-1.0,alpha)*
                                    std::pow(2,a-alpha)*
                                    BinomialFactor(a, alpha)*
                                    V(a-alpha, b, c+2*alpha);
                    }
                }

    //W.save("w.dat",arma::arma_ascii);
    
    //std::cout<<"Calculating X"<<std::endl;
    //Calculate X
    for(int a=0; a<=m_Order; a++)
        for(int b=0; b<=m_Order; b++)
            for(int c=0; c<=m_Order; c++)
                if(2*a+b+c<=m_Order)
                {
                    X(a,b,c)=0;
                    for(int alpha=0; alpha<=a; alpha++)
                    {
                        //assert( W(a-alpha, b+2*alpha, c)!=-1 );
                        Complex t1 = BinomialFactor(a, alpha); 
                        X(a,b,c) +=  t1*W(a-alpha, b+2*alpha, c);
                    }
                }
    

    //X.save("X.dat",arma::arma_ascii);
    
    //std::cout<<"Calculating Y"<<std::endl;
    //calculate Y
    arma::cube Y(N,N,N);
    Y.fill(-1);
    for(int l=0; l<=m_Order; l++)
    {
//        if( (m_Order-l)%2==0 ) //continue only if even
//        {
            const double pow2l = std::pow(2,l);

            for(int j=0; j<=(m_Order-l)/2; j++)
            {
                const double pow1j = std::pow(-1,j); 
                
                for(int m=0; m<=m_Order; m++)
                {
                    Y(m,l,j) = pow1j*
                            (std::sqrt(2.0*l+1.0)/pow2l)*
                            TrinomialFactor( m, j, l-m-2*j )*
                            BinomialFactor( 2*(l-j), l-j )/
                            std::sqrt( double(TrinomialFactor( m, m, l-m )) );                        
                }
            }
//        }
    }

    //Y.save("Y.dat",arma::arma_ascii);

    
    //std::cout<<"Calculating Yhat"<<std::endl;
    //calcuate Yhat
    for(int l=0; l<=m_Order; l++)
//        if( (m_Order-l)%2==0 ) //continue only if even
//        {
            for(int v=0; v<=(m_Order-l)/2; v++)
            {
                for(int m=0; m<=l; m++)
                {
//                    if( (l-m)%2==0 ) //continue if even
//                    {
                        Yhat(m,l,v) = 0;

                        for(int j=0; j<=(l-m)/2; j++)
                        {
                            assert(Y(m,l,j)!=-1);
                            //assert(X(v+j, l-m-2*j, m)!=-1);
                            Yhat(m,l,v) = Y(m,l,j)*X(v+j, l-m-2*j, m);
//                        std::cout<<"mlv="<<m<<" "<<l<<" "<<v<<std::endl;
//                        std::cout<<"Yhat(m,l,v)="<<Yhat(m,l,v)<<std::endl;
                        }
//                    }
                }
//            }
        }

    //Yhat.save("Yhat.dat",arma::arma_ascii);
    
    
    //std::cout<<"Calculating Q"<<std::endl;
    //calculate Q
    arma::cube Q(N,N,N);
    Q.fill(-1);
    for(int l=0; l<=m_Order; l++)
    {
//        if((m_Order-l)%2==0)
//        {
            for(int n=0; n<=m_Order; n++)
                if( (n-l)%2==0 )
                {
                    const int k=(n-l)/2;
                    for(int v=0; v<=(m_Order-l)/2; v++)
                    {
                        Q(k,l,v) = std::pow(-1, k+v)/std::pow(4,k) * 
                                std::sqrt( (2.0*l+4.0*k+3.0)/3.0 ) *
                                TrinomialFactor(v, k-v, l+v+1)*
                                BinomialFactor(2*(l+v+1+k), l+v+1+k)/
                                BinomialFactor(2*(l+v+1), l+v+1);
//                        std::cout<<"klv="<<k<<" "<<l<<" "<<v<<std::endl;
//                        std::cout<<"Q(k,l,v)="<<Q(k,l,v)<<std::endl;
                        
                    }
                }
//        }
    }
    
    //Q.save("Q.dat",arma::arma_ascii);
    
    //std::cout<<"Calculating Z"<<std::endl;    
    m_Z.set_size(N,N,N);
    m_Z.fill(-1);
    
    
    //calculate Z
    for(int n=0; n<=m_Order; n++)
        for(int l=0; l<=n; l++)
            for(int m=0; m<=l; m++)
            {
                if( (n-l)%2==0 )
                {
//                    std::cout<<n<<" "<<l<<" "<<m<<" "<<std::endl;
                    const int k=(n-l)/2;

                    m_Z(m,n,l)=0;
                    for(int v=0; v<=(n-l)/2; v++)
                    {
                        //assert(Yhat(m,l,v)!=-1);
                        assert(Q(k,l,v)!=-1);
                        m_Z(m,n,l) += std::conj(Yhat(m,l,v))*Q(k,l,v);
//                        std::cout<<"Yhat(m,l,v)="<<Yhat(m,l,v)<<std::endl;
//                        std::cout<<"m_Z(m,n,l)="<<m_Z(m,n,l)<<std::endl;
//k                        std::cout<<"klv="<<k<<" "<<l<<" "<<v<<std::endl;

                    }
                    m_Z(m,n,l) = 3.0/(4.0*M_PI) * m_Z(m,n,l);
//                  std::cout<<"mnl="<<m<<" "<<n<<" "<<l<<std::endl;
//                  std::cout<<"m_Z(m,n,l)="<<m_Z(m,n,l)<<std::endl;
                }
            }

}




vtkPolyData* Zernike::Reconstruct(double x_min, double x_max, 
                            double y_min, double y_max,
                            double z_min, double z_max, double step) const
{
    vtkPolyData* result = vtkPolyData::New();
    
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> f = vtkSmartPointer<vtkFloatArray>::New();
    f->SetName("field");
    f->SetNumberOfComponents(1);
    
    for(double x=x_min; x<=x_max; x+=step)
        for(double y=y_min; y<=y_max; y+=step)
            for(double z=z_min; z<=z_max; z+=step)
            {
                pts->InsertNextPoint(x,y,z);
                
                double indic = 0;

               
                for(int n=0; n<=m_Order; n++)
                {
                    for(int l=0; l<=n; l++)
                        if((n-l)%2==0)
                        {
                            for(int m=-l; m<=l; m++)
                            {
                                std::complex<double> Zmnl;
                                
                                if(m>=0)
                                    Zmnl = m_Z(m,n,l);
                                else
                                    Zmnl = std::conj(m_Z(-m,n,l))*std::pow(-1,-m);
                                
                                const int k = (n-l)/2;
                                
                                const int m1 = m;
                                m = std::abs(m);
                                
                                std::complex<double> pZmnl = 0.0+0.0i;                                
                
                                for(int v=0; v<=k; v++)
                                {
                                    const double Qklv =  std::pow(-1, k+v)/std::pow(4,k) * 
                                            std::sqrt( (2.0*l+4.0*k+3.0)/3.0 ) *
                                            TrinomialFactor(v, k-v, l+v+1)*
                                            BinomialFactor(2*(l+v+1+k), l+v+1+k)/
                                            BinomialFactor(2*(l+v+1), l+v+1);

                                    for(int j=0; j<=(l-m)/2; j++)
                                    {
                                        const double pow2l = std::pow(2,l);
                                        const double pow1j = std::pow(-1,j); 
                                        const double Ymlj = pow1j*
                                                (std::sqrt(2.0*l+1.0)/pow2l)*
                                                TrinomialFactor( m, j, l-m-2*j )*
                                                BinomialFactor( 2*(l-j), l-j )/
                                                std::sqrt( double(TrinomialFactor( m, m, l-m )) );  

                                        std::complex<double> pp = x+y*1i;

                                        pZmnl += std::pow( pp, m )*Qklv * Ymlj * 
                                                std::pow( std::pow(x,2)+std::pow(y,2)+std::pow(z,2), v+j )*
                                                std::pow( z, l-m-2*j );
                                    }
                                }
                                
                                if(m1<0)
                                    pZmnl = std::conj(pZmnl)*std::pow(-1,m);

                                indic += std::abs(pZmnl*Zmnl);
                            }
                        }                    
                }
                
                f->InsertNextTuple1(indic);
            }
        
    result->SetPoints(pts);
    result->GetPointData()->AddArray(f);
    return result;
}














bool TestGeometricMoments(const char* meshname);
bool TestZernikeMoments(const char* meshname);


int main(int argc, char** argv)
{
    Zernike z;
    
    const char* meshname = argv[1];
    const int order = atoi(argv[2]);
    
    std::cout<<"Version 0.99"<<std::endl;
    std::cout<<"Filename: "<<meshname<<std::endl;
    std::cout<<"Order: "<<order<<std::endl;
    
    
    
    if (!TestZernikeMoments(meshname))
    {
        std::cout<<"Error in geometric moments code"<<std::endl;
        return -1;
    }
  
    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(meshname);
    rdr->Update();
    
    //center the mesh, just in case
    double cntr[] = {0,0,0};
    for(int i=0; i<rdr->GetOutput()->GetNumberOfPoints(); i++)
    {
        double *pt = rdr->GetOutput()->GetPoint(i);
        cntr[0] += pt[0];
        cntr[1] += pt[1];
        cntr[2] += pt[2];
    }
    cntr[0] /= rdr->GetOutput()->GetNumberOfPoints();
    cntr[1] /= rdr->GetOutput()->GetNumberOfPoints();
    cntr[2] /= rdr->GetOutput()->GetNumberOfPoints();

    //get the bounding box
    double *bounds = rdr->GetOutput()->GetBounds();
    double max_size = 2*std::max(bounds[5]-bounds[4], std::max(bounds[1]-bounds[0],bounds[3]-bounds[2]));
    
    for(int i=0; i<rdr->GetOutput()->GetNumberOfPoints(); i++)
    {
        double pt[3];
        rdr->GetOutput()->GetPoints()->GetPoint(i,pt);
        pt[0] = (pt[0]-cntr[0])/max_size;
        pt[1] = (pt[1]-cntr[1])/max_size;
        pt[2] = (pt[2]-cntr[2])/max_size;
        rdr->GetOutput()->GetPoints()->SetPoint(i,pt);
    }

/*    
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetInputData( rdr->GetOutput() );
    wr->SetFileName("shifted.vtk");
    wr->Write();
*/    
    
    
    z.SetMesh(rdr->GetOutput());
    z.SetOrder(order); //calculate mesh area
    z.CalculateSurfaceGeometricMoments();
    z.SaveGeometricMoments("geom.dat",false);
    z.CalculateZernikeMoments();
    z.SaveZernikeMoments("zern.dat",false);
    
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::Take( z.Reconstruct(-1,1,-1,1,-1,1,0.1) );
    
    vtkSmartPointer<vtkPolyDataWriter> wr = vtkSmartPointer<vtkPolyDataWriter>::New();
    wr->SetInputData( pd );
    wr->SetFileName("reconstr.vtk");
    wr->Write();
    
    
    
    return 0;   
}




bool TestZernikeMoments(const char* meshname)
{ 
    bool result = TestGeometricMoments(meshname);
    
    std::cout<<"testing zernike moments"<<std::endl;
    
    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(meshname);
    rdr->Update();
    
    Zernike z;
    z.SetMesh(rdr->GetOutput());
    z.SetOrder(1); 
    z.CalculateSurfaceGeometricMoments();
    z.CalculateZernikeMoments();

    
    double G000=z.GetGeometricMoment(0,0,0);
    double G100=z.GetGeometricMoment(1,0,0);
    double G010=z.GetGeometricMoment(0,1,0);
    std::complex<double> Z111=z.GetZernikeMoment(1,1,1);
    std::complex<double> Z000=z.GetZernikeMoment(0,0,0);
    
    double G001 = z.GetGeometricMoment(0,0,1);
    std::complex<double> Z011 = z.GetZernikeMoment(0,1,1);
    
    
    std::complex<double> GG = G100-G010*1.0i;
    result = result && std::abs(Z000 - (3.0/(4.0*M_PI)*G000))<1.e-10 &&
        std::abs( Z111 - (3.0/(4.0*M_PI)*std::sqrt(2.5)*GG) )<1.e-10 &&
        std::abs( Z011 - (3.0/(4.0*M_PI)*std::sqrt(2.5)*G001) )<1.e-10;
  
    if( result )
        std::cout<<"Testing zernike moments successful"<<std::endl;
    else
        std::cout<<"Testing zernike moments failed"<<std::endl;

    return result;
}


bool TestGeometricMoments(const char* meshname)
{ 
    std::cout<<"testing geometric moments"<<std::endl;
    
    Zernike z0; //0 moment -- area
    Zernike z1; //1st moment -- centroid
    
    
    z0.RunTests();
    
    vtkSmartPointer<vtkPolyDataReader> rdr = vtkSmartPointer<vtkPolyDataReader>::New();
    rdr->SetFileName(meshname);
    rdr->Update();
    
    z0.SetMesh(rdr->GetOutput());
    z0.SetOrder(0); //calculate mesh area
    z0.CalculateSurfaceGeometricMoments();

    z1.SetMesh(rdr->GetOutput());
    z1.SetOrder(1); //calculate centroid
    z1.CalculateSurfaceGeometricMoments();


    vtkPolyData* mesh =rdr->GetOutput();
    double total_area = 0;
    for(int i=0; i<mesh->GetNumberOfCells(); i++)
    {
        vtkCell* cell = mesh->GetCell(i);
        const double area = dynamic_cast<vtkTriangle*>(cell)->ComputeArea();
        total_area += area;
    }
    
    double centroid[] = {0,0,0};
    for(int i=0; i<mesh->GetNumberOfPoints(); i++)
    {
        double *p = mesh->GetPoint(i);
        centroid[0] += p[0];
        centroid[1] += p[1];
        centroid[2] += p[2];
    }
    
    centroid[0] /= mesh->GetNumberOfPoints();
    centroid[1] /= mesh->GetNumberOfPoints();
    centroid[2] /= mesh->GetNumberOfPoints();
    

    const double area_diff = total_area-z0.GetGeometricMoment(0,0,0);
    const double c1_diff = centroid[0]-z1.GetGeometricMoment(1,0,0)/total_area;
    const double c2_diff = centroid[1]-z1.GetGeometricMoment(0,1,0)/total_area;
    const double c3_diff = centroid[2]-z1.GetGeometricMoment(0,0,1)/total_area;
    
    bool result= area_diff<1e-15 && c1_diff<1e-15 && c2_diff<1e-15 && c3_diff<1e-15;
    
    if( result )
        std::cout<<"Testing geometric moments successful"<<std::endl;
    else
        std::cout<<"Testing geometric moments failed"<<std::endl;

    
    return result;
}

/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkImageData.h>
#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char** argv)
{
    const char* infile = argv[1];
    const char* outfile = argv[2];
    
    vtkSmartPointer<vtkDataSetReader> rdr = vtkSmartPointer<vtkDataSetReader>::New();
    rdr->SetFileName(infile);
    rdr->Update();
    vtkImageData* image = dynamic_cast<vtkImageData*>(rdr->GetOutput());
    
    std::ofstream file;
    file.open(outfile,ios::binary);
    
    std::stringstream header;
    
    header <<"#INRIMAGE-4#{"<<std::endl
        <<"XDIM="<<image->GetDimensions()[0]<<std::endl
        <<"YDIM="<<image->GetDimensions()[1]<<std::endl
        <<"ZDIM="<<image->GetDimensions()[2]<<std::endl
        <<"VDIM=1"<<std::endl
        <<"TYPE=unsigned fixed"<<std::endl
        <<"PIXSIZE=8 bits"<<std::endl
        <<"SCALE=2**0"<<std::endl
        <<"CPU=decm"<<std::endl
        <<"VX="<<image->GetSpacing()[0]<<std::endl
        <<"VY="<<image->GetSpacing()[1]<<std::endl
        <<"VZ="<<image->GetSpacing()[2]<<std::endl
        <<"#GEOMETRY=CARTESIAN"<<std::endl;

                
    header  <<std::setw(256-header.str().size())
            <<std::setfill('\n')             
                <<"##}\n";
        
    file<<header.str();
    
    int *dim = image->GetDimensions();
    file.write((char*)image->GetScalarPointer(),dim[0]*dim[1]*dim[2]);
    
    file.close();
    //std::cout<<header.str();
}   

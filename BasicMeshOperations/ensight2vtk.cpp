/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Reads Ensight Gold (ASCII) files and saves them as one VTK per timestep
 */

#include <vtkSmartPointer.h>
#include <vtkEnSightGoldReader.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkFieldData.h>
#include <vtkDataArrayCollection.h>
#include <vtkDataArray.h>
#include <vtkDataSetWriter.h>
#include <vtkEnSightReader.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>

int main(int argc, char** argv)
{
    if( argc<2 )
    {
        std::cout<<"Options: "<<std::endl;
        std::cout<<"-case <filename> -- filename of the ensight casefile"<<std::endl;
        std::cout<<"-o <prefix> -- prefix for the output filenames (will add numbers)"<<std::endl;
        std::cout<<"-t [ASCII|ASCII] -- binary or ascii"<<std::endl;
        return -1;
    }
    
    std::string filename;
    std::string prefix = "data_";
    bool isbinary = false;
    
    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-case") == 0) 
            filename.assign( argv[++c] );
        else if (strcmp(argv[c], "-o") == 0) 
            prefix.assign( argv[++c] );
        else if (strcmp(argv[c], "-t") == 0) 
            if (strcmp(argv[++c], "BIN") == 0) 
                isbinary = true;
    }    
    
        
    vtkSmartPointer<vtkEnSightReader> rdr;
    
    if (!isbinary)
        rdr = vtkSmartPointer<vtkEnSightGoldReader>::New();
    else
        rdr = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
    
    //const char* filename = "/home/costa/Downloads/lvsmooth-nastin_normals/lvsmooth-nastin.ensi.case";
    rdr->SetCaseFileName(filename.c_str());
    rdr->ReadAllVariablesOn();
    
    std::cout<<"Verify if the file is valid: "<<filename<<std::endl;
    if (!rdr->CanReadFile(filename.c_str()))
    {
        std::cout<<"Invalid file: "<<filename<<std::endl;        
    }

    //std::cout<<"Format type:"<<rdr->DetermineEnSightVersion()<<std::endl;
    //if (!rdr->DetermineEnSightVersion()!=vtkEnSightGoldReader::ENSIGHT_GOLD)
    //{
    //    std::cout<<"Not insight gold ASCII format"<<filename<<std::endl;        
    //}

    
    rdr->DebugOff();
    rdr->Update();

    std::cout<<"Number of timesets: "<<rdr->GetTimeSets()->GetNumberOfItems()<<std::endl;
    //std::cout<<"Number of tuples, timeset 0: "<<rdr->GetTimeSets()->GetItem(0)->GetNumberOfTuples()<<std::endl;
    //std::cout<<"Number of components/tuple, timeset 0: "<<rdr->GetTimeSets()->GetItem(0)->GetNumberOfComponents()<<std::endl;
    std::cout<<"Time values for timeset 0: ";
    
    const int n_timesteps = rdr->GetTimeSets()->GetItem(0)->GetNumberOfTuples();
    for(int i=0; i<n_timesteps; i++)
    {
        std::cout<<rdr->GetTimeSets()->GetItem(0)->GetTuple1(i)<<" ";
    }
    std::cout<<std::endl;
    int n = rdr->GetNumberOfOutputPorts();
    
    //std::cout<<"Output class: "<< rdr->GetOutput(0)->GetClassName()<<std::endl;
    std::cout<<"Number of outputs: "<<n<<std::endl;

    vtkMultiBlockDataSet *ds = rdr->GetOutput();
    
    std::cout<<"Number of blocks: "<<ds->GetNumberOfBlocks()<<std::endl;
    std::cout<<"Number of points: "<<ds->GetNumberOfPoints()<<std::endl;
    vtkFieldData* fd = ds->GetFieldData();
    
    std::cout<<"Number of arrays: "<<fd->GetNumberOfArrays()<<std::endl;
    std::cout<<"Number of tuples: "<<fd->GetNumberOfTuples()<<std::endl;

    //try to save each timestep in a separate vtk
    std::cout<<std::endl<<"Saving timesteps in separate vtk files"<<std::endl;
    for(int i=0; i<n_timesteps; i++)
    {
        rdr->SetTimeValue(rdr->GetTimeSets()->GetItem(0)->GetTuple1(i));
        rdr->Update();

        std::cout<<"Number of blocks: "<<rdr->GetOutput()->GetNumberOfBlocks()<<std::endl;
        
        char ending[20];
        sprintf(ending,"%03d.vtk",i);
        std::string name = prefix+ending;
        
        vtkSmartPointer<vtkDataSetWriter> wr = vtkSmartPointer<vtkDataSetWriter>::New();
        wr->SetFileName(name.c_str());
        wr->SetInputData(rdr->GetOutput()->GetBlock(0));
        wr->Write();
    }
    
    //write the times of the files into a text file
    std::string name = prefix+"times.txt";
    
    std::ofstream file(name.c_str());
    for(int i=0; i<n_timesteps; i++)
    {
        file<<rdr->GetTimeSets()->GetItem(0)->GetTuple1(i)<<" ";
    }

    
    return 0;
}


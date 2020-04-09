/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief Reads Ensight Gold (ASCII) files and calcultes uniformity index for each timestep 
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
#include <vtkPointDataToCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkXMLUnstructuredGridReader.h>

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>

double UniformityIndex(vtkUnstructuredGrid* mesh, const char* cellarrayname);


int main(int argc, char** argv)
{
    //auto rd1 = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    //rd1->SetFileName("box.vtu");
    //rd1->Update();
    //vtkUnstructuredGrid *volmesh1 = dynamic_cast<vtkUnstructuredGrid*>( rd1->GetOutput() ); 
    //cout<<"Ui = "<< UniformityIndex( volmesh1,  "CellIds" )<<std::endl;
    //exit(0);

    if( argc<2 )
    {
        std::cout<<"Options: "<<std::endl;
        std::cout<<"-case <filename> -- filename of the ensight casefile of TETRA mesh"<<std::endl;
        std::cout<<"-a <name> -- point data array on which uniformity is calculated"<<std::endl;
        std::cout<<"-l <int> -- timeline of the variable in ensi.case (>=1)"<<std::endl;
        std::cout<<"-o <filename.csv> -- filename for the output index table"<<std::endl;
        std::cout<<"-t [BIN|ASCII] -- binary or ascii ensight case"<<std::endl;
        return -1;
    }
    
    std::string filename;
    std::string output_filename = "index.csv";
    const char *arrayname = nullptr;
    bool isbinary = false;
    int timeline = 0;

    for (int c = 1; c < argc; c++) {
        if (strcmp(argv[c], "-case") == 0) 
            filename.assign( argv[++c] );
        else if (strcmp(argv[c], "-a") == 0) 
            arrayname = argv[++c];
        else if (strcmp(argv[c], "-l") == 0) 
            timeline = atoi(argv[++c])-1;
        else if (strcmp(argv[c], "-o") == 0) 
            output_filename.assign( argv[++c] );
        else if (strcmp(argv[c], "-t") == 0) 
            if (strcmp(argv[++c], "BIN") == 0) 
                isbinary = true;
    }    
    
        
    if(arrayname == nullptr)
    {
        std::cout<<"No array specified"<<std::endl;
        exit(-1);
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
    std::cout<<"Time values for timeset "<<timeline<<": ";
    
    const int n_timesteps = rdr->GetTimeSets()->GetItem(timeline)->GetNumberOfTuples();
    for(int i=0; i<n_timesteps; i++)
    {
        std::cout<<rdr->GetTimeSets()->GetItem(timeline)->GetTuple1(i)<<" ";
    }
    std::cout<<std::endl;
    int n = rdr->GetNumberOfOutputPorts();
    
    std::cout<<"Number of outputs: "<<n<<std::endl;

    vtkMultiBlockDataSet *ds = rdr->GetOutput();
    
    std::cout<<"Number of blocks: "<<ds->GetNumberOfBlocks()<<std::endl;
    std::cout<<"Number of points: "<<ds->GetNumberOfPoints()<<std::endl;
    vtkFieldData* fd = ds->GetFieldData();
            vtkUnstructuredGrid *volmesh = dynamic_cast<vtkUnstructuredGrid*>( rdr->GetOutput()->GetBlock(0) ); 


    std::cout<<"Number of arrays: "<<fd->GetNumberOfArrays()<<std::endl;
    std::cout<<"Number of tuples: "<<fd->GetNumberOfTuples()<<std::endl;


    std::ofstream table(output_filename);
    table<<"Time, UniformityIndex"<<std::endl;

    std::cout<<"Number of timesteps: "<<n_timesteps<<std::endl;

    for(int i=0; i<n_timesteps; i++)
    {
        if ( rdr->GetTimeSets()->GetItem(timeline)->GetTuple1(i)<0 ) continue;

        std::cout<<"==================================================="<<std::endl;
        std::cout<<"Timestep "<<i<<"/"<<n_timesteps<<std::endl;
        std::cout<<"Time "<<rdr->GetTimeSets()->GetItem(timeline)->GetTuple1(i)<<std::endl;
        rdr->SetTimeValue(rdr->GetTimeSets()->GetItem(timeline)->GetTuple1(i));
        rdr->Update();

        //std::cout<<"Number of blocks: "<<rdr->GetOutput()->GetNumberOfBlocks()<<std::endl;
        
        vtkUnstructuredGrid *volmesh = dynamic_cast<vtkUnstructuredGrid*>( rdr->GetOutput()->GetBlock(0) ); 

        std::cout<<"Pointdata to celldata"<<std::endl;
        auto p2c = vtkSmartPointer<vtkPointDataToCellData>::New();
        p2c->AddPointDataArray(arrayname);
        p2c->SetInputData(volmesh);
        p2c->Update();

        std::cout<<"Calculating"<<std::endl;
        double ui = UniformityIndex( dynamic_cast<vtkUnstructuredGrid*>(p2c->GetOutput()),  arrayname);

        table<<rdr->GetTimeSets()->GetItem(timeline)->GetTuple1(i)<<", "<<ui<<std::endl;
    }
    
    
    return 0;
}

double UniformityIndex(vtkUnstructuredGrid* mesh, const char* cellarrayname)
{
    /*
    vs = SUM(Vc)
    fa = SUM( fc*Vc )/vs
    UniformityIndex = 1 - SUM( |fc-fa|*Vc )/ ( 2 |fa| vs )
    */
    std::cout<<"Calculating averages"<<std::endl;

    std::vector<double> volume(mesh->GetNumberOfCells());

    std::cout<<"Array "<<cellarrayname <<std::endl;
    auto scalars = mesh->GetCellData()->GetArray( cellarrayname );
    std::cout<<"Tuples:"<<scalars->GetNumberOfTuples()<<std::endl;

    double vs = 0;
    double fa = 0;
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    for(uint64_t cellid=0; cellid<mesh->GetNumberOfCells(); cellid++)
    {
        mesh->GetCellPoints( cellid,  ptIds );

        if( ptIds->GetNumberOfIds()!=4 )
        {
            std::cout<<"The mesh has non TETRA. Stopping."<<std::endl;
            throw;
        }

        double p1[3], p2[3], p3[3], p4[3];
        mesh->GetPoint(ptIds->GetId(0), p1);
        mesh->GetPoint(ptIds->GetId(1), p2);
        mesh->GetPoint(ptIds->GetId(2), p3);
        mesh->GetPoint(ptIds->GetId(3), p4);
        volume[cellid] = vtkTetra::ComputeVolume ( p1, p2, p3, p4 );
        vs += volume[cellid];
        fa += volume[cellid] * scalars->GetTuple1(cellid);
    }
    fa /= vs;

    cout<<"Total volume: "<<vs<<std::endl;
    cout<<"fa: "<<fa<<std::endl;

    std::cout<<"Calculating index"<<std::endl;
    double ui = 0;
    for(uint64_t cellid=0; cellid<mesh->GetNumberOfCells(); cellid++)
    {
        ui += abs(scalars->GetTuple1(cellid) - fa)*volume[cellid];
    }
    ui = 1 - ui / (2*abs(fa)*vs);
    std::cout<<"UI = "<<ui<<std::endl;

    return ui;
}
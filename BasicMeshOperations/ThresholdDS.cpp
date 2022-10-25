/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ThresholdDS.cpp
 * Author: costa
 *
 * Created on April 3, 2018, 4:59 PM
 */

#include <vtkDataSetReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkDataSet.h>
#include <vtkThreshold.h>
#include <vtkDataArray.h>
#include <vtkSmartPointer.h>

#include <cstdlib>

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    if(argc<5)
    {
        cout<<"Usage: ThresholdDS dataset.vtk output.vtk cellarrayname thresholdlow thresholdhigh"<<endl;
        return -1;
    }
    
    int c=1;
    const char* infile = argv[c++];
    const char* outfile = argv[c++];
    const char* arrayname = argv[c++];
    double tlow = atof(argv[c++]);
    double thigh = atof(argv[c++]);
    
    cout<<"Input file: "<<infile<<endl;
    cout<<"Output: "<<outfile<<endl;
    cout<<"Variable: "<<arrayname<<endl;
    cout<<"T Low: "<<tlow<<endl;
    cout<<"T High: "<<thigh<<endl;
    
    cout<<"Reading mesh"<<endl;
    vtkSmartPointer<vtkDataSetReader> rd =    vtkSmartPointer<vtkDataSetReader> ::New();
    rd->SetFileName(infile);
    rd->Update();
    
    cout<<"Thresholding"<<endl;
    vtkSmartPointer<vtkThreshold> th =     vtkSmartPointer<vtkThreshold> ::New();
    th->SetInputData(rd->GetOutput());
    th->SetInputArrayToProcess(0,0,0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayname);
    th->SetLowerThreshold(tlow);
    th->SetUpperThreshold(thigh);
    th->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);

    th->Update();
    
    cout<<"Saving"<<endl;
    vtkSmartPointer<vtkUnstructuredGridWriter> wr =     vtkSmartPointer<vtkUnstructuredGridWriter> ::New();
    wr->SetFileName(outfile);
    wr->SetInputData( (vtkDataObject*)th->GetOutput());
    wr->SetFileTypeToBinary();
    wr->Write();
    
    
    
    
    return 0;
}


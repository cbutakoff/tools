/*=========================================================================
Copyright (c) Constantine Butakoff
All rights reserved.
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.
=========================================================================*/

/*! \file
    \brief copies scalars from one polydata to another
 */
#include "vtkPolyDataReader.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkDataArray.h"
#include "vtkType.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkCell.h"
#include "VTKCommonTools.h"

int main(int argc, char *argv[]) {

    if (argc < 5) {
        std::cout << "Usage: " << std::endl;
        std::cout << argv[0] << " <source> <target> <property name> <output>" << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << " reading source " << argv[1] << std::endl;
    vtkSmartPointer<vtkPolyData> sourcePd = vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(argv[1]));

    std::cout << " reading target " << argv[2] << std::endl;
    vtkSmartPointer<vtkPolyData> targetPd = vtkSmartPointer<vtkPolyData>::Take(
            CommonTools::LoadShapeFromFile(argv[2]));


    char* property_name = argv[3];
    std::cout << " property name: " << property_name << std::endl;

    bool use_point_data = sourcePd->GetPointData()->GetArray(property_name) != NULL;

    if (use_point_data) {
        std::cout << "Copying pointdata" << std::endl;

        auto target_scalars = targetPd -> GetPointData() -> GetArray(property_name);
        auto source_scalars = sourcePd -> GetPointData() -> GetArray(property_name);

        vtkSmartPointer<vtkPointLocator> in_target_locator = vtkSmartPointer<vtkPointLocator>::New();
        in_target_locator->SetDataSet(targetPd);
        in_target_locator->BuildLocator();

        for (int i = 0; i < sourcePd -> GetNumberOfPoints(); i++) {
            double *p1 = sourcePd -> GetPoint(i);
            double p2[] = {0.0, 0.0, 0.0};
            vtkIdType pointId = in_target_locator->FindClosestPoint(p1); //(p1, p2, cellId, subId, dist);

            target_scalars->SetTuple1(pointId, source_scalars->GetTuple1(i));
        }
    } else {
        std::cout << "Copying celldata" << std::endl;

        auto target_scalars = targetPd -> GetCellData() -> GetArray(property_name);
        auto source_scalars = sourcePd -> GetCellData() -> GetArray(property_name);

        vtkSmartPointer<vtkCellLocator> in_target_locator = vtkSmartPointer<vtkCellLocator>::New();
        in_target_locator->SetDataSet(targetPd);
        in_target_locator->BuildLocator();

        cout << "Built locator" << endl;

        for (int i = 0; i < sourcePd ->GetNumberOfCells(); i++) {
            
            if(i%10000==0) cout<<"Cell "<<i<<"/"<<sourcePd ->GetNumberOfCells()<<"\r"<<std::flush;
                
            double cell_center_p[3];
            double pt[3];
            double weights[3];
            int not_used;
            sourcePd->GetCell(i)->GetParametricCenter(cell_center_p);
            sourcePd->GetCell(i)->EvaluateLocation(not_used, cell_center_p, pt, weights);

            double closestpoint[3];
            vtkIdType cellid;
            int subid;
            double dist2;
            in_target_locator->FindClosestPoint(pt, closestpoint, cellid, subid, dist2);

            target_scalars->SetTuple1(cellid, source_scalars->GetTuple1(i));
        }
        cout<<endl;
    }

    char outputFileName[255];
    strcpy(outputFileName, argv[4]);
    std::cout << " save output to " << outputFileName << std::endl;

    CommonTools::SaveShapeToFile(targetPd, outputFileName);

    return 0;
}

/*

author  Rubén Cárdenes, Nov. 2012 
adapted from Michel Couprie's pink software
*/


#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <assert.h>
#include <mcimage.h>
#include <mccodimage.h>
//#include <mcutil.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#define VERBOSE

/* =============================================================== */
int main(int argc, char **argv)
/* =============================================================== */
{
  index_t rs, cs, ds, N;
  int32_t datatype;
  struct xvimage * image;

  if (argc < 4)
  {
    fprintf(stderr, "usage: %s in.vtk/mhd out.pgm [datatype] \n", argv[0]);
    fprintf(stderr, "       datatype: 1 for byte, 2 for short int, 4 for long int, 5 for float\n");
    exit(1);
  }

  char* input_file = argv[1];
  char* output_file = argv[2];
	
  typedef float                        PixelType;
  typedef itk::Image< PixelType,  3 >  ImageType3D;
  typedef itk::ImageFileReader< ImageType3D >  ReaderType3D;
  typedef ImageType3D::IndexType             IndexType3D;
  typedef itk::ImageRegionIterator<ImageType3D> IteratorType;

  ReaderType3D::Pointer reader3D = ReaderType3D::New();
  reader3D->SetFileName( input_file );	
  reader3D->Update();

  int	max1 = reader3D->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
  int	max2 = reader3D->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  int	max3 = reader3D->GetOutput()->GetLargestPossibleRegion().GetSize()[2];

  ImageType3D::SpacingType spacing = reader3D->GetOutput()->GetSpacing();
  ImageType3D::PointType origin = reader3D->GetOutput()->GetOrigin();
  IteratorType It( reader3D->GetOutput(), reader3D->GetOutput()->GetLargestPossibleRegion() );
	
  printf("input file: %s, dims: %d %d %d\n",input_file,max1,max2,max3);
  	
  rs = max2;
  cs = max1;
  ds = max3;
  N = rs * cs * ds;
  datatype = atoi(argv[3]);
	
//---------------------------------------------------	

  if ((datatype != 1) && (datatype != 2) && (datatype != 4) && (datatype != 5))
  {
    fprintf(stderr, "%s: bad value for pix size: %d\n", argv[0], datatype);
	fprintf(stderr, "usage: %s in.vtk/mhd out.pgm [datatype] \n", argv[0]);
    fprintf(stderr, "       datatype: 1 for byte, 2 for u_short, 4 for long int, 5 for float\n");
    exit(1);
  }

  if (datatype == 1)
  {
	unsigned char* I;
    image = allocimage(NULL, rs, cs, ds, VFF_TYP_1_BYTE);
    if (image == NULL)
    {   fprintf(stderr,"%s : allocimage failed\n", argv[0]);
        exit(1);
    }
	int i=0;
	I = UCHARDATA(image);
	for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		I[i] = (unsigned char)It.Get();
		i++;
	}
  }
  else if (datatype == 2)
  {
    uint16_t * I;
    image = allocimage(NULL, rs, cs, ds, VFF_TYP_2_BYTE);
    if (image == NULL)
    {   fprintf(stderr,"%s : allocimage failed\n", argv[0]);
        exit(1);
    }
	int i=0;
	I = USHORTDATA(image);
	for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		I[i] = (int32_t)It.Get();
		i++;
	}
  }
  else if (datatype == 4)
  {
	int32_t *I;
    image = allocimage(NULL, rs, cs, ds, VFF_TYP_4_BYTE);
    if (image == NULL)
    {   fprintf(stderr,"%s : allocimage failed\n", argv[0]);
        exit(1);
    }
	int i=0;
	I = SLONGDATA(image);
	for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		  I[i] = (int)It.Get();
		  i++;
	}
  }
  else if (datatype == 5)
  {
	  float* I; 
    image = allocimage(NULL, rs, cs, ds, VFF_TYP_FLOAT);
    if (image == NULL)
    {   fprintf(stderr,"%s : allocimage failed\n", argv[0]);
        exit(1);
    }
	int i=0;
	  I = FLOATDATA(image);
	for( It.GoToBegin(); !It.IsAtEnd(); ++It  ) {
		I[i] = (float)It.Get();
		i++;
	}
  }

  image->xdim = spacing[0];
  image->ydim = spacing[1];
  image->zdim = spacing[2];
	
  image->origin_x = origin[0];
  image->origin_y = origin[1];
  image->origin_z = origin[2];
	
  writeimage(image, output_file);
  freeimage(image);

  return 0;
} /* main */



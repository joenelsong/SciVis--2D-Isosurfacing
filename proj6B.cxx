/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SpecularSpheres.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//
// This examples demonstrates the effect of specular lighting.
//



#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include "tricase.cxx"
#include "TriangleList.h"


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D

    return dims[0]*dims[1]*dims[2];
    // 2D
    //return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    //3D
    return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    //return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    //3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    //return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    //3D
    return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    //return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    idx[0] = pointId%dims[0];
    idx[1] = (pointId/dims[0])%dims[1];
    idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    //idx[0] = pointId%dims[0];
    //idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    //idx[0] = cellId%(dims[0]-1);
    //idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: EvaluateVectorFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a vector field defined on the mesh.  Its size is 2*dims[0]*dims[1].
//        The first value in the field is the x-component for the first point.
//        The second value in the field is the y-component for the first point.
//
//     rv (output): the interpolated field value. (0,0) if the location is out of bounds.
//
// ****************************************************************************

class SegmentList
{
   public:
                   SegmentList() { maxSegments = 10000; segmentIdx = 0; pts = new float[4*maxSegments]; };
     virtual      ~SegmentList() { delete [] pts; };

     void          AddSegment(float X1, float Y1, float X2, float Y2);
     vtkPolyData  *MakePolyData(void);

   protected:
     float        *pts;
     int           maxSegments;
     int           segmentIdx;
};

void
SegmentList::AddSegment(float X1, float Y1, float X2, float Y2)
{
    pts[4*segmentIdx+0] = X1;
    pts[4*segmentIdx+1] = Y1;
    pts[4*segmentIdx+2] = X2;
    pts[4*segmentIdx+3] = Y2;
    segmentIdx++;
}

vtkPolyData *
SegmentList::MakePolyData(void)
{
    int nsegments = segmentIdx;
    int numPoints = 2*(nsegments);
    vtkPoints *vtk_pts = vtkPoints::New();
    vtk_pts->SetNumberOfPoints(numPoints);
    int ptIdx = 0;
    vtkCellArray *lines = vtkCellArray::New();
    lines->EstimateSize(numPoints,2);
    for (int i = 0 ; i < nsegments ; i++)
    {
        double pt[3];
        pt[0] = pts[4*i];
        pt[1] = pts[4*i+1];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx, pt);
        pt[0] = pts[4*i+2];
        pt[1] = pts[4*i+3];
        pt[2] = 0.;
        vtk_pts->SetPoint(ptIdx+1, pt);
        vtkIdType ids[2] = { ptIdx, ptIdx+1 };
        lines->InsertNextCell(2, ids);
        ptIdx += 2;
    }

    vtkPolyData *pd = vtkPolyData::New();
    pd->SetPoints(vtk_pts);
    pd->SetLines(lines);
    lines->Delete();
    vtk_pts->Delete();

    return pd;
}

unsigned int IdentifyCase(unsigned int v0, unsigned int v1, unsigned int v2, unsigned int v3, 
    unsigned int v4, unsigned int v5, unsigned int v6, unsigned int v7)
{
	unsigned int caseNumber = v0 + (v1<<1) + (v2<<2) + (v3<<3) + (v4<<4) + (v5<<5) + (v6<<6) + (v7<<7);
    //printf("case number: %d\n",caseNumber);
    return caseNumber;
}
bool InterpolatePointLocation(int edgeIdx, float *A, float *B, float valA, float valB, float isoVal, float *OutputCoords)
{
	float t; //d
    t = (isoVal - valA) / (valB-valA) ;

	if(edgeIdx == 0 || edgeIdx == 2 || edgeIdx == 4 || edgeIdx == 6) { // X interpolation
		OutputCoords[0] = A[0] + t*(B[0]-A[0]);
		if (A[1] != B[1]) { fprintf(stderr,"Error: A[1] != B[1]"); throw 20; }
		OutputCoords[1] = A[1]; // A[1] and B[1] should be the same...
        if (A[2] != B[2]) { fprintf(stderr,"Error: A[2] != B[2]"); throw 20; }
        OutputCoords[2] = A[2]; // A[1] and B[1] should be the same...
	}
	else if(edgeIdx == 1 || edgeIdx == 3 || edgeIdx == 5 || edgeIdx == 7) // Y Interpolation
	{
		if (A[0] != B[0]) { fprintf(stderr,"Error: A[0] != B[0]"); throw 20; }
		OutputCoords[0] = A[0]; // A[0] and B[0] should be the same...
        if (A[2] != B[2]) { fprintf(stderr,"Error: A[2] != B[2]"); throw 20; }
        OutputCoords[2] = A[2]; // A[0] and B[0] should be the same...
		OutputCoords[1] = A[1] + t*(B[1]-A[1]);
	}
    else if(edgeIdx == 8 || edgeIdx == 9 || edgeIdx == 10 || edgeIdx == 11) // Z Interpolation
    {
        if (A[0] != B[0]) { fprintf(stderr,"Error: A[0] != B[0]"); throw 20; }
        OutputCoords[0] = A[0]; // A[1] and B[1] should be the same...
        if (A[1] != B[1]) { fprintf(stderr,"Error: A[1] != B[1]"); throw 20; }
        OutputCoords[1] = A[1]; // A[1] and B[1] should be the same...
        OutputCoords[2] = A[2] + t*(B[2]-A[2]);
    }
	else {
        fprintf(stderr,"INVALID EDGE!\n");
		return false; // Invalid Edge Error
        
	}
	return true;
}

int main()
{
    int  i, j, k;
    //throw 21;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj6B.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims); // wee

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    TriangleList tl;
    //tl.AddTriangle(-10, -10, +10, +10, -10, +10, +10, +10, +10); // Add segment (-10,-10) -> (+10, -10)
    //tl.AddTriangle(-10, -10, +10, -10, +10, +10, +10, +10, +10);
    
// YOUR CODE TO GENERATE ISOLINES SHOULD GO HERE!
    // Use lookup table from tricase.cxx
    const float ISOVALUE = 3.2; // Given by instructions
    //printf("tricase[0][0]:%d\n", triCase[0][0]);
    
    // Iterate through cells starting at 0 and going to numCells
    int x = dims[0];    int y = dims[1];     int z= dims[2];
    int numCells = (x-1) * (y-1) * (z-1);
    float v0, v1, v2, v3, v4, v5, v6, v7 = -99; // Corners 
        		
    unsigned int v0Flag, v1Flag, v2Flag, v3Flag, v4Flag, v5Flag, v6Flag, v7Flag; // flags = t true if corner is above ISOVALUE
    v0Flag=v1Flag=v2Flag=v3Flag=v4Flag=v5Flag=v6Flag=v7Flag=-99; // Initialize to crazy value to detect errors
    // DEBUG CODE //
    int FlagCounts[8] = {0, 0, 0, 0, 0, 0, 0, 0}; // Initialize to crazy value to detect errors
    int caseCount[256];
    int cc;
    for (cc = 0; cc<256; cc++) // Initialize case counts to 0
        caseCount[cc] = 0;

    //printf("before loop \n");
    for (k = 0; k < (z-1); k++) {
        for (j = 0; j < (y-1); j++) { // j = y index
            for (i = 0; i < (x-1); i++ )    { // i = x index}
                // Get Fvalues for each vertex of the cell

                v0 =  F [ x*j + i + k*x*y ];
                v1 =  F [ x*j + i + k*x*y +1 ];
                v2 =  F [ x*(j+1) + i + k*x*y ];
                v3 =  F [ x*(j+1) + i + k*x*y + 1];
                v4 =  F [ x*j + i + ((k+1)*x*y) ];
                v5 =  F [ x*j + i + ((k+1)*x*y) +1 ];
                v6 =  F [ x*(j+1) + i + ((k+1)*x*y) ];
                v7 =  F [ x*(j+1) + i + ((k+1)*x*y) + 1];

        		// Set Corner Flags for Case identification
        		v0Flag = (v0 >= ISOVALUE) ? 1 : 0;
        		v1Flag = (v1 >= ISOVALUE) ? 1 : 0;
        		v2Flag = (v2 >= ISOVALUE) ? 1 : 0;
        		v3Flag = (v3 >= ISOVALUE) ? 1 : 0;
                v4Flag = (v4 >= ISOVALUE) ? 1 : 0;
                v5Flag = (v5 >= ISOVALUE) ? 1 : 0;
                v6Flag = (v6 >= ISOVALUE) ? 1 : 0;
                v7Flag = (v7 >= ISOVALUE) ? 1 : 0;
        		
        		// Look up table for precomputed answers to all cases
        		unsigned int icase = IdentifyCase(v0Flag, v1Flag, v2Flag, v3Flag, v4Flag, v5Flag, v6Flag, v7Flag);
                //if (icase > 200) printf("icase: %d", icase);
    		//if (brFlag && !blFlag && !tlFlag && !trFlag) {
        			//printf("Vertice Flags: %d %d %d %d\n", blFlag, brFlag, tlFlag, trFlag);
        		if (icase <0) printf("icase ERROR: = %d\n", icase);
                if (v0Flag) FlagCounts[0]++; if (v1Flag) FlagCounts[1]++; if (v2Flag) FlagCounts[2]++; if (v3Flag) FlagCounts[3]++;
                if (v4Flag) FlagCounts[4]++; if (v5Flag) FlagCounts[5]++; if (v6Flag) FlagCounts[6]++; if (v7Flag) FlagCounts[7]++;
        		//}
        		
        		
    	//if (icase != 10)	{continue; } // Test Case by Cases
        		// Working: {}
        		// Probably Working: {}
        		// Not Working: {}
    						
        		caseCount[icase]++;

    		    int edge = -99;
    		    //int edge2 = -99;
    		    float pt1[3][3]; // interpolated position on edge 3 coordinates and 3 points to make a triangle

    		    int s=0;
    		    while( triCase[icase][s*3] != -1)
    		    {
    		    	float A[3];
    		    	float B[3];
    		    	float valA =-99;
    		    	float valB =-99;

    		    	int p;
    		    	for (p =0; p <3; p++) // Iterate Three Times to Interpolate 3 points so we can draw a triangle
    		    	{
    			    	edge = triCase[icase][3*s+p]; // change
    		    	//printf("edge = %d\n", edge);
    			    	switch (edge) {
    			    		case 0:
    			    			A[0] = X[i];      A[1] = Y[j];      A[2] = Z[k];
    			    			B[0] = X[i+1];    B[1] = Y[j];      B[2] = Z[k];
    			    			valA=v0;
    			    			valB=v1;
    			    			break;
    			    		case 1:
    			    			A[0] = X[i+1];    A[1] = Y[j];      A[2] = Z[k];
                                B[0] = X[i+1];    B[1] = Y[j+1];    B[2] = Z[k];
                                valA=v1;
                                valB=v3;
                                break;
    			    		case 2:
    			    			A[0] = X[i];      A[1] = Y[j+1];    A[2] = Z[k];
                                B[0] = X[i+1];    B[1] = Y[j+1];    B[2] = Z[k];
                                valA=v2;
                                valB=v3;
                                break;
    			    		case 3:
    			    			A[0] = X[i];      A[1] = Y[j];      A[2] = Z[k];
                                B[0] = X[i];      B[1] = Y[j+1];    B[2] = Z[k];
                                valA=v0;
                                valB=v2;
                                break;
                            case 4:
                                A[0] = X[i];      A[1] = Y[j];      A[2] = Z[k+1];
                                B[0] = X[i+1];    B[1] = Y[j];      B[2] = Z[k+1];
                                valA=v4;
                                valB=v5;
                                break;
                            case 5:
                                A[0] = X[i+1];    A[1] = Y[j];      A[2] = Z[k+1];
                                B[0] = X[i+1];    B[1] = Y[j+1];    B[2] = Z[k+1];
                                valA=v5;
                                valB=v7;
                                break;
                            case 6:
                                A[0] = X[i];      A[1] = Y[j+1];    A[2] = Z[k+1];
                                B[0] = X[i+1];    B[1] = Y[j+1];    B[2] = Z[k+1];
                                valA=v6;
                                valB=v7;
                                break;
                            case 7:
                                A[0] = X[i];      A[1] = Y[j];      A[2] = Z[k+1];
                                B[0] = X[i];      B[1] = Y[j+1];    B[2] = Z[k+1];
                                valA=v4;
                                valB=v6;
                                break;
                            case 8:
                                A[0] = X[i];      A[1] = Y[j];      A[2] = Z[k];
                                B[0] = X[i];      B[1] = Y[j];      B[2] = Z[k+1];
                                valA=v0;
                                valB=v4;
                                break;
                            case 9:
                                A[0] = X[i+1];    A[1] = Y[j];      A[2] = Z[k];
                                B[0] = X[i+1];    B[1] = Y[j];      B[2] = Z[k+1];
                                valA=v1;
                                valB=v5;
                                break;
                            case 10:
                                A[0] = X[i];      A[1] = Y[j+1];    A[2] = Z[k];
                                B[0] = X[i];      B[1] = Y[j+1];    B[2] = Z[k+1];
                                valA=v2;
                                valB=v6;
                                break;
                            case 11:
                                A[0] = X[i+1];    A[1] = Y[j+1];    A[2] = Z[k];
                                B[0] = X[i+1];    B[1] = Y[j+1];    B[2] = Z[k+1];
                                valA=v3;
                                valB=v7;
                                break;
                        }

    		    		InterpolatePointLocation(edge, A, B, valA, valB, ISOVALUE, pt1[p]);
    		    	}

                    int t, t2;
                    for (t =0; t<3; t++) {
                        for (t2=0; t2<3; t2++) {
                            if ( pt1[t][t2] > 10 ) 
                                printf("Large Value for dim[%d] = %f\n",t, pt1[t][t2]);
                        }
                    }
    		    	tl.AddTriangle(pt1[0][0], pt1[0][1], pt1[0][2], pt1[1][0], pt1[1][1], pt1[1][2], pt1[2][0], pt1[2][1], pt1[2][2]);
                   // printf("Triangle Added!\n");
    		    	//sl.AddSegment(-10, -10, +10, -10);
                    s++; // increment While Loop
    		    }

        		// ~TestCode- should create diagonal lines from bottom left to top right of every cell
        		//printf("bl=%f4, br=%f4, tl=%f4, tr=%f4",X[i], Y[j], X[i+1], Y[j+1] );
        		//sl.AddSegment(X[i], Y[j], X[i+1], Y[j+1]);

        	}
        }
    }
    
    
    printf("Case Counts\n");
    int pc = 0;
    for (pc = 0; pc < 256; pc++) {
        //if (caseCount[pc] < 400)
            printf("%d:%d\n  ",pc, caseCount[pc]);
    }
    printf("Flag Counts\n");
    for (pc = 0; pc < 8; pc++)
        printf("%d:%d\n  ",pc, FlagCounts[pc]);
    
    
    vtkPolyData *pd = tl.MakePolyData();

    //This can be useful for debugging
/*
    vtkDataSetWriter *writer = vtkDataSetWriter::New();
    writer->SetFileName("paths.vtk");
    writer->SetInputData(pd);
    writer->Write();
 */

    vtkSmartPointer<vtkDataSetMapper> win1Mapper =
      vtkSmartPointer<vtkDataSetMapper>::New();
    win1Mapper->SetInputData(pd);
    win1Mapper->SetScalarRange(0, 0.15);

    vtkSmartPointer<vtkActor> win1Actor =
      vtkSmartPointer<vtkActor>::New();
    win1Actor->SetMapper(win1Mapper);

    vtkSmartPointer<vtkRenderer> ren1 =
      vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renWin =
      vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(ren1);

    vtkSmartPointer<vtkRenderWindowInteractor> iren =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);
    ren1->AddActor(win1Actor);
    ren1->SetBackground(0.0, 0.0, 0.0);
    renWin->SetSize(800, 800);

    ren1->GetActiveCamera()->SetFocalPoint(0,0,0);
    ren1->GetActiveCamera()->SetPosition(0,0,50);
    ren1->GetActiveCamera()->SetViewUp(0,1,0);
    ren1->GetActiveCamera()->SetClippingRange(20, 120);
    ren1->GetActiveCamera()->SetDistance(30);

    // This starts the event loop and invokes an initial render.
    //
    iren->Initialize();
    iren->Start();

    pd->Delete();
}

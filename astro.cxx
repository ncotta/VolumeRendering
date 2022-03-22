/*
Author(s): Niklaas Cotta, Hank Childs
All code without a specified author was written by Niklaas Cotta
Date Created: 3/07/22
Description: Volume Rendering using VTK
*/

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
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkContourFilter.h>
#include <vtkRectilinearGrid.h>
#include <vtkDataSetWriter.h>
#include <vtkRectilinearGridToTetrahedra.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <iostream>
#include <cmath>

using namespace std;

/* ================== Camera Implementation ================== */

// written by Hank Childs
int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
}

// written by Hank Childs
struct Camera
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};

// written by Hank Childs
Camera
SetupCamera(void)
{
    Camera rv;
    rv.focus[0] = 0;
    rv.focus[1] = 0;
    rv.focus[2] = 0;
    rv.up[0] = 0;
    rv.up[1] = -1;
    rv.up[2] = 0;
    rv.angle = 30;
    rv.near = 7.5e+7;
    rv.far = 1.4e+8;
    rv.position[0] = -8.25e+7;
    rv.position[1] = -3.45e+7;
    rv.position[2] = 3.35e+7;

    return rv;
}

/* ================== TransferFunction Implementation ================== */

struct TransferFunction
{
    double          min;
    double          max;
    int             numBins;
    unsigned char  *colors;    // size is 3*numBins
    double         *opacities; // size is numBins

    int GetBin(double value)
    {
        int result = -1;

        if ((min <= value) && (value <= max))
        {
            result = numBins * (value - min) / (max - min);
        }
        
        return result;
    }

    void ApplyTransferFunction(double value, unsigned char *RGB, double &opacity)
    {
        int bin = GetBin(value);
        if (bin != -1)  // bin is valid
        {
            // written by Hank Childs
            RGB[0] = colors[3*bin+0];
            RGB[1] = colors[3*bin+1];
            RGB[2] = colors[3*bin+2];
            opacity = opacities[bin];
        }
        else 
        {
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
            opacity = 0;
        }
    }
};

// written by Hank Childs
TransferFunction
SetupTransferFunction(void)
{
    int  i;

    TransferFunction rv;
    rv.min = 10;
    rv.max = 15;
    rv.numBins = 256;
    rv.colors = new unsigned char[3*256];
    rv.opacities = new double[256];
    unsigned char charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    for (i = 0 ; i < 256 ; i++)
        rv.opacities[i] = charOpacity[i]/255.0;
    const int numControlPoints = 8;
    unsigned char controlPointColors[numControlPoints*3] = { 
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0, 
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76 
       };
    double controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
    for (i = 0 ; i < numControlPoints-1 ; i++)
    {
        int start = controlPointPositions[i]*rv.numBins;
        int end   = controlPointPositions[i+1]*rv.numBins+1;
        //cerr << "Working on " << i << "/" << i+1 << ", with range " << start << "/" << end << endl;
        if (end >= rv.numBins)
            end = rv.numBins-1;
        for (int j = start ; j <= end ; j++)
        {
            double proportion = (j/(rv.numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (int k = 0 ; k < 3 ; k++)
                rv.colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }    

    return rv;
}

// Helper function to find a ray for a given pixel
void findRay(Camera rv, double raydir[3], int W, int H)
{
    double look[3], u[3], v[3], dx[3], dy[3];
    double radAngle = (M_PI * rv.angle ) / 180.0;  // convert degrees to radians
    double length = 500.0;

    // Obtain look vector
    look[0] = rv.focus[0] - rv.position[0];
    look[1] = rv.focus[1] - rv.position[1];
    look[2] = rv.focus[2] - rv.position[2];

    // Normalize look
    double lookMag = sqrt((look[0] * look[0]) + (look[1] * look[1]) + (look[2] * look[2]));
    double lookNorm[3];
    lookNorm[0] = look[0] / lookMag;
    lookNorm[1] = look[1] / lookMag;
    lookNorm[2] = look[2] / lookMag;

    // Normalize u
    double lookUp[3];
    lookUp[0] = look[1] * rv.up[2] - look[2] * rv.up[1];
    lookUp[1] = look[2] * rv.up[0] - look[0] * rv.up[2];
    lookUp[2] = look[0] * rv.up[1] - look[1] * rv.up[0];

    double lookUpMag = sqrt((lookUp[0] * lookUp[0]) + (lookUp[1] * lookUp[1]) + (lookUp[2] * lookUp[2]));
    u[0] = lookUp[0] / lookUpMag;
    u[1] = lookUp[1] / lookUpMag;
    u[2] = lookUp[2] / lookUpMag;

    // Normalize v
    double lookU[3];
    lookU[0] = look[1] * u[2] - look[2] * u[1];
    lookU[1] = look[2] * u[0] - look[0] * u[2];
    lookU[2] = look[0] * u[1] - look[1] * u[0];

    double lookUMag = sqrt((lookU[0] * lookU[0]) + (lookU[1] * lookU[1]) + (lookU[2] * lookU[2]));
    v[0] = lookU[0] / lookUMag;
    v[1] = lookU[1] / lookUMag;;
    v[2] = lookU[2] / lookUMag;;

    // Get dx, dy
    dx[0] = (2.0 * tan(radAngle / 2.0)) / length * u[0];
    dx[1] = (2.0 * tan(radAngle / 2.0)) / length * u[1];
    dx[2] = (2.0 * tan(radAngle / 2.0)) / length * u[2];

    dy[0] = (2.0 * tan(radAngle / 2.0)) / length * v[0];
    dy[1] = (2.0 * tan(radAngle / 2.0)) / length * v[1];
    dy[2] = (2.0 * tan(radAngle / 2.0)) / length * v[2];

    // Get d(i, j)
    raydir[0] = lookNorm[0] + ((((2.0 * (double)W) + 1.0 - length) / 2.0) * dx[0]) + ((((2.0 * (double)H) + 1.0 - length) / 2.0) * dy[0]);
    raydir[1] = lookNorm[1] + ((((2.0 * (double)W) + 1.0 - length) / 2.0) * dx[1]) + ((((2.0 * (double)H) + 1.0 - length) / 2.0) * dy[1]);
    raydir[2] = lookNorm[2] + ((((2.0 * (double)W) + 1.0 - length) / 2.0) * dx[2]) + ((((2.0 * (double)H) + 1.0 - length) / 2.0) * dy[2]);
}

// Helper function to intersect the volume data with a given ray
void intersectVolume(Camera rv, double raydir[3], double samples[][3], int nsamples)
{
    int stepsize = round((rv.far - rv.near) / (nsamples - 1));  // 255 samples
    double distance = rv.near;

    for (int i = 0; i < nsamples; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            samples[i][j] = rv.position[j] + (raydir[j] * distance);
        }
        distance += (double)stepsize;
    }
}

// Helper function to interpolate 3D data given a point pt (via trilinear interpolation)
double interpolate(const double *pt, const int *dims,
                   const float *X, const float *Y, 
                   const float *Z, const float *F)
{
    int idx[3] = {-1, -1, -1};
    int x0, x1, y0, y1, z0, z1;
    double c000, c001, c010, c011, c100, c101, c110, c111;
    double c00, c01, c10, c11;
    double c0, c1, c;
    int indexArray[3];
    double xd, yd, zd;

    // Step 1) Find cell that contains point pt
    for (int i = 0; i < dims[0]; i++)
    {   
        if ((pt[0] >= X[i]) && (pt[0] <= X[i+1]))
        {
            idx[0] = i;
            break;
        }
    }

    for (int j = 0; j < dims[1]; j++)
    {
        if ((pt[1] >= Y[j]) && (pt[1] <= Y[j+1]))
        {
            idx[1] = j;
            break;
        }
    }

    for (int k = 0; k < dims[2]; k++)
    {
        if ((pt[2] >= Z[k]) && (pt[2] <= Z[k+1]))
        {
            idx[2] = k;
            break;
        } 
    }

    if ((idx[0] == -1) || (idx[1] == -1) || (idx[2] == -1))
    {
        return 0;
    }

    // Step 2) Get coords
    x0 = idx[0];
    x1 = idx[0] + 1;
    y0 = idx[1];
    y1 = idx[1] + 1;
    z0 = idx[2];
    z1 = idx[2] + 1;

    // Step 3) Get points (xzy)
    indexArray[0] = x0;
    indexArray[1] = y0;
    indexArray[2] = z0;
    c000 = F[GetPointIndex(indexArray, dims)];
    
    indexArray[0] = x1;
    indexArray[1] = y0;
    indexArray[2] = z0;
    c100 = F[GetPointIndex(indexArray, dims)];

    indexArray[0] = x1;
    indexArray[1] = y0;
    indexArray[2] = z1;
    c101 = F[GetPointIndex(indexArray, dims)];

    indexArray[0] = x0;
    indexArray[1] = y0;
    indexArray[2] = z1;
    c001 = F[GetPointIndex(indexArray, dims)];

    indexArray[0] = x0;
    indexArray[1] = y1;
    indexArray[2] = z0;
    c010 = F[GetPointIndex(indexArray, dims)];

    indexArray[0] = x1;
    indexArray[1] = y1;
    indexArray[2] = z0;
    c110 = F[GetPointIndex(indexArray, dims)];

    indexArray[0] = x1;
    indexArray[1] = y1;
    indexArray[2] = z1;
    c111 = F[GetPointIndex(indexArray, dims)];

    indexArray[0] = x0;
    indexArray[1] = y1;
    indexArray[2] = z1;
    c011 = F[GetPointIndex(indexArray, dims)];

    xd = (pt[0] - X[x0]) / (X[x1] - X[x0]);
    yd = (pt[1] - Y[y0]) / (Y[y1] - Y[y0]);
    zd = (pt[2] - Z[z0]) / (Z[z1] - Z[z0]);

    c00 = c000 * (1 - xd) + c100 * xd;
    c01 = c001 * (1 - xd) + c101 * xd;
    c10 = c010 * (1 - xd) + c110 * xd;
    c11 = c011 * (1 - xd) + c111 * xd;

    c0 = c00 * (1 - yd) + c10 * yd;
    c1 = c01 * (1 - yd) + c11 * yd;

    c = c0 * (1 - zd) + c1 * zd;

    return c;
}

// written by Hank Childs
void WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

// written by Hank Childs
vtkImageData *NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return image;
}


/* ================== Main ================== */

int main()
{
    // Read from file
    vtkDataSetReader *reader = vtkDataSetReader::New();
    reader->SetFileName("astro512.vtk");
    reader->Update();

    // Rectilinear grid
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) reader->GetOutput();

    int dims[3];
    rgrid->GetDimensions(dims);

    // 3D Data 
    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *Z = (float *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);

    int ncells = rgrid->GetNumberOfCells();

    Camera rv = SetupCamera();
    TransferFunction tf = SetupTransferFunction();

    // Change these two values to change window size and number of samples per ray
    int length = 500;  // window size is length x length
    int nsamples = 1024;  // number of samples per ray
    double raydir[3];
    double samples[nsamples][3];
    double samplevals[nsamples];

    vtkImageData *image = NewImage(length, length);
    unsigned char *pixel;

    // for each pixel
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            // Step 1: Find ray for pixel
            findRay(rv, raydir, i, j);                

            // Step 2: Intersect volume with ray (coordinates)       
            intersectVolume(rv, raydir, samples, nsamples); 
            for (int k = 0; k < nsamples; k++)
            {
                // Get the value at those coordinates
                samplevals[k] = interpolate(samples[k], dims, X, Y, Z, F);               
            }
            
            // Step 3: Calculate color from intersection
            unsigned char s0[3], back[3];
            double s0a, backa;
            double total[3];
            double totala;

            tf.ApplyTransferFunction(samplevals[0], s0, s0a);  // seed value
            s0a = 1.0 - pow(1.0 - s0a, (500.0/nsamples));  // opacity correction
            
            total[0] = double(s0[0]/256.0);
            total[1] = double(s0[1]/256.0);
            total[2] = double(s0[2]/256.0);
            totala = s0a;

            for (int l = 1; l < nsamples; l++)
            {
                tf.ApplyTransferFunction(samplevals[l], back, backa);
                backa = 1.0 - pow(1.0 - backa, (500.0/nsamples));

                // Apply Composition
                total[0] = total[0] + (1.0 - totala) * backa * double(back[0]/256.0);
                total[1] = total[1] + (1.0 - totala) * backa * double(back[1]/256.0);
                total[2] = total[2] + (1.0 - totala) * backa * double(back[2]/256.0);
                totala = totala + (1.0 - totala) * backa;
            }

            // Step 4: Assign value to pixel
            pixel = (unsigned char *) image->GetScalarPointer(i,j,0);
            pixel[0] = (unsigned char)(total[0]*256.0);
            pixel[1] = (unsigned char)(total[1]*256.0);
            pixel[2] = (unsigned char)(total[2]*256.0);
        }
    }

    WriteImage(image, "astroOutput");
}

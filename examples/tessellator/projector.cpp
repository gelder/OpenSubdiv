#include "Projector.h"

#include "lava/df/triangleMesh.h"
#include "lava/df/hbv.h"

#include "bedrock/tf/stopwatch.h"

#include "bedrock/gf/vec3d.h"
#include "bedrock/gf/vec3f.h"
#include "bedrock/gf/size3.h"
#include "bedrock/gf/triangle.h"


#include "bedrock/tf/refPtr.h"
#include "bedrock/tf/weakPtr.h"

#include <iostream>
#include <stdio.h>

using namespace std;

class DfHbv;

bool
WeightObject_Projector::Initialize(
    const OpenSubdivTessellatorSharedPtr &tessellator,    
    int subdivLevel) 
{
    _tessellator = tessellator;

    _subdivLevel = subdivLevel;
    
    return true;
}


bool
WeightObject_Projector::ProjectPoints(
    std::vector<WeightObject_Projector::EvalCoords> *pointsProjectedOntoSubdiv,
    const std::vector<GfVec3d> &coarsePoints,    
    const std::vector<GfVec3d> &pointsToProject)
{
    TfStopwatch updateTessellatorTime, createTriangleMeshTime,
        closestPointProjectTime;

    updateTessellatorTime.Start();

    // Two buffers used to store computed values for OpenSubdiv
    OpenSubdiv::OsdCpuVertexBuffer*    vertexBuffer =
        _tessellator->CreateVertexBuffer();
    OpenSubdiv::OsdCpuVertexBuffer*    vvBuffer =
        _tessellator->CreateVVBuffer();

    pointsProjectedOntoSubdiv->clear();
    
    std::string errMsg;

    { // Initialize point positions in tessellator
        std::vector<float> points(coarsePoints.size() * 3);
        for (size_t i=0; i<coarsePoints.size(); ++i) {
            points[i*3] = coarsePoints[i][0];
            points[i*3+1] = coarsePoints[i][1];
            points[i*3+2] = coarsePoints[i][2];            
        }
        _tessellator->SetCoarsePositions(vertexBuffer, points);
    }

    if (not _tessellator->Refine(vertexBuffer, vvBuffer, false, &errMsg)) {
        delete vertexBuffer;
        delete vvBuffer;                            
        TF_WARN(errMsg);
        return false;
    }

    updateTessellatorTime.Stop();

    createTriangleMeshTime.Start();
    
    std::vector<int> quadVerts;
    std::vector<size_t> triVerts;
    std::vector<GfVec3d> quadPoints;
    DfHbv hbv;
    
    if (not _tessellator->GetQuads(&quadVerts, &errMsg)) {
        delete vertexBuffer;
        delete vvBuffer;                    
        TF_WARN(errMsg);
        return false;
    }

    int numQuads = quadVerts.size()/4;
    
    // Six triangle indices for each four quad indices
    triVerts.resize(quadVerts.size() * 1.5);

    // Triangulate quad mesh
    int currentIndex = 0;
    for (size_t i=0; i<numQuads; ++i) {
        int quadStart = i*4;

        triVerts[currentIndex++] = quadVerts[quadStart+0];
        triVerts[currentIndex++] = quadVerts[quadStart+1];
        triVerts[currentIndex++] = quadVerts[quadStart+2];

        triVerts[currentIndex++] = quadVerts[quadStart+2];
        triVerts[currentIndex++] = quadVerts[quadStart+3];
        triVerts[currentIndex++] = quadVerts[quadStart+0];
    }        

    // Extract refined points 
    { 
        std::vector<float> refinedPoints;
        if (not _tessellator->GetPositions(vertexBuffer,
                                          &refinedPoints, &errMsg)) {
            delete vertexBuffer;
            delete vvBuffer;            
            TF_WARN(errMsg);
            return false;
        }
        quadPoints.resize(refinedPoints.size()/3);
        for (size_t i=0; i<quadPoints.size(); ++i) {
            quadPoints[i] = GfVec3d(
                refinedPoints[i*3], refinedPoints[i*3+1], refinedPoints[i*3+2]);
        }
    }


    // Create triangle mesh
    DfTriangleMesh triangleMesh;
    if (not triangleMesh.Initialize(quadPoints.size(), triVerts, &errMsg)) {
        
        TF_WARN(errMsg);
        delete vertexBuffer;
        delete vvBuffer;                            
        return false;
    }

    triangleMesh.UpdateGeometry(quadPoints);
    hbv.Update(triangleMesh.ComputeFaceRanges());


    createTriangleMeshTime.Stop();    


    std::vector<int> coords(pointsToProject.size());
    closestPointProjectTime.Start();
    int numTris = triVerts.size()/3;
    for (unsigned int i = 0; i < pointsToProject.size(); ++i) {
        GfVec3d hit;
        size_t faceIndex;
        triangleMesh.FindClosestFace(
            hbv, pointsToProject[i], &faceIndex, &hit);
        
        if (faceIndex >= numTris) {
            TF_WARN("Bad indexing in Projector");
            delete vertexBuffer;
            delete vvBuffer;                            
            return false;            
        }

        // Create a Gf triangle from the three
        // points on the triangle hit
        int vertIndex = faceIndex*3;
        GfTriangle gfTriangle(
            quadPoints[triVerts[vertIndex]],
            quadPoints[triVerts[vertIndex+1]],
            quadPoints[triVerts[vertIndex+2]]);;


        double bn;
        GfVec3d baryCentric = gfTriangle.GetBarycentric(hit, &bn);

        // The closest point found on the triangle is recorded with
        // the indices of the triangle points and the floating point
        // coefficients that define the point's parametric location
        // inside the triangle
        pointsProjectedOntoSubdiv->push_back(
            WeightObject_Projector::EvalCoords (
                triVerts[vertIndex],
                triVerts[vertIndex+1],
                triVerts[vertIndex+2],
                baryCentric[0], baryCentric[1], baryCentric[2]));
    }
    
    closestPointProjectTime.Stop();

/*    
    printf("Hbv results for %d point projections:\n\ttotal #pointsInTree=%d \n\tupdate tessellator time = %fs,  \n\tcreateTriangleMeshTime = %fs,  \n\tclosestPoint project time = %fs,\n",
           (int)pointsToProject.size(), 
           (int)coarsePoints.size(),
           updateTessellatorTime.GetSeconds(),           
           createTriangleMeshTime.GetSeconds(),
           closestPointProjectTime.GetSeconds());
*/

    delete vertexBuffer;        
    delete vvBuffer;
    
    return true;
}





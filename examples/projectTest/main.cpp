//
//     Copyright (C) Pixar. All rights reserved.
//
//     This license governs use of the accompanying software. If you
//     use the software, you accept this license. If you do not accept
//     the license, do not use the software.
//
//     1. Definitions
//     The terms "reproduce," "reproduction," "derivative works," and
//     "distribution" have the same meaning here as under U.S.
//     copyright law.  A "contribution" is the original software, or
//     any additions or changes to the software.
//     A "contributor" is any person or entity that distributes its
//     contribution under this license.
//     "Licensed patents" are a contributor's patent claims that read
//     directly on its contribution.
//
//     2. Grant of Rights
//     (A) Copyright Grant- Subject to the terms of this license,
//     including the license conditions and limitations in section 3,
//     each contributor grants you a non-exclusive, worldwide,
//     royalty-free copyright license to reproduce its contribution,
//     prepare derivative works of its contribution, and distribute
//     its contribution or any derivative works that you create.
//     (B) Patent Grant- Subject to the terms of this license,
//     including the license conditions and limitations in section 3,
//     each contributor grants you a non-exclusive, worldwide,
//     royalty-free license under its licensed patents to make, have
//     made, use, sell, offer for sale, import, and/or otherwise
//     dispose of its contribution in the software or derivative works
//     of the contribution in the software.
//
//     3. Conditions and Limitations
//     (A) No Trademark License- This license does not grant you
//     rights to use any contributor's name, logo, or trademarks.
//     (B) If you bring a patent claim against any contributor over
//     patents that you claim are infringed by the software, your
//     patent license from such contributor to the software ends
//     automatically.
//     (C) If you distribute any portion of the software, you must
//     retain all copyright, patent, trademark, and attribution
//     notices that are present in the software.
//     (D) If you distribute any portion of the software in source
//     code form, you may do so only under this license by including a
//     complete copy of this license with your distribution. If you
//     distribute any portion of the software in compiled or object
//     code form, you may only do so under a license that complies
//     with this license.
//     (E) The software is licensed "as-is." You bear the risk of
//     using it. The contributors give no express warranties,
//     guarantees or conditions. You may have additional consumer
//     rights under your local laws which this license cannot change.
//     To the extent permitted under your local laws, the contributors
//     exclude the implied warranties of merchantability, fitness for
//     a particular purpose and non-infringement.
//

    #include <stdlib.h>

#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>
#include <osd/cpuEvalLimitContext.h>
#include <osd/cpuEvalLimitController.h>
#include <osd/cpuVertexBuffer.h>
#include <osd/error.h>
#include <osd/mesh.h>
#include <osd/vertex.h>

#include "shape.h"

#include "../common/stopwatch.h"

#include <cfloat>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#ifdef OPENSUBDIV_HAS_OPENMP
    #include <omp.h>
#endif

using namespace OpenSubdiv;

//------------------------------------------------------------------------------
typedef HbrMesh<OsdVertex>     OsdHbrMesh;
typedef HbrVertex<OsdVertex>   OsdHbrVertex;
typedef HbrFace<OsdVertex>     OsdHbrFace;
typedef HbrHalfedge<OsdVertex> OsdHbrHalfedge;

typedef FarMesh<OsdVertex>              OsdFarMesh;
typedef FarMeshFactory<OsdVertex>       OsdFarMeshFactory;
typedef FarSubdivisionTables<OsdVertex> OsdFarMeshSubdivision;


#ifdef NOT_NEEDED
//------------------------------------------------------------------------------
static int
GetNumPtexFaces( OsdHbrMesh const * hmesh, int nfaces ) {

    OsdHbrFace * lastface = hmesh->GetFace( nfaces-1 );
    assert(lastface);
    
    int result = lastface->GetPtexIndex();
    
    result += (hmesh->GetSubdivision()->FaceIsExtraordinary(hmesh, lastface) ? 
                  lastface->GetNumVertices() : 0);

    return ++result;
}
#endif


//------------------------------------------------------------------------------
static void
createOsdMesh(int level)
{
    float points[] = { 0.000000, -1.414214, 1.000000,
                       1.414214, 0.000000, 1.000000,    
                       -1.414214, 0.000000, 1.000000,
                       0.000000, 1.414214, 1.000000,
                       -1.414214, 0.000000, -1.000000,
                       0.000000, 1.414214, -1.000000,
                       0.000000, -1.414214, -1.000000,
                       1.414214, 0.000000, -1.000000};

    int nverts[] = { 4, 4, 4, 4, 4, 4};

    int indices[] = { 0, 1, 3, 2,
                      2, 3, 5, 4,
                      4, 5, 7, 6,
                      6, 7, 1, 0,
                      1, 7, 5, 3,
                      6, 0, 2, 4};
    

//    Scheme scheme = kCatmark;

    int numVertices = sizeof(points)/(3*sizeof(float));
    OpenSubdivShape *shape = new OpenSubdivShape();
    std::string errorMessage;
    if (not shape->Initialize(numVertices, 
                              nverts, sizeof(nverts)/sizeof(int),
                              indices, sizeof(indices)/sizeof(int), 2,
                              &errorMessage)) {
        std::cout  << "Shape initialization failed with " << errorMessage << std::endl;
    }


    // Push the vertex data:
    std::vector<float> pointsVec;
    pointsVec.resize(sizeof(points));
    for (int i=0; i<(int)sizeof(points); ++i) {
        pointsVec[i] = points[i];
        
    }
    shape->SetCoarsePositions(pointsVec);

    std::vector<float> refinedPositions;
    std::vector<int> refinedQuads;

    if (not (shape->Refine(2)                                      and
             shape->GetPositions(&refinedPositions, &errorMessage) and
             shape->GetQuads(&refinedQuads, &errorMessage))) {
        std::cout << errorMessage << std::endl;
    } else {
        std::cout << "Hot damn, it worked.\n";
        for (int i=0; i<(int)refinedPositions.size(); i+=3)  {
            std::cout << "(" << refinedPositions[i] <<
                ", " << refinedPositions[i+1] <<
                "," << refinedPositions[i+2] << ")\n";
        }
        for (int i=0; i<(int)refinedQuads.size(); i+=4)  {
            std::cout << "(" << refinedQuads[i] <<
                ", " << refinedQuads[i+1] <<
                "," << refinedQuads[i+2] << ")\n";
        }
        
    }
         
    
        
    
//    coarseMesh->SetVVData(vvInterleaved);

    
}

//------------------------------------------------------------------------------
static void
callbackError(OpenSubdiv::OsdErrorType err, const char *message)
{
    printf("OsdError: %d\n", err);
    printf("%s", message);
}


//------------------------------------------------------------------------------
int main(int, char**) {


    OsdSetErrorCallback(callbackError);

    createOsdMesh(1);

}

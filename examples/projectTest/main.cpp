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

    int nvertsPerFace[] = { 4, 4, 4, 4, 4, 4};

    int faceverts[] = { 1, 2, 4, 3,
                        3, 4, 6, 5,
                        5, 6, 8, 7,
                        7, 8, 2, 1,
                        2, 8, 6, 4,
                        7, 1, 3, 5};
    

//    Scheme scheme = kCatmark;

    
    OpenSubdivShape shape(points, sizeof(points)/(3*sizeof(float)),
                          nvertsPerFace, sizeof(nvertsPerFace)/sizeof(int),
                          faceverts, sizeof(faceverts)/sizeof(int));

#ifdef NOT_YET    
    OsdHbrMesh *hmesh = shape.GetHbrMesh();
    
    std::cout << "Made the shape\n";

    g_positions.resize(g_orgPositions.size(),0.0f);

    int nfaces = hmesh->GetNumFaces(),
        nptexfaces = GetNumPtexFaces(hmesh, nfaces);
    
    // Create FAR mesh
    OsdFarMeshFactory factory( hmesh, level, /*adaptive*/ true);    
    
    delete g_fmesh;
    g_fmesh = factory.Create( /*ptex*/ true, /*fvar*/ false);
    
    int nverts = g_fmesh->GetNumVertices();
    

    
    // Create v-buffer & populate w/ positions
    delete g_samplesVB;
    g_samplesVB = OsdCpuVertexBuffer::Create(3, nverts);

        
    // Create a Compute context, used to "pose" the vertices
    delete g_computeCtx;
    g_computeCtx = OsdCpuComputeContext::Create(g_fmesh);
    
    g_computeCtrl.Refine( g_computeCtx, g_fmesh->GetKernelBatches(), g_samplesVB );

    
    // Create eval context & data buffers
    delete g_evalCtx;
    g_evalCtx = OsdCpuEvalLimitContext::Create(g_fmesh);

    delete g_Q;
    g_Q = OsdCpuGLVertexBuffer::Create(6,nsamples);

    delete g_dQu;
    g_dQu = OsdCpuGLVertexBuffer::Create(6,nsamples);

    delete g_dQv;
    g_dQv = OsdCpuGLVertexBuffer::Create(6,nsamples);
        
    updateGeom();

    // Bind g_Q as a GL_POINTS VBO
    glBindVertexArray(g_samplesVAO);
    
    glBindBuffer(GL_ARRAY_BUFFER, g_Q->BindVBO());

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof (GLfloat) * 6, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof (GLfloat) * 6, (float*)12);

    glBindVertexArray(0);
#endif
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

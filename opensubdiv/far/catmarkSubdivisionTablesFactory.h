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

#ifndef FAR_CATMARK_SUBDIVISION_TABLES_FACTORY_H
#define FAR_CATMARK_SUBDIVISION_TABLES_FACTORY_H

#include <cassert>
#include <vector>

#include "../version.h"

#include "../far/catmarkSubdivisionTables.h"
#include "../far/meshFactory.h"
#include "../far/subdivisionTablesFactory.h"

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

template <class T, class U> class FarMeshFactory;

/// \brief A specialized factory for FarCatmarkSubdivisionTables
///
/// Separating the factory allows us to isolate Far data structures from Hbr dependencies.
///
template <class T, class U> class FarCatmarkSubdivisionTablesFactory {
protected:
    template <class X, class Y> friend class FarMeshFactory;

    /// Creates a FarCatmarkSubdivisiontables instance.
    static FarCatmarkSubdivisionTables<U> * Create( FarMeshFactory<T,U> * meshFactory, FarMesh<U> * farMesh, FarKernelBatchVector *batches  );
};

// This factory walks the Hbr vertices and accumulates the weights and adjacency
// (valance) information specific to the catmark subdivision scheme. The results
// are stored in a FarCatmarkSubdivisionTable<U>.
template <class T, class U> FarCatmarkSubdivisionTables<U> * 
FarCatmarkSubdivisionTablesFactory<T,U>::Create( FarMeshFactory<T,U> * meshFactory, FarMesh<U> * farMesh, FarKernelBatchVector *batches ) {

    assert( meshFactory and farMesh );

    int maxlevel = meshFactory->GetMaxLevel();
    
    std::vector<int> & remap = meshFactory->getRemappingTable();
    
    FarSubdivisionTablesFactory<T,U> tablesFactory( meshFactory->GetHbrMesh(), maxlevel, remap );

    FarCatmarkSubdivisionTables<U> * result = new FarCatmarkSubdivisionTables<U>(farMesh, maxlevel);

    // Allocate memory for the indexing tables
    result->_F_ITa.resize(tablesFactory.GetNumFaceVerticesTotal(maxlevel)*2);
    result->_F_IT.resize(tablesFactory.GetFaceVertsValenceSum());

    result->_E_IT.resize(tablesFactory.GetNumEdgeVerticesTotal(maxlevel)*4);
    result->_E_W.resize(tablesFactory.GetNumEdgeVerticesTotal(maxlevel)*2);

    result->_V_ITa.resize((tablesFactory.GetNumVertexVerticesTotal(maxlevel)
                           - tablesFactory.GetNumVertexVerticesTotal(0))*5); // subtract corase cage vertices
    result->_V_IT.resize(tablesFactory.GetVertVertsValenceSum()*2);
    result->_V_W.resize(tablesFactory.GetNumVertexVerticesTotal(maxlevel)
                        - tablesFactory.GetNumVertexVerticesTotal(0));

    // Prepare batch table
    batches->reserve(maxlevel*5);

    int vertexOffset = 0;
    int F_IT_offset = 0;
    int V_IT_offset = 0;
    int faceTableOffset = 0;
    int edgeTableOffset = 0;
    int vertTableOffset = 0;

    unsigned int * F_IT = result->_F_IT.empty() ? 0 : &result->_F_IT[0];
    int * F_ITa = result->_F_ITa.empty() ? 0 : &result->_F_ITa[0];
    int * E_IT = result->_E_IT.empty() ? 0 : &result->_E_IT[0];
    float * E_W = result->_E_W.empty() ? 0 : &result->_E_W[0];
    int * V_ITa = result->_V_ITa.empty() ? 0 : &result->_V_ITa[0];
    unsigned int * V_IT = result->_V_IT.empty() ? 0 : &result->_V_IT[0];
    float * V_W = result->_V_W.empty() ? 0 : &result->_V_W[0];

    for (int level=1; level<=maxlevel; ++level) {

        // pointer to the first vertex corresponding to this level
        vertexOffset = tablesFactory._vertVertIdx[level-1] +
            (int)tablesFactory._vertVertsList[level-1].size();
        result->_vertsOffsets[level] = vertexOffset;

        // Face vertices
        // "For each vertex, gather all the vertices from the parent face."
        int nFaceVertices = (int)tablesFactory._faceVertsList[level].size();

        // add a batch for face vertices
        if (nFaceVertices > 0)  // in torus case, nfacevertices could be zero
            batches->push_back(FarKernelBatch(level,
                                              CATMARK_FACE_VERTEX,
                                              0,
                                              0,
                                              nFaceVertices,
                                              faceTableOffset,
                                              vertexOffset));

        vertexOffset += nFaceVertices;
        faceTableOffset += nFaceVertices;

        for (int i=0; i < nFaceVertices; ++i) {

            HbrVertex<T> * v = tablesFactory._faceVertsList[level][i];
            assert(v);

            HbrFace<T> * f=v->GetParentFace();
            assert(f);

            int valence = f->GetNumVertices();

            F_ITa[2*i+0] = F_IT_offset;
            F_ITa[2*i+1] = valence;

            for (int j=0; j<valence; ++j)
                F_IT[F_IT_offset++] = remap[f->GetVertex(j)->GetID()];
        }
        F_ITa += nFaceVertices * 2;

        // Edge vertices

        // Triangular interpolation mode :
        // see "smoothtriangle" tag introduced in prman 3.9 and HbrCatmarkSubdivision<T>
        typename HbrCatmarkSubdivision<T>::TriangleSubdivision triangleMethod =
            dynamic_cast<HbrCatmarkSubdivision<T> *>(meshFactory->GetHbrMesh()->GetSubdivision())->GetTriangleSubdivisionMethod();

        // "For each vertex, gather the 2 vertices from the parent edege and the
        // 2 child vertices from the faces to the left and right of that edge.
        // Adjust if edge has a crease or is on a boundary."
        int nEdgeVertices = (int)tablesFactory._edgeVertsList[level].size();

        // add a batch for edge vertices
        if (nEdgeVertices > 0)
            batches->push_back(FarKernelBatch(level,
                                              CATMARK_EDGE_VERTEX,
                                              0,
                                              0,
                                              nEdgeVertices,
                                              edgeTableOffset,
                                              vertexOffset));

        vertexOffset += nEdgeVertices;
        edgeTableOffset += nEdgeVertices;

        for (int i=0; i < nEdgeVertices; ++i) {

            HbrVertex<T> * v = tablesFactory._edgeVertsList[level][i];
            assert(v);
            HbrHalfedge<T> * e = v->GetParentEdge();
            assert(e);

            float esharp = e->GetSharpness();

            // get the indices 2 vertices from the parent edge
            E_IT[4*i+0] = remap[e->GetOrgVertex()->GetID()];
            E_IT[4*i+1] = remap[e->GetDestVertex()->GetID()];

            float faceWeight=0.5f, vertWeight=0.5f;

            // in the case of a fractional sharpness, set the adjacent faces vertices
            if (!e->IsBoundary() && esharp <= 1.0f) {

                float leftWeight, rightWeight;
                HbrFace<T>* rf = e->GetRightFace();
                HbrFace<T>* lf = e->GetLeftFace();

                leftWeight = ( triangleMethod == HbrCatmarkSubdivision<T>::k_New && lf->GetNumVertices() == 3) ? HBR_SMOOTH_TRI_EDGE_WEIGHT : 0.25f;
                rightWeight = ( triangleMethod == HbrCatmarkSubdivision<T>::k_New && rf->GetNumVertices() == 3) ? HBR_SMOOTH_TRI_EDGE_WEIGHT : 0.25f;

                faceWeight = 0.5f * (leftWeight + rightWeight);
                vertWeight = 0.5f * (1.0f - 2.0f * faceWeight);

                faceWeight *= (1.0f - esharp);

                vertWeight = 0.5f * esharp + (1.0f - esharp) * vertWeight;

                E_IT[4*i+2] = remap[lf->Subdivide()->GetID()];
                E_IT[4*i+3] = remap[rf->Subdivide()->GetID()];
            } else {
                E_IT[4*i+2] = -1;
                E_IT[4*i+3] = -1;
            }
            E_W[2*i+0] = vertWeight;
            E_W[2*i+1] = faceWeight;
        }
        E_IT += 4 * nEdgeVertices;
        E_W += 2 * nEdgeVertices;

        // Vertex vertices

        FarVertexKernelBatchFactory batchFactory((int)tablesFactory._vertVertsList[level].size(), 0);

        int nVertVertices = (int)tablesFactory._vertVertsList[level].size();
        for (int i=0; i < nVertVertices; ++i) {

            HbrVertex<T> * v = tablesFactory._vertVertsList[level][i],
                         * pv = v->GetParentVertex();
            assert(v and pv);

            // Look at HbrCatmarkSubdivision<T>::Subdivide for more details about
            // the multi-pass interpolation
            unsigned char masks[2];
            int npasses;
            float weights[2];
            masks[0] = pv->GetMask(false);
            masks[1] = pv->GetMask(true);

            // If the masks are identical, only a single pass is necessary. If the
            // vertex is transitioning to another rule, two passes are necessary,
            // except when transitioning from k_Dart to k_Smooth : the same
            // compute kernel is applied twice. Combining this special case allows
            // to batch the compute kernels into fewer calls.
            if (masks[0] != masks[1] and (
                not (masks[0]==HbrVertex<T>::k_Smooth and
                     masks[1]==HbrVertex<T>::k_Dart))) {
                weights[1] = pv->GetFractionalMask();
                weights[0] = 1.0f - weights[1];
                npasses = 2;
            } else {
                weights[0] = 1.0f;
                weights[1] = 0.0f;
                npasses = 1;
            }

            int rank = FarSubdivisionTablesFactory<T,U>::GetMaskRanking(masks[0], masks[1]);

            V_ITa[5*i+0] = V_IT_offset;
            V_ITa[5*i+1] = 0;
            V_ITa[5*i+2] = remap[ pv->GetID() ];
            V_ITa[5*i+3] = -1;
            V_ITa[5*i+4] = -1;

            for (int p=0; p<npasses; ++p)
                switch (masks[p]) {
                    case HbrVertex<T>::k_Smooth :
                    case HbrVertex<T>::k_Dart : {
                        HbrHalfedge<T> *e = pv->GetIncidentEdge(),
                                       *start = e;
                        while (e) {
                            V_ITa[5*i+1]++;

                            V_IT[V_IT_offset++] = remap[ e->GetDestVertex()->GetID() ];

                            V_IT[V_IT_offset++] = remap[ e->GetLeftFace()->Subdivide()->GetID() ];

                            e = e->GetPrev()->GetOpposite();

                            if (e==start) break;
                        }
                        break;
                    }
                    case HbrVertex<T>::k_Crease : {

                        class GatherCreaseEdgesOperator : public HbrHalfedgeOperator<T> {
                        public:
                            HbrVertex<T> * vertex; int eidx[2]; int count; bool next;

                            GatherCreaseEdgesOperator(HbrVertex<T> * v, bool n) : vertex(v), count(0), next(n) { eidx[0]=-1; eidx[1]=-1; }

                            virtual void operator() (HbrHalfedge<T> &e) {
                                if (e.IsSharp(next) and count < 2) {
                                    HbrVertex<T> * a = e.GetDestVertex();
                                    if (a==vertex)
                                        a = e.GetOrgVertex();
                                    eidx[count++]=a->GetID();
                                }
                            }
                        };

                        GatherCreaseEdgesOperator op( pv, p==1 );
                        pv->ApplyOperatorSurroundingEdges( op );

                        assert(V_ITa[5*i+3]==-1 and V_ITa[5*i+4]==-1);
                        assert(op.eidx[0]!=-1 and op.eidx[1]!=-1);
                        V_ITa[5*i+3] = remap[op.eidx[0]];
                        V_ITa[5*i+4] = remap[op.eidx[1]];
                        break;
                    }
                    case HbrVertex<T>::k_Corner :
                        // in the case of a k_Crease / k_Corner pass combination, we
                        // need to set the valence to -1 to tell the "B" Kernel to
                        // switch to k_Corner rule (as edge indices won't be -1)
                        if (V_ITa[5*i+1]==0)
                            V_ITa[5*i+1] = -1;

                    default : break;
                }

            if (rank>7)
                // the k_Corner and k_Crease single-pass cases apply a weight of 1.0
                // but this value is inverted in the kernel
                V_W[i] = 0.0;
            else
                V_W[i] = weights[0];

            batchFactory.AddVertex( i, rank );
        }
        V_ITa += nVertVertices*5;
        V_W += nVertVertices;

        // add batches for vert vertices
        if (nVertVertices > 0)
            batchFactory.AppendCatmarkBatches(batches, level, vertTableOffset, vertexOffset);
        vertexOffset += nVertVertices;
        vertTableOffset += nVertVertices;
    }
    result->_vertsOffsets[maxlevel+1] = vertexOffset;
    return result;
}

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;

} // end namespace OpenSubdiv

#endif /* FAR_CATMARK_SUBDIVISION_TABLES_FACTORY_H */

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


#include "shape.h"

#include <hbr/bilinear.h>
#include <hbr/loop.h>
#include <hbr/catmark.h>
#include <hbr/vertexEdit.h>
#include <hbr/cornerEdit.h>
#include <hbr/holeEdit.h>

using namespace OpenSubdiv;
using namespace std;

//------------------------------------------------------------------------------
template <class T> std::string   
hbrToObj( HbrMesh<T> * mesh ) {

    std::stringstream sh;

    sh<<"# This file uses centimeters as units for non-parametric coordinates.\n\n";
    
    int nv = mesh->GetNumVertices();
    for (int i=0; i<nv; ++i) {
       const float * pos = mesh->GetVertex(i)->GetData().GetPos();
       sh << "v " << pos[0] << " " << pos[1] << " " << pos[2] <<"\n";
    }
    
    int nf = mesh->GetNumFaces();
    for (int i=0; i<nf; ++i) {
    
        sh << "f ";
        
        HbrFace<T> * f = mesh->GetFace(i);
        
        for (int j=0; j<f->GetNumVertices(); ++j) {
            int vert = f->GetVertex(j)->GetID()+1;
            sh << vert << "/" << vert << "/" << vert << " ";
        }
        sh << "\n";
    }

    sh << "\n";

    return sh.str();
}

//------------------------------------------------------------------------------
template <class T> void
createVertices( OpenSubdivShape const * sh, HbrMesh<T> * mesh, std::vector<float> * verts ) {

    T v;
    for(int i=0;i<sh->GetNverts(); i++ ) {
        v.SetPosition( sh->GetVerts()[i*3], sh->GetVerts()[i*3+1], sh->GetVerts()[i*3+2] );
        mesh->NewVertex( i, v );
    }
}

//------------------------------------------------------------------------------
template <class T> void
createVertices( OpenSubdivShape const * sh, HbrMesh<T> * mesh, std::vector<float> & verts ) {

    T v;
    for(int i=0;i<sh->GetNverts(); i++ )
        mesh->NewVertex( i, v );
}
#ifdef NOT_YET
//------------------------------------------------------------------------------
template <class T> void
copyVertexPositions( OpenSubdivShape const * sh, HbrMesh<T> * mesh, std::vector<float> & verts ) {

    int nverts = mesh->GetNumVertices();
    
    verts.resize( nverts * 3 );
    
    std::copy(sh->GetVerts.begin(), sh->GetVerts.end(), verts.begin());
    
    // Sometimes Hbr dupes some vertices during Mesh::Finish()
    if (nverts > sh->GetNverts()) {
    
        for (int i=sh->GetNverts(); i<nverts; ++i) {
        
            HbrVertex<T> * v = mesh->GetVertex(i);
            
            HbrFace<T> * f = v->GetIncidentEdge()->GetFace();
            
            int vidx = -1;
            for (int j=0; j<f->GetNumVertices(); ++j)
                if (f->GetVertex(j)==v) {
                    vidx = j;
                    break;
                }
            assert(vidx>-1);
        
            const int * shfaces = &sh->faceverts[0];
            for (int j=0; j<f->GetID(); ++j)
                shfaces += sh->nvertsPerFace[j];
        
            int shvert = shfaces[vidx];
            
            verts[i*3+0] = sh->GetVerts[shvert*3+0];
            verts[i*3+1] = sh->GetVerts[shvert*3+1];
            verts[i*3+2] = sh->GetVerts[shvert*3+2];
        }
    }
}


//------------------------------------------------------------------------------
template <class T> void
createTopology( OpenSubdivShape const * sh, HbrMesh<T> * mesh, Scheme scheme) {
}

//------------------------------------------------------------------------------
template <class T> HbrMesh<T> *
simpleHbr(char const * shapestr, Scheme scheme, std::vector<float> * verts=0) {

    OpenSubdivShape * sh = OpenSubdivShape::parseShape( shapestr );

    HbrMesh<T> * mesh = createMesh<T>(scheme);

    createVertices<T>(sh, mesh, verts);

    createTopology<T>(sh, mesh, scheme);

    if(verts)
        copyVertexPositions<T>(sh,mesh,*verts);

    delete sh;

    return mesh;
}

//------------------------------------------------------------------------------
template <class T> HbrMesh<T> *
simpleHbr(char const * shapestr, Scheme scheme, std::vector<float> & verts) {

    OpenSubdivShape * sh = OpenSubdivShape::parseShape( shapestr );

    HbrMesh<T> * mesh = createMesh<T>(scheme);

    createVertices<T>(sh, mesh, verts);

    createTopology<T>(sh, mesh, scheme);

    copyVertexPositions<T>(sh,mesh,verts);

    delete sh;

    return mesh;
}
#endif



//------------------------------------------------------------------------------
// The simplest constructor, only point positions and polygonal
// mesh topology
/*
OpenSubdivShape::OpenSubdivShape(const std::vector<float>  &verts,
                                 const std::vector<int>    &nvertsPerFace,
                                 const std::vector<int>    &faceverts) :
    _verts(verts),
    _nvertsPerFace(nvertsPerFace),
    _faceverts(faceverts),
    _hbrMesh(NULL)
{

}
*/

OpenSubdivShape::OpenSubdivShape(
    string name,
    int numVertices,
    int maxLevels,
    vector<int> indices,
    vector<int> nverts,
    vector<string> vvNames,
    vector<string> fvNames,
    const vector<float>& fvData,
    tag& tagData,
    Scheme scheme)
{
    _scheme = scheme;
    
    HbrMesh<OsdVertex> * mesh = 0;
    
    static HbrBilinearSubdivision<OsdVertex> _bilinear;
    static HbrLoopSubdivision<OsdVertex>     _loop;
    static HbrCatmarkSubdivision<OsdVertex>  _catmark;
    
    switch (scheme) {
    case kBilinear : mesh = new HbrMesh<OsdVertex>( &_bilinear ); break;
    case kLoop     : mesh = new HbrMesh<OsdVertex>( &_loop     ); break;
    case kCatmark  : mesh = new HbrMesh<OsdVertex>( &_catmark  ); break;
    }

    
    //_name = name;
    if (fvNames.empty()) {
        _hbrMesh = new HbrMesh<OsdVertex>(&_catmark);
    } else {
        int fvarcount = (int) fvNames.size();

        // For now we only handle 1 float per FV variable.
        _fvarwidths.assign(fvarcount, 1);

        int startIndex = 0;
        for (int fvarindex = 0; fvarindex < fvarcount; ++fvarindex) {
            _fvarindices.push_back(startIndex);
            _fvaroffsets[fvNames[fvarindex]] = startIndex;
            startIndex += _fvarwidths[fvarindex];
        }

        _hbrMesh = new HbrMesh<OsdVertex>(&_catmark, fvarcount, &_fvarindices[0],
                               &_fvarwidths[0], fvarcount);
    }
    _farMesh = 0;
    _computeContext = 0;
    _vvNames = vvNames;
    _fvNames = fvNames;

    OsdVertex v;
    for (int i = 0; i < numVertices; ++i) {
        HbrVertex* hvert = _hbrMesh->NewVertex(i, v);
        if (!hvert) {
            TF_FATAL_ERROR("Unable to create OSD vertex from tidscene");
        }
    }

    // Sanity check
    int fvarWidth = _hbrMesh->GetTotalFVarWidth();
    if (fvData.size() < nverts.size() * fvarWidth ||
        fvarWidth != fvNames.size()) {
        TF_FATAL_ERROR("Incorrectly sized face data: name count = %d, "
                       "data width = %d, face count = %d, total data size = %d.",
                       (int) fvNames.size(),
                       fvarWidth,
                       (int) nverts.size(),
                       (int) fvData.size());
    }

    // Create faces
    const int * fv=&(_faceverts[0]);
    // face-vertex count offset
    int fvcOffset = 0;    
    for(int f=0, ptxidx=0;f< GetNfaces(); f++ ) {

        int nv = _nvertsPerFace[f];

        if ((_scheme==kLoop) and (nv!=3)) {
            printf("Trying to create a Loop surbd with non-triangle face\n");
            exit(1);
        }

        for(int j=0;j<nv;j++) {
            HbrVertex<OsdVertex> * origin      = _hbrMesh->GetVertex( fv[j] );
            HbrVertex<OsdVertex> * destination = _hbrMesh->GetVertex( fv[ (j+1)%nv] );
            HbrHalfedge<OsdVertex> * opposite  = destination->GetEdge(origin);

            if(origin==NULL || destination==NULL) {
                printf(" An edge was specified that connected a nonexistent vertex\n");
                exit(1);
            }

            if(origin == destination) {
                printf(" An edge was specified that connected a vertex to itself\n");
                exit(1);
            }

            if(opposite && opposite->GetOpposite() ) {
                printf(" A non-manifold edge incident to more than 2 faces was found\n");
                exit(1);
            }

            if(origin->GetEdge(destination)) {
                printf(" An edge connecting two vertices was specified more than once."
                       " It's likely that an incident face was flipped\n");
                exit(1);
            }
        }

        HbrFace<OsdVertex> * face = _hbrMesh->NewFace(nv, (int *)fv, 0);

        face->SetPtexIndex(ptxidx);


        if (!hface) {
            TF_FATAL_ERROR("Unable to create OSD face from tidscene");
        }

        // prideout: 3/21/2013 - Inspired by "GetFVarData" in examples/mayaViewer/hbrUtil.cpp
        if (!fvNames.empty()) {
            const float* faceData = &(fvData[fvcOffset*fvarWidth]);
            for (int fvi = 0; fvi < nv; ++fvi) {
                int vindex = indices[fvi + fvcOffset];
                HbrVertex* v = _hbrMesh->GetVertex(vindex);
                HbrFVarData& fvarData = v->GetFVarData(hface);
                if (!fvarData.IsInitialized()) {
                    fvarData.SetAllData(fvarWidth, faceData); 
                } else if (!fvarData.CompareAll(fvarWidth, faceData)) {

                    // If data exists for this face vertex, but is different
                    // (e.g. we're on a UV seam) create another fvar datum
                    HbrFVarData& fvarData = v->NewFVarData(hface);
                    fvarData.SetAllData(fvarWidth, faceData);
                }

                // Advance pointer to next set of face-varying data
                faceData += fvarWidth;
            }
        }

        fvcOffset += nv;
        
        if ( (_scheme==kCatmark or _scheme==kBilinear) and nv != 4 )
            ptxidx+=nv;
        else
            ptxidx++;

        fv+=nv;
    }


    _ProcessTagsAndFinishMesh(
        _hbrMesh, tagData.tags, tagData.numArgs, tagData.intArgs,
        tagData.floatArgs, tagData.stringArgs);

    FarMeshFactory meshFactory(_hbrMesh, maxLevels, false);

    _farMesh = meshFactory.Create();
    _computeContext = ComputeContext::Create(_farMesh);
    _vertexBuffer = VertexBuffer::Create(
        3, _farMesh->GetNumVertices());

    if (vvNames.empty()) {
        _vvBuffer = 0;
    } else {
        _vvBuffer = VertexBuffer::Create(
            vvNames.size(), _farMesh->GetNumVertices());
    }
}

void
OpenSubdivShape::SetCoarsePositions(const vector<float>& coords)
{
    const float* pFloats = &coords.front();
    int numFloats = (int) coords.size();
    _vertexBuffer->UpdateData(pFloats, 0, numFloats / 3);
}

void
OpenSubdivShape::SetVVData(const vector<float>& data)
{
    if (!_vvBuffer) {
        if (!data.empty()) {
            TF_FATAL_ERROR("Mesh was not constructed with VV variables.");
        }
        return;
    }

    int numElements = _vvBuffer->GetNumElements();
    int numVertices = (int) data.size() / numElements;
    const float* pFloats = &data.front();
    _vvBuffer->UpdateData(pFloats, 0, numVertices);
}

// Interleave the face-varying sets specified by "names", adding
// floats into the "fvdata" vector.  The number of added floats is:
//    names.size() * NumRefinedFaces * 4
void
OpenSubdivShape::GetFVData(
    int level, const vector<string>& names, vector<float>* outdata)
{
    // First some sanity checking.
    if (!outdata) {
        return;
    }
    TF_FOR_ALL(nameit, names) {
        string name = *nameit;
        if (_fvaroffsets.find(name) == _fvaroffsets.end()) {
            printf("Can't find facevarying variable %s\n", name.c_str());
            return;
        }
    }

    // Fetch *all* faces; this includes all subdivision levels.
    vector<HbrFace *> faces;
    _hbrMesh->GetFaces(std::back_inserter(faces));

    // Iterate through all faces, filtering on the requested subdivision level.
    TF_FOR_ALL(faceit, faces) {
        HbrFace* face = *faceit;
        if (face->GetDepth() != level) {
            continue;
        }
        int ncorners = face->GetNumVertices();
        for (int corner = 0; corner < ncorners; ++corner) {
            HbrFVarData& fvariable = face->GetFVarData(corner);
            TF_FOR_ALL(nameit, names) {
                string name = *nameit;
                int offset = _fvaroffsets[name];
                const float* data = fvariable.GetData(offset);
                outdata->push_back(*data);
            }
        }
    }
}


OpenSubdivShape::OpenSubdivShape(const float* points,
                                 int pointsLen /*= # of 3-float points*/, 
                                 const int *nvertsPerFace, int numFaces,
                                 const int *faceverts, int facevertsLen)
{
    _verts.resize(pointsLen*3);
    for (int i=0; i<pointsLen*3; ++i) {
        _verts[i]= points[i];
    }
    
    _nvertsPerFace.resize(numFaces);
    int totalNumFaceVerts=0;
    for (int i=0; i<numFaces; ++i) {
        _nvertsPerFace[i] = nvertsPerFace[i];
        totalNumFaceVerts += nvertsPerFace[i];
    }

    if (totalNumFaceVerts != facevertsLen) {
        // XXX real error needed
        printf("Bad mesh passed to OpenSubdivShape\n");
        _verts.clear();
        _nvertsPerFace.clear();
        return;
    }

    _faceverts.resize(totalNumFaceVerts);
    int badIndices=0;
    for (int i=0; i<facevertsLen; ++i) {
        if (faceverts[i] >= pointsLen) {
            badIndices++;
            _faceverts[i] = 0;
        } else {
            _faceverts[i] = faceverts[i];
        }
    }
}
    

//------------------------------------------------------------------------------
OpenSubdivShape::~OpenSubdivShape() {
    for (int i=0; i<(int)_tags.size(); ++i)
        delete _tags[i];

    delete _hbrMesh;
    delete _farMesh;
    delete _computeContext;
    delete _vertexBuffer;
    delete _vvBuffer;
}


HbrMesh<OsdVertex>*
OpenSubdivShape::GetHbrMesh()
{

    if (_hbrMesh) {
        return _hbrMesh;
    }

    _hbrMesh = CreateMesh();
    
    
    return _hbrMesh;
}


FarMesh<OsdVertex>*
OpenSubdivShape::GetFarMesh()
{

    if (_farMesh) {
        return _farMesh;
    }

    // XXX how to pass in level and adaptive?
    int level = 2;
    
    // Create FAR mesh
    OsdFarMeshFactory factory( GetHbrMesh(), level, /*adaptive*/ false);    
    
    _farMesh = factory.Create( /*ptex*/ true, /*fvar*/ false);

}

//------------------------------------------------------------------------------
HbrMesh<OsdVertex>*
OpenSubdivShape::CreateMesh(Scheme scheme)
{       
    return mesh;
}



//------------------------------------------------------------------------------
std::string OpenSubdivShape::genRIB() const {
    std::stringstream rib;
    
    rib << "HierarchicalSubdivisionMesh \"catmull-clark\" ";
    
    rib << "[";
    std::copy(_nvertsPerFace.begin(), _nvertsPerFace.end(), std::ostream_iterator<int>(rib));
    rib << "] ";

    rib << "[";
    std::copy(_faceverts.begin(), _faceverts.end(), std::ostream_iterator<int>(rib));
    rib << "] ";
    
    std::stringstream names, nargs, intargs, floatargs, strargs;
    for (int i=0; i<(int)_tags.size();) {
        tag * t = _tags[i];
        
        names << t->name;

        nargs << t->intargs.size() << " " << t->floatargs.size() << " " << t->stringargs.size();
        
        std::copy(t->intargs.begin(), t->intargs.end(), std::ostream_iterator<int>(intargs));

        std::copy(t->floatargs.begin(), t->floatargs.end(), std::ostream_iterator<float>(floatargs));

        std::copy(t->stringargs.begin(), t->stringargs.end(), std::ostream_iterator<std::string>(strargs));
        
        if (++i<(int)_tags.size()) {
            names << " ";
            nargs << " ";
            intargs << " ";
            floatargs << " ";
            strargs << " ";
        }
    }
    
    rib << "["<<names<<"] " << "["<<nargs<<"] " << "["<<intargs<<"] " << "["<<floatargs<<"] " << "["<<strargs<<"] ";

    rib << "\"P\" [";
    std::copy(_verts.begin(), _verts.end(), std::ostream_iterator<float>(rib));
    rib << "] ";
    
    return rib.str();
}


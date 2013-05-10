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


//------------------------------------------------------------------------------
template <class T> std::string   
hbrToObj( OpenSubdiv::HbrMesh<T> * mesh ) {

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
        
        OpenSubdiv::HbrFace<T> * f = mesh->GetFace(i);
        
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
createVertices( OpenSubdivShape const * sh, OpenSubdiv::HbrMesh<T> * mesh, std::vector<float> * verts ) {

    T v;
    for(int i=0;i<sh->getNverts(); i++ ) {
        v.SetPosition( sh->verts[i*3], sh->verts[i*3+1], sh->verts[i*3+2] );
        mesh->NewVertex( i, v );
    }
}

//------------------------------------------------------------------------------
template <class T> void
createVertices( OpenSubdivShape const * sh, OpenSubdiv::HbrMesh<T> * mesh, std::vector<float> & verts ) {

    T v;
    for(int i=0;i<sh->getNverts(); i++ )
        mesh->NewVertex( i, v );
}

//------------------------------------------------------------------------------
template <class T> void
copyVertexPositions( OpenSubdivShape const * sh, OpenSubdiv::HbrMesh<T> * mesh, std::vector<float> & verts ) {

    int nverts = mesh->GetNumVertices();
    
    verts.resize( nverts * 3 );
    
    std::copy(sh->verts.begin(), sh->verts.end(), verts.begin());
    
    // Sometimes Hbr dupes some vertices during Mesh::Finish()
    if (nverts > sh->getNverts()) {
    
        for (int i=sh->getNverts(); i<nverts; ++i) {
        
            OpenSubdiv::HbrVertex<T> * v = mesh->GetVertex(i);
            
            OpenSubdiv::HbrFace<T> * f = v->GetIncidentEdge()->GetFace();
            
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
            
            verts[i*3+0] = sh->verts[shvert*3+0];
            verts[i*3+1] = sh->verts[shvert*3+1];
            verts[i*3+2] = sh->verts[shvert*3+2];
        }
    }
}

//------------------------------------------------------------------------------
template <class T> void
createTopology( OpenSubdivShape const * sh, OpenSubdiv::HbrMesh<T> * mesh, Scheme scheme) {
}

//------------------------------------------------------------------------------
template <class T> OpenSubdiv::HbrMesh<T> *
simpleHbr(char const * shapestr, Scheme scheme, std::vector<float> * verts=0) {

    OpenSubdivShape * sh = OpenSubdivShape::parseShape( shapestr );

    OpenSubdiv::HbrMesh<T> * mesh = createMesh<T>(scheme);

    createVertices<T>(sh, mesh, verts);

    createTopology<T>(sh, mesh, scheme);

    if(verts)
        copyVertexPositions<T>(sh,mesh,*verts);

    delete sh;

    return mesh;
}

//------------------------------------------------------------------------------
template <class T> OpenSubdiv::HbrMesh<T> *
simpleHbr(char const * shapestr, Scheme scheme, std::vector<float> & verts) {

    OpenSubdivShape * sh = OpenSubdivShape::parseShape( shapestr );

    OpenSubdiv::HbrMesh<T> * mesh = createMesh<T>(scheme);

    createVertices<T>(sh, mesh, verts);

    createTopology<T>(sh, mesh, scheme);

    copyVertexPositions<T>(sh,mesh,verts);

    delete sh;

    return mesh;
}



//------------------------------------------------------------------------------
// The simplest constructor, only point positions and polygonal
// mesh topology
OpenSubdivShape(const std::vector<float>  &verts,
                const std::vector<int>    &nvertsPerFace,
                const std::vector<int>    &faceverts):
    _verts[verts],
    _nvertsPerFace[nvertsPerFace],
    _faceverts[faceverts]
{


}

OpenSubdivShape(const float* points,
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

    int badIndices=0;
    for (int i=0; i<facevertsLen) {
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
    for (int i=0; i<(int)tags.size(); ++i)
        delete tags[i];

    if (_hbrMesh)
        delete _hbrMesh;
}


OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>*
OpenSubdivShape::GetHbrMesh()
{

    if (_hbrMesh) {
        return _hbrMesh;
    }

    _hbrMesh = CreateMesh();
    
    const int * fv=&(_faceverts[0]);
    for(int f=0, ptxidx=0;f< GetNfaces(); f++ ) {

        int nv = _nvertsPerFace[f];

        if ((scheme==kLoop) and (nv!=3)) {
            printf("Trying to create a Loop surbd with non-triangle face\n");
            exit(1);
        }

        for(int j=0;j<nv;j++) {
            OpenSubdiv::HbrVertex<T> * origin      = _hbrMesh->GetVertex( fv[j] );
            OpenSubdiv::HbrVertex<T> * destination = _hbrMesh->GetVertex( fv[ (j+1)%nv] );
            OpenSubdiv::HbrHalfedge<T> * opposite  = destination->GetEdge(origin);

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

        OpenSubdiv::HbrFace<T> * face = _hbrMesh->NewFace(nv, (int *)fv, 0);

        face->SetPtexIndex(ptxidx);

        if ( (scheme==kCatmark or scheme==kBilinear) and nv != 4 )
            ptxidx+=nv;
        else
            ptxidx++;

        fv+=nv;
    }

    _hbrMesh->SetInterpolateBoundaryMethod( OpenSubdiv::HbrMesh<T>::k_InterpolateBoundaryEdgeOnly );

    applyTags( _hbrMesh, sh );
    
    _hbrMesh->Finish();
    
    return _hbrMesh;
}


//------------------------------------------------------------------------------
OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>*
OpenSubdivShape::CreateMesh(OpenSubdiv::Scheme scheme)
{       
    OpenSubdiv::HbrMesh<T> * mesh = 0;

    static OpenSubdiv::HbrBilinearSubdivision<T> _bilinear;
    static OpenSubdiv::HbrLoopSubdivision<T>     _loop;
    static OpenSubdiv::HbrCatmarkSubdivision<T>  _catmark;

    switch (scheme) {
    case kBilinear : mesh = new OpenSubdiv::HbrMesh<T>( &_bilinear ); break;
    case kLoop     : mesh = new OpenSubdiv::HbrMesh<T>( &_loop     ); break;
    case kCatmark  : mesh = new OpenSubdiv::HbrMesh<T>( &_catmark  ); break;
    }

    return mesh;
}



//------------------------------------------------------------------------------
std::string OpenSubdivShape::genRIB() const {
    std::stringstream rib;
    
    rib << "HierarchicalSubdivisionMesh \"catmull-clark\" ";
    
    rib << "[";
    std::copy(nvertsPerFace.begin(), nvertsPerFace.end(), std::ostream_iterator<int>(rib));
    rib << "] ";

    rib << "[";
    std::copy(faceverts.begin(), faceverts.end(), std::ostream_iterator<int>(rib));
    rib << "] ";
    
    std::stringstream names, nargs, intargs, floatargs, strargs;
    for (int i=0; i<(int)tags.size();) {
        tag * t = tags[i];
        
        names << t->name;

        nargs << t->intargs.size() << " " << t->floatargs.size() << " " << t->stringargs.size();
        
        std::copy(t->intargs.begin(), t->intargs.end(), std::ostream_iterator<int>(intargs));

        std::copy(t->floatargs.begin(), t->floatargs.end(), std::ostream_iterator<float>(floatargs));

        std::copy(t->stringargs.begin(), t->stringargs.end(), std::ostream_iterator<std::string>(strargs));
        
        if (++i<(int)tags.size()) {
            names << " ";
            nargs << " ";
            intargs << " ";
            floatargs << " ";
            strargs << " ";
        }
    }
    
    rib << "["<<names<<"] " << "["<<nargs<<"] " << "["<<intargs<<"] " << "["<<floatargs<<"] " << "["<<strargs<<"] ";

    rib << "\"P\" [";
    std::copy(verts.begin(), verts.end(), std::ostream_iterator<float>(rib));
    rib << "] ";
    
    return rib.str();
}


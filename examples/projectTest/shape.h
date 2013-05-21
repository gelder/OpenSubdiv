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

#ifndef OSD_SHAPE_H
#define OSD_SHAPE_H


#include <osd/vertex.h>
#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>
#include <osd/cpuVertexBuffer.h>

#include <far/meshFactory.h>
#include <hbr/mesh.h>

#include <stdio.h>
#include <string.h>

#include <list>
#include <map>
#include <string>
#include <sstream>
#include <vector>



//------------------------------------------------------------------------------

class OpenSubdivShape {

    enum Scheme {
        kBilinear,
        kCatmark,
        kLoop
    };

  public:
    struct tag {

        static tag * parseTag( char const * stream );
        
        std::string genTag() const;

        std::string              name;
        std::vector<int>         intargs;
        std::vector<float>       floatargs;
        std::vector<std::string> stringargs;
    };


    OpenSubdivShape(
        std::string name,
        int numVertices,
        int maxLevels,
        std::vector<int> indices,
        std::vector<int> nverts,
        std::vector<std::string> vvNames,
        std::vector<std::string> fvNames,
        const std::vector<float>& fvData,
        tag& tagData,
        Scheme scheme);


    // Constructor using raw types.
    //
    OpenSubdivShape(const float* points,
                    int pointsLen /*= # of 3-float points*/, 
                    const int *nvertsPerFace, int numFaces,
                    const int *faceverts, int facevertsLen);

    std::string genRIB() const;

    ~OpenSubdivShape();

    int GetNverts() const { return (int)_verts.size()/3; }

    int GetNfaces() const { return (int)_nvertsPerFace.size(); }

    const std::vector<float>  &GetVerts() { return _verts;}
    const std::vector<int>    &GetNVertsPerFace() { return _nvertsPerFace;}
    const std::vector<int>    &GetFaceVerts() { return _faceverts;}
    const std::vector<tag *>  &GetTags() { return _tags;}
    
    Scheme                    &GetScheme() { return _scheme;}

    void SetCoarsePositions(const std::vector<float>& coords);
    void SetVVData(const std::vector<float>& data);
    void GetFVData(int subdivisionLevel,
                   const std::vector<std::string>& names,
                   std::vector<float>* fvdata);
    
    
  private:


    void _ProcessTagsAndFinishMesh(
        OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex> *mesh,
        std::vector<std::string> &tags,
        std::vector<int> &numArgs,
        std::vector<int> &intArgs,
        std::vector<float> &floatArgs,
        std::vector<std::string> &stringArgs);


    // Positions of points in the coarse mesh, 3 floats per point
    std::vector<float>  _verts;

    // A vector with one entry for every face (length = number of faces) in
    // the coarse mesh. Each entry is the number of verts in that face
    std::vector<int>    _nvertsPerFace;

    // Index for each point in each face into vertex buffers.  The
    // length is the sum of all the entries in nvertsPerFace
    std::vector<int>    _faceverts;

    // Annotation on the coarse mesh for things like creases. 
    std::vector<tag *>  _tags;
    Scheme              _scheme;


    // Names of vertex varying (vv) attributes
    std::vector<std::string> _vvNames;

    // Names of face varying (fv) attributes
    std::vector<std::string> _fvNames;

    // Size and indexing information for face varying data
    std::vector<int> _fvarwidths;
    std::vector<int> _fvarindices;
    std::map<std::string, int> _fvaroffsets;
    
    // The next block of member variables are the OpenSubdiv meshes,
    // buffers, and computational objects used to generate the refined
    // subdivision surface
    //

    // The lowest level mesh, it definites the polygonal topology
    // and is used for refinement
    OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>*  _hbrMesh;

    // A mesh of patches (adaptive), or quads (uniform) generated
    // by performing feature adaptive or uniform subdivision on the hbrMesh
    OpenSubdiv::FarMesh<OpenSubdiv::OsdVertex>*  _farMesh;

    // An object used to compute values based on the farMesh
    OpenSubdiv::OsdCpuComputeContext*            _computeContext;

    // Two buffers used to store computed values
    OpenSubdiv::OsdCpuVertexBuffer*              _vertexBuffer;
    OpenSubdiv::OsdCpuVertexBuffer*              _vvBuffer;
};




#endif /* OSD_SHAPE_H */

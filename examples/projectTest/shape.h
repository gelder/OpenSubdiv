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


    OpenSubdivShape();


    // returns false on error.  If errorMessage is non-NULL it'll
    // be populated upon error.
    bool Initialize(
        const std::string &name,
        int numVertices,
        int maxLevels,
        const std::vector<int> &nverts,        
        const std::vector<int> &indices,
        const std::vector<std::string> &vvNames,
        const std::vector<std::string> &fvNames,
        const std::vector<float> &fvData,
        const tag& tagData,
        Scheme scheme,
        std::string *errorMessage = NULL);



    // Initialize using raw types. 
    //
    // This is useful for automated tests initializing with data like:
    // int nverts[] = { 4, 4, 4, 4, 4, 4};
    //
    bool Initialize(
        int numVertices,
        const int *nverts, int numFaces,
        const int *indices, int indicesLen,
        int maxLevels,
        std::string *errorMessage = NULL);
    
    std::string genRIB() const;

    ~OpenSubdivShape();

    int GetNverts() const { return (int)_verts.size()/3; }

    int GetNfaces() const { return (int)_nverts.size(); }

    const std::vector<float>  &GetVerts() { return _verts;}
    const std::vector<int>    &GetNVerts() { return _nverts;}
    const std::vector<int>    &GetIndices() { return _indices;}
    const std::vector<tag *>  &GetTags() { return _tags;}
    
    Scheme                    &GetScheme() { return _scheme;}

    void SetCoarsePositions(const std::vector<float>& coords);
    void SetVVData(const std::vector<float>& data);
    void GetFVData(int subdivisionLevel,
                   const std::vector<std::string>& names,
                   std::vector<float>* fvdata);

    
    // The Refine method uses OpenSubdiv to refine all subdiv nodes in the held
    //  tidscene and cache the results.
    //
    //  If an upperLimit is set, this halts subdivision of each gprim so that
    //  the number of vertices never exceeds it.  This feature is useful
    //  for renderers like OpenGL ES that only support 16-bit indices.
    //
    bool Refine(int levels, int upperLimit = (int)0x7fffffff,
                std::string *errorMessage = NULL);

    // Fetch the XYZ coordinates of the post-refined vertices. 
    bool GetPositions(std::vector<float>* positions,
                      std::string *errorMessage = NULL);

    // Fetch the topology of the post-refined mesh. 
    bool GetQuads(std::vector<int>* quads,
                  std::string *errorMessage = NULL);

    // Write the refined quad mesh to given filename, return false on error
    bool WriteRefinedObj( const std::string &filename,
                          std::string *errorMessage = NULL);


private:
    
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
    std::vector<int>    _nverts;

    // Index for each point in each face into vertex buffers.  The
    // length is the sum of all the entries in nverts
    std::vector<int>    _indices;

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

    // Objects used to compute values based on the farMesh
    OpenSubdiv::OsdCpuComputeContext*            _computeContext;
    OpenSubdiv::OsdCpuComputeController*         _computeController;

    // Two buffers used to store computed values
    OpenSubdiv::OsdCpuVertexBuffer*              _vertexBuffer;
    OpenSubdiv::OsdCpuVertexBuffer*              _vvBuffer;



    int _firstVertexOffset;
    int _numRefinedVerts;
    int _numRefinedQuads;
    int _refinedLevel;
    int _maxLevels;

};




#endif /* OSD_SHAPE_H */

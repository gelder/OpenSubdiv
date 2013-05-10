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

#include <hbr/mesh.h>

#include <stdio.h>
#include <string.h>

#include <list>
#include <string>
#include <sstream>
#include <vector>

//------------------------------------------------------------------------------
static char const * sgets( char * s, int size, char ** stream ) {
    for (int i=0; i<size; ++i) {
        if ( (*stream)[i]=='\n' or (*stream)[i]=='\0') {

            memcpy(s, *stream, i);
            s[i]='\0';

            if ((*stream)[i]=='\0')
                return 0;
            else {
                (*stream) += i+1;
                return s;
            }
        }
    }
    return 0;
}


//------------------------------------------------------------------------------
enum Scheme {
  kBilinear,
  kCatmark,
  kLoop
};

//------------------------------------------------------------------------------

class OpenSubdivShape {

  public:
    struct tag {

        static tag * parseTag( char const * stream );
        
        std::string genTag() const;

        std::string              name;
        std::vector<int>         intargs;
        std::vector<float>       floatargs;
        std::vector<std::string> stringargs;
    };

    // The simplest constructor, only point positions and polygonal
    // mesh topology
    OpenSubdivShape(const std::vector<float>  &points,
                    const std::vector<int>    &nvertsPerFace,
                    const std::vector<int>    &faceverts);

    // Constructor using raw types.
    //
    OpenSubdivShape(const float* points,
                    int pointsLen /*= # of 3-float points*/, 
                    const int *nvertsPerFace, int numFaces,
                    const int *faceverts, int facevertsLen);

    std::string genRIB() const;

    ~OpenSubdivShape();

    // Return the hbr mesh, allocating and filling it if needed
    OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex> *GetHbrMesh();
    
    int GetNverts() const { return (int)_verts.size()/3; }

    int GetNfaces() const { return (int)_nvertsPerFace.size(); }

    
    
  private:

    OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>  *CreateMesh(
        OpenSubdiv::Scheme scheme = kCatmark);
    
    std::vector<float>  _verts;
    std::vector<float>  _uvs;
    std::vector<float>  _normals;
    std::vector<int>    _nvertsPerFace;
    std::vector<int>    _faceverts;
    std::vector<int>    _faceuvs;
    std::vector<int>    _facenormals;
    std::vector<tag *>  _tags;
    Scheme              _scheme;


    // Initialized to NULL, only allocated if feature adaptive refinement
    // is required, not if subdivision/patch tables are loaded.
    //
    OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>  *_hbrMesh;
    
};




#endif /* OSD_SHAPE_H */

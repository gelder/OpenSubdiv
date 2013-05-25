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


#include <far/meshFactory.h>
#include <osd/vertex.h>

#include <fstream>

using namespace OpenSubdiv;
using namespace std;


//------------------------------------------------------------------------------
// The simplest constructor, only point positions and polygonal
// mesh topology
/*
OpenSubdivShape::OpenSubdivShape(const std::vector<float>  &verts,
                                 const std::vector<int>    &nverts,
                                 const std::vector<int>    &indices) :
    _verts(verts),
    _nverts(nverts),
    _indices(indices),
    _hbrMesh(NULL)
{

}
*/

OpenSubdivShape::OpenSubdivShape():
    _scheme(kCatmark),
    _hbrMesh(NULL),
    _farMesh(NULL),
    _computeContext(NULL),
    _computeController(NULL),
    _vertexBuffer(NULL),
    _vvBuffer(NULL),
    _firstVertexOffset(0),
    _numRefinedVerts(0),
    _numRefinedQuads(0),    
    _refinedLevel(0),
    _maxLevels(0)
{
}


bool
OpenSubdivShape::Initialize(
    int numVertices,
    const int *nverts, int numFaces,
    const int *indices, int indicesLen,
    int maxLevels,
    string *errorMessage)
{

    std::vector<int> nvertsVec, indicesVec;
    nvertsVec.resize(numFaces);
    for (int i=0; i<(int)nvertsVec.size(); ++i) {
        nvertsVec[i] = nverts[i];
    }
    indicesVec.resize(indicesLen);    
    for (int i=0; i<(int)indicesVec.size(); ++i) {
        indicesVec[i] = indices[i];
    }    

    return Initialize("test", numVertices, maxLevels, nvertsVec, indicesVec,
                      vector<string>(), vector<string>(), vector<float>(), tag(),
                      kCatmark, errorMessage);
}
    


bool
OpenSubdivShape::Initialize(
    const string &name,
    int numVertices,
    int maxLevels,
    const vector<int> &nverts,
    const vector<int> &indices,
    const vector<string> &vvNames,
    const vector<string> &fvNames,
    const vector<float>& fvData,
    const tag& tagData,
    Scheme scheme,
    std::string *errorMessage)
{

    // Check the polygonal topology for errors before copying 
    int totalNumIndices=0;
    for (int i=0; i<(int)nverts.size(); ++i) {
        totalNumIndices += nverts[i];
    }

    if (totalNumIndices != (int)indices.size()) {
        if (errorMessage)
            *errorMessage = "Bad indexing for face topology";
        return false;
    }

    // nverts is OK
    _nverts = nverts;

    for (int i=0; i<(int)indices.size(); ++i) {
        if ((indices[i] < 0) or (indices[i] >= numVertices)) {
            if (errorMessage)            
                *errorMessage = "Bad index in face topology";
            _nverts.clear();
            return false;
        } 
    }

    // indices is OK, copy the rest of the data and begin
    // initalization
    _indices = indices;
    _scheme = scheme;
    _vvNames = vvNames;
    _fvNames = fvNames;
    _maxLevels = maxLevels;

    // Decide on a subdivision scheme
    HbrSubdivision<OsdVertex> *schemeToUse = NULL;
    static HbrBilinearSubdivision<OsdVertex> _bilinear;
    static HbrLoopSubdivision<OsdVertex>     _loop;
    static HbrCatmarkSubdivision<OsdVertex>  _catmark;
    
    switch (scheme) {
    case kBilinear : schemeToUse = &_bilinear; break;
    case kLoop     : schemeToUse = &_loop; break;
    case kCatmark  : schemeToUse = &_catmark; break;
    }

    
    //_name = name;

    // Construct the Hbr mesh representation that will hold the
    // coarse mesh and be refined by Far 
    if (fvNames.empty()) {
        _hbrMesh = new HbrMesh<OsdVertex>(schemeToUse);
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

        _hbrMesh = new HbrMesh<OsdVertex>(schemeToUse, fvarcount,
                                          &_fvarindices[0],
                                          &_fvarwidths[0], fvarcount);
    }
    

    // Create vertices in the hbr mesh
    OsdVertex v;
    for (int i = 0; i < numVertices; ++i) {
        HbrVertex<OsdVertex>* hvert = _hbrMesh->NewVertex(i, v);
        if (!hvert) {
        if (errorMessage)            
            printf("Unable to create OSD vertex from tidscene");
        }
    }

    // Sanity check on face varying data
    int fvarWidth = _hbrMesh->GetTotalFVarWidth();
    if (fvData.size() < nverts.size() * fvarWidth ||
        fvarWidth != (int)fvNames.size()) {
        if (errorMessage)        
        printf("Incorrectly sized face data: name count = %d, "
                       "data width = %d, face count = %d, total data size = %d.",
                       (int) fvNames.size(),
                       fvarWidth,
                       (int) nverts.size(),
                       (int) fvData.size());
    }

    // Create faces
    const int * fv=&(_indices[0]);
    // face-vertex count offset
    int fvcOffset = 0;    
    for(int f=0, ptxidx=0;f< GetNfaces(); f++ ) {

        int nv = _nverts[f];

        if ((_scheme==kLoop) and (nv!=3)) {
        if (errorMessage)            
            printf("Trying to create a Loop surbd with non-triangle face\n");
            exit(1);
        }

        for(int j=0;j<nv;j++) {
            HbrVertex<OsdVertex> * origin      = _hbrMesh->GetVertex( fv[j] );
            HbrVertex<OsdVertex> * destination = _hbrMesh->GetVertex( fv[ (j+1)%nv] );
            HbrHalfedge<OsdVertex> * opposite  = destination->GetEdge(origin);

            if(origin==NULL || destination==NULL) {
        if (errorMessage)                
                printf(" An edge was specified that connected a nonexistent vertex\n");
                exit(1);
            }

            if(origin == destination) {
        if (errorMessage)                
                printf(" An edge was specified that connected a vertex to itself\n");
                exit(1);
            }

            if(opposite && opposite->GetOpposite() ) {
        if (errorMessage)                
                printf(" A non-manifold edge incident to more than 2 faces was found\n");
                exit(1);
            }

            if(origin->GetEdge(destination)) {
        if (errorMessage)                
                printf(" An edge connecting two vertices was specified more than once."
                       " It's likely that an incident face was flipped\n");
                exit(1);
            }
        }

        HbrFace<OsdVertex> * face = _hbrMesh->NewFace(nv, (int *)fv, 0);

        face->SetPtexIndex(ptxidx);

        if (!face) {
        if (errorMessage)            
            printf("Unable to create OSD face from tidscene");
            exit(1);
        }

        // prideout: 3/21/2013 - Inspired by "GetFVarData" in examples/mayaViewer/hbrUtil.cpp
        if (!fvNames.empty()) {
            const float* faceData = &(fvData[fvcOffset*fvarWidth]);
            for (int fvi = 0; fvi < nv; ++fvi) {
                int vindex = indices[fvi + fvcOffset];
                HbrVertex<OsdVertex>* v = _hbrMesh->GetVertex(vindex);
                HbrFVarData<OsdVertex>& fvarData = v->GetFVarData(face);
                if (!fvarData.IsInitialized()) {
                    fvarData.SetAllData(fvarWidth, faceData); 
                } else if (!fvarData.CompareAll(fvarWidth, faceData)) {

                    // If data exists for this face vertex, but is different
                    // (e.g. we're on a UV seam) create another fvar datum
                    HbrFVarData<OsdVertex>& fvarData = v->NewFVarData(face);
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


//    _ProcessTagsAndFinishMesh(
//        _hbrMesh, tagData.tags, tagData.numArgs, tagData.intargs,
//        tagData.floatargs, tagData.stringargs);
    _hbrMesh->Finish();
    
    // create the quad tables to include all levels by specifying firstLevel as 1
    FarMeshFactory<OsdVertex> meshFactory(_hbrMesh, maxLevels, false, /*firstLevel=*/1);

    _farMesh = meshFactory.Create();
    _computeContext = OsdCpuComputeContext::Create(_farMesh);
    _vertexBuffer = OsdCpuVertexBuffer::Create(
        3, _farMesh->GetNumVertices());

    if (vvNames.empty()) {
        _vvBuffer = NULL;
    } else {
        _vvBuffer = OsdCpuVertexBuffer::Create(
            vvNames.size(), _farMesh->GetNumVertices());
    }

    _computeController = new OpenSubdiv::OsdCpuComputeController();

    return true;
}

void
OpenSubdivShape::SetCoarsePositions(const vector<float>& coords)
{
    //XXX: add sanity check on size
    
    const float* pFloats = &coords.front();
    int numFloats = (int) coords.size();
    _vertexBuffer->UpdateData(pFloats, 0, numFloats / 3);
}

void
OpenSubdivShape::SetVVData(const vector<float>& data)
{
    if (!_vvBuffer) {
        if (!data.empty()) {
            printf("Mesh was not constructed with VV variables.");
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
    for (int i=0; i<(int)names.size(); ++i) {
        if (_fvaroffsets.find(names[i]) == _fvaroffsets.end()) {
            printf("Can't find facevarying variable %s\n", names[i].c_str());
            return;
        }
    }

    // Fetch *all* faces; this includes all subdivision levels.
    vector<HbrFace<OsdVertex> *> faces;
    _hbrMesh->GetFaces(std::back_inserter(faces));

    // Iterate through all faces, filtering on the requested subdivision level.
    for (int faceIndex = 0; faceIndex < (int)faces.size(); ++faceIndex) {
        if (faces[faceIndex]->GetDepth() != level) {
            continue;
        }
        int ncorners = faces[faceIndex]->GetNumVertices();
        for (int corner = 0; corner < ncorners; ++corner) {
            HbrFVarData<OsdVertex>& fvariable = faces[faceIndex]->GetFVarData(corner);
            for (int i=0; i<(int)names.size(); ++i) {            
                int offset = _fvaroffsets[names[i]];
                const float* data = fvariable.GetData(offset);
                outdata->push_back(*data);
            }
        }
    }
}


bool
OpenSubdivShape::Refine(int level, int upperLimit, string *errorMessage)
{
  if (level > _maxLevels) {
      printf("Refinement level exceeds the maximum.");
      return false;
  }   

    _computeController->Refine(_computeContext,
                               _farMesh->GetKernelBatches(),
                               _vertexBuffer,        
                               _vvBuffer);   
    
    // Find a subdivision level that is less dense than the upper limit:
    const FarSubdivisionTables<OsdVertex>* ftable = _farMesh->GetSubdivisionTables();
    _refinedLevel = level;
    while (_refinedLevel > 0) {
        // GetNumVertices() returns number of vertices at a given level.
        // GetNumVerticesTotal() returns total number of vertices at
        //   a given level including lower levels.
        int levelSize = ftable->GetNumVertices(_refinedLevel);
        if (levelSize < upperLimit) {
            break;
        }
        _refinedLevel--;
    }

    std::cout << "_refinedLevel = " << _refinedLevel << endl;
    
    // Check if we can't can squeeze it into the upper bound:
    if (_refinedLevel == 0) {
        printf("Upper limit on vertex density is too small\n");
        return false;
    }       
    
    // Find quads array at _refinedLevel
    const FarPatchTables * ptables = _farMesh->GetPatchTables();
    const FarPatchTables::PatchArrayVector & parrays = ptables->GetPatchArrayVector();
    if (_refinedLevel > (int)parrays.size()) {
        printf("Invalid size of patch array %d %d\n", _refinedLevel, (int)parrays.size());;
        return false;
    }

    // parrays doesn't contain base mesh, so it starts with level==1
    const FarPatchTables::PatchArray & parray = parrays[_refinedLevel-1];

    // Populate the metadata for the aggregate API:
    _firstVertexOffset = ftable->GetFirstVertexOffset(_refinedLevel);
    std::cout << "_firstVertexOffset = " << _firstVertexOffset << endl;
    
    _numRefinedQuads = (int) parray.GetNumPatches();
    _numRefinedVerts = (int) ftable->GetNumVertices(_refinedLevel);

    // Keep track of the float offsets for this mesh, for all UV data.
    /*
    int numVVFloats = 0;
    TF_FOR_ALL(it, _vvNames) {
        string name = *it;
        meta["vv_offset_" + name] = numVVFloats;
        numVVFloats++;
    }
    TF_FOR_ALL(it, _vvNames) {
        string name = *it;
        meta["vv_stride_" + name] = numVVFloats;
    }
    int numFVFloats = 0;
    TF_FOR_ALL(it, _fvNames) {
        string name = *it;
        meta["fv_offset_" + name] = numVVFloats;
        numFVFloats++;
    }
    TF_FOR_ALL(it, _fvNames) {
        string name = *it;
        meta["fv_stride_" + name] = numVVFloats;
    }
    
    _numRefinedQuads += meta["quadsCount"].Get<int>();
    _numRefinedVerts += meta["vertsCount"].Get<int>();
    */
    return true;
}




bool
OpenSubdivShape::GetPositions(vector<float>* coords,
                              string *errorMessage)
{
    if (!coords || (_numRefinedVerts == 0)) {
        return false;
    }
    
    coords->resize(_numRefinedVerts * 3); // 3 floats/point
    const float* cBegin = _vertexBuffer->BindCpuBuffer();
    for (int i=0; i<_numRefinedVerts*3; ++i) {
        (*coords)[i] = cBegin[i + (_firstVertexOffset*3)];
    }

    return true;
}


bool
OpenSubdivShape::GetQuads(vector<int>* quads,
                          string *errorMessage)
{
    if (!quads || (_numRefinedQuads == 0)) {
        return false;
    }
    
    quads->resize(_numRefinedQuads * 4);

    const FarPatchTables * ptables = _farMesh->GetPatchTables();
    const FarPatchTables::PatchArrayVector & parrays = ptables->GetPatchArrayVector();
    if (_refinedLevel > (int)parrays.size()) {
        printf("Invalid size of patch array\n");
        return false;
        
    }
    const FarPatchTables::PatchArray & parray = parrays[_refinedLevel-1];

    const unsigned int *pQuads = &ptables->GetPatchTable()[parray.GetVertIndex()];
    size_t numInts = parray.GetNumPatches() * 4;
    for (int i=0; i<(int)numInts; ++i) {
        (*quads)[i] = pQuads[i] - _firstVertexOffset;
    }

    return true;
}


bool
OpenSubdivShape::WriteRefinedObj( const std::string &filename,
                                std::string *errorMessage)
{

    // open a file stream
    std::ofstream fout(filename.c_str());
    if (!fout.is_open()) {
        if (errorMessage)
            printf("Could not open file: (%s)\n", filename.c_str());
        return false;
    }


    std::vector<float> refinedPositions;
    std::vector<int> refinedQuads;
    if (not (GetPositions(&refinedPositions, errorMessage) and
             GetQuads(&refinedQuads, errorMessage))) {
        return false;
    }
    
    fout << "# Exported by OpenSubdivShape\n";

    fout << "#Positions = " << refinedPositions.size()/3 << std::endl;        
    for (int i=0; i<(int)refinedPositions.size(); i+=3)  {
        fout << "v " << refinedPositions[i] <<
            " " << refinedPositions[i+1] <<
            " " << refinedPositions[i+2] << "\n";               
    }

    fout << "# Quads = " << refinedQuads.size()/4 << std::endl;        
    for (int i=0; i<(int)refinedQuads.size(); i+=4)  {
        // OBJ face indices are 1 based
        fout << "f " << refinedQuads[i]+1 <<
            " " << refinedQuads[i+1]+1 <<
            " " << refinedQuads[i+2]+1 <<
            " " << refinedQuads[i+3]+1 <<
            "\n";
    }


    fout.close();
    
    return true;
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


// ProcessTagsAndFinishMesh(...)
// This translates prman-style lists of tags into OSD method calls.
//
// prideout: 3/19/2013 - since tidSceneRenderer has a similar
//           function, we should factor this into an amber utility, or
//           into osd itself.  I'd vote for the latter.  It already has 
//           a shapeUtils in its regression suite that almost fits the bill.
//
// prideout: 3/19/2013 - edits are not yet supported.
void
OpenSubdivShape::_ProcessTagsAndFinishMesh(
    HbrMesh<OsdVertex> *mesh,
    vector<string> &tags,
    vector<int> &numArgs,
    vector<int> &intArgs,
    vector<float> &floatArgs,
    vector<string> &stringArgs)
{
    mesh->SetInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryEdgeOnly);

    int* currentInt = &intArgs[0];
    float* currentFloat = &floatArgs[0];
    string* currentString = &stringArgs[0];

    // TAGS (crease, corner, hole, smooth triangles, edits(vertex,
    // edge, face), creasemethod, facevaryingpropagatecorners, interpolateboundary
    for(int i = 0; i < (int)tags.size(); ++i){
	const char * tag = tags[i].c_str();
	int nint = numArgs[3*i];
	int nfloat = numArgs[3*i+1];
	int nstring = numArgs[3*i+2];

	// XXX could use tokens here to reduce string matching overhead
	if(strcmp(tag, "interpolateboundary") == 0) {
	    // Interp boundaries
	    assert(nint == 1);
	    switch(currentInt[0]) {
            case 0:
                mesh->SetInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryNone);
                break;
            case 1:
                mesh->SetInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryEdgeAndCorner);
                break;
            case 2:
                mesh->SetInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryEdgeOnly);
                break;
            default:
                printf("Subdivmesh contains unknown interpolate boundary method: %d\n",
                       currentInt[0]);
		break;
	    }
	    // Processing of this tag is done in mesh->Finish()
	} else if(strcmp(tag, "crease") == 0) {
	    for(int j = 0; j < nint-1; ++j) {
		// Find the appropriate edge
                HbrVertex<OsdVertex>* v = mesh->GetVertex(currentInt[j]);
                HbrVertex<OsdVertex>* w = mesh->GetVertex(currentInt[j+1]);
                HbrHalfedge<OsdVertex>* e = NULL;
		if(v && w) {
		    e = v->GetEdge(w);
		    if(!e) {
			// The halfedge might be oriented the other way
			e = w->GetEdge(v);
		    }
		}
		if(!e) {
		    printf("Subdivmesh has non-existent sharp edge (%d,%d).\n",
                            currentInt[j], currentInt[j+1]);
		} else {
		    e->SetSharpness(std::max(0.0f, ((nfloat > 1) ? currentFloat[j] : currentFloat[0])));
		}
	    }
	} else if(strcmp(tag, "corner") == 0) {
	    for(int j = 0; j < nint; ++j) {
                HbrVertex<OsdVertex>* v = mesh->GetVertex(currentInt[j]);
		if(v) {
		    v->SetSharpness(std::max(0.0f, ((nfloat > 1) ? currentFloat[j] : currentFloat[0])));
		} else {
		    printf("Subdivmesh has non-existent sharp vertex %d.\n", currentInt[j]);
		}
	    }
	} else if(strcmp(tag, "hole") == 0) {
	    for(int j = 0; j < nint; ++j) {
                HbrFace<OsdVertex>* f = mesh->GetFace(currentInt[j]);
		if(f) {
		    f->SetHole();
		} else {
		    printf("Subdivmesh has hole at non-existent face %d.\n",
                            currentInt[j]);
		}
	    }
	} else if(strcmp(tag, "facevaryinginterpolateboundary") == 0) {
	    switch(currentInt[0]) {
            case 0:
                mesh->SetFVarInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryNone);
                break;
            case 1:
		mesh->SetFVarInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryEdgeAndCorner);
		break;
            case 2:
		mesh->SetFVarInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryEdgeOnly);
		break;
            case 3:
		mesh->SetFVarInterpolateBoundaryMethod(HbrMesh<OsdVertex>::k_InterpolateBoundaryAlwaysSharp);
		break;
            default:
		printf("Subdivmesh contains unknown facevarying interpolate "
                        "boundary method: %d.\n", currentInt[0]);
		break;
	    }
	} else if(strcmp(tag, "smoothtriangles") == 0) {
	    // Do nothing - CatmarkMesh should handle it
	} else if(strcmp(tag, "creasemethod") == 0) {
	    if(nstring < 1) {
		printf("Creasemethod tag missing string argument on SubdivisionMesh.\n");
	    } else {
                OpenSubdiv::HbrSubdivision<OpenSubdiv::OsdVertex>* subdivisionMethod = mesh->GetSubdivision();
		if(strcmp(currentString->c_str(), "normal") == 0) {
		    subdivisionMethod->SetCreaseSubdivisionMethod(
                        OpenSubdiv::HbrSubdivision<OpenSubdiv::OsdVertex>::k_CreaseNormal);
		} else if(strcmp(currentString->c_str(), "chaikin") == 0) {
		    subdivisionMethod->SetCreaseSubdivisionMethod(
                        OpenSubdiv::HbrSubdivision<OpenSubdiv::OsdVertex>::k_CreaseChaikin);		    
		} else {
		    printf("Creasemethod tag specifies unknown crease "
                            "subdivision method '%s' on SubdivisionMesh.\n",
                            currentString->c_str());
		}
	    }
	} else if(strcmp(tag, "facevaryingpropagatecorners") == 0) {
	    if(nint != 1) {
		printf("Expecting single integer argument for "
                        "\"facevaryingpropagatecorners\" on SubdivisionMesh.\n");
	    } else {
		mesh->SetFVarPropagateCorners(currentInt[0] != 0);
	    }
	} else if(strcmp(tag, "vertexedit") == 0
		  || strcmp(tag, "edgeedit") == 0) {
	    // XXX DO EDITS
            printf("vertexedit and edgeedit not yet supported.\n");
	} else {
	    // Complain
            printf("Unknown tag: %s.\n", tag);
	}

	// update the tag data pointers
	currentInt += nint;
	currentFloat += nfloat;
	currentString += nstring;
    }
    mesh->Finish();
}


//------------------------------------------------------------------------------
std::string OpenSubdivShape::genRIB() const {
    std::stringstream rib;
    
    rib << "HierarchicalSubdivisionMesh \"catmull-clark\" ";
    
    rib << "[";
    std::copy(_nverts.begin(), _nverts.end(), std::ostream_iterator<int>(rib));
    rib << "] ";

    rib << "[";
    std::copy(_indices.begin(), _indices.end(), std::ostream_iterator<int>(rib));
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


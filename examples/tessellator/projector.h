#ifndef WEIGHTOBJECT_PROJECTOR_H
#define WEIGHTOBJECT_PROJECTOR_H

#include "Tessellator.h"
#include <osd/evalLimitContext.h>

#include "bedrock/gf/vec3d.h"
#include <vector>
#include <boost/shared_ptr.hpp>



/// Shared pointer type for storing a OpenSubdivTessellator object.
typedef boost::shared_ptr<OpenSubdivTessellator> OpenSubdivTessellatorSharedPtr;

class WeightObject_Projector  
{
  public:

    // Simple struct to hold three indices and three floats to indicate
    // parametric position inside a triangle
    struct EvalCoords {
    public:
        EvalCoords() {}
        EvalCoords(int index0, int index1, int index2,
                   float b0, float b1, float b2) {
            indices[0] = index0; indices[1] = index1; indices[2] = index2;
            baryCentric[0] = b0; baryCentric[1] = b1; baryCentric[2] = b2; 
        }
        
        int indices[3];
        float baryCentric[3];
    };


   bool Initialize(
       const OpenSubdivTessellatorSharedPtr &tessellator,       
       int subdivLevel);


   // Project points onto the uniformly subdivided surface given by
   // the tessellator
   bool ProjectPoints(
       std::vector<WeightObject_Projector::EvalCoords> *pointsProjectedOntoSubdiv,
       const std::vector<GfVec3d> &coarsePoints,    
       const std::vector<GfVec3d> &pointsToProject);

  private:
   OpenSubdivTessellatorSharedPtr _tessellator;
   int _subdivLevel;
};


#endif // WEIGHTOBJECT_PROJECTOR_H

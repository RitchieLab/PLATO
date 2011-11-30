// MDRPDTRepository.h

#ifndef __MDRPDTREPOSITORY_H__
#define __MDRPDTREPOSITORY_H__

#include "familyrepository.h"
#include "../method_lib/DataSet.h"
// #include "fake.h"


namespace MdrPDT {
// template<class gtConv>
// class FamilyRepository

class MDRPDTRepository: public MdrPDT::FamilyRepository<MdrPDT::GenotypeConversion>{

  public:
    
    /// Load data into repository
    void Load(Methods::DataSet* ds);
    
    /// Returns the index of the locus in the repository
    inline int getNewMarkerIndex(int l){return marker_map_conversion[l];}
    
    /// Returns number of DSPs
    unsigned int PostLoad();
    
  private:

    map<int, int> marker_map_conversion;

};

}
#endif

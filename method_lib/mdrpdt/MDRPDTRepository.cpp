//MDRPDTRepository.cpp

#include "MDRPDTRepository.h"
#include "pdtmodel.h"

namespace MdrPDT{

///
/// Loads data into repository
/// @param ds DataSet to load from
///
void MDRPDTRepository::Load(Methods::DataSet* ds){
  
  Methods::Sample* set_samp;
  
  // need to store 
  vector<int>* mark_map = ds->get_marker_map();
  // map to store 
  marker_map_conversion.clear();
  // vector map included 
  vector<int> mark_included;
  
  string al1, al2;
  
  // store indexes of included markers here in marker map order
  // to match expectation of processing
  for(unsigned int markindx=0; markindx < ds->num_loci(); markindx++){
//     if(ds->get_locus(mark_map[markindx])->isEnabled()){
    if(ds->get_locus(markindx)->isEnabled()){
      mark_included.push_back(markindx);
      marker_map_conversion[(*mark_map)[markindx]] = mark_included.size()-1; // store index in map for future use
    }
  }
  
  this->snpCount = mark_included.size();
  
  Individual* newIndividual;
  char status, gender;
  int genotype;
  
  for(int indindx=0; indindx < ds->num_inds(); indindx++){
    set_samp = ds->get_sample(indindx);
    if(!(set_samp->isEnabled()))
      continue; // skip this ind
    newIndividual = GetIndividual(set_samp->getFamID().c_str(), set_samp->getInd().c_str(), true);
// cout << "fam = " << set_samp->getFamID().c_str() << " id = " << set_samp->getInd().c_str();
// cout << " mom = " << set_samp->getMomID() << " dad=" << set_samp->getDadID() << endl;
    newIndividual->SetParentIDs(set_samp->getDadID(), set_samp->getMomID());
    if(set_samp->getAffected())
      status = '2';
    else
      status = '1';
    
    newIndividual->SetStatus(status);
    
    if(set_samp->getSex())
      gender = '1';
    else
      gender = '2';
    
    newIndividual->SetGender(gender);
    assert(newIndividual);
    
    // add genotypes, using only those in marker list
    vector<int>::iterator iter;
    for(iter = mark_included.begin(); iter != mark_included.end(); ++iter){
      genotype = set_samp->get_genotype(*iter);
      
      switch(genotype){
        case 0:
          al1 = "1";
          al2 = "1";
          break;
        case 1:
          al1 = "1";
          al2 = "2";
          break;
        case 2:
          al1 = "2";
          al2 = "2";
          break;
        case 3:
          al1 = "0";
          al2 = "0";
          break;
      };  
      newIndividual->AddGenotype(convertGenotypes(al1.c_str(), al2.c_str(), 1));
    }
     
  }
  
  // need to reconcile the pedigrees and add snp info
	Iterator itr = GetIterator();
	Pedigree *ped = itr.GetNext();

	while (ped) {
		ped->Reconcile();
		ped = itr.GetNext();
	}

	char label[64];
	unsigned int maxWidth = PdtModel::MaxLabelLength;
	//Unless we have a dat file, we'll just keep a simple integer based id
	for (int i=0; i<this->snpCount; i++) {
		sprintf(label, "%d", i+1);
		snpLabels.push_back(label);
		if (strlen(label) > maxWidth)
			maxWidth = strlen(label);
	}
	PdtModel::MaxLabelLength = maxWidth;  
  
  
}

///
/// Same as base class function except no output
/// @return Returns number of DSPs
/// 
unsigned int MDRPDTRepository::PostLoad(){
	std::map<std::string, Pedigree*>::iterator itr = pedigrees.begin();
	std::map<std::string, Pedigree*>::iterator end = pedigrees.end();

	dspCount = 0;
	while (itr != end) {
		itr->second->PostLoad();
		dspCount += itr->second->GetDSPCount();
//		itr->second->PostLoad();
		itr++;
	}
// 	if (dspCount < 1) {
// 		cerr<<"There are not enough affected individuals present to perform an analysis. Please check that affected status is designated according to: A)"
// 			<<Individual::AffectedValue<<"\tU)"<<Individual::UnaffectedValue<<". See the manual for instructions on how to correctly designate affected status using the configuration file.\n";
// 		exit(1);
// 	}
	return dspCount;  
}


}

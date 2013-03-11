//FilterFactory.cpp
#include "FilterFactory.h"
#include "Filters.h"
namespace Filters{
std::map<std::string, FilterFactory::FilterType> FilterFactory::FilterMap;

///
/// Function that creates a filter based on the name
/// Used by system code in managing analysis pipeline
/// @param filterName name of filter
/// @return pointer to new Filter object
/// @throws FilterExcept if not a valid filter name
///
Filter * FilterFactory::create_filter(std::string filterName){
  if(FilterMap.empty()){
    SetFilterMap();
  }

  Filter * newFilter;
  switch(FilterMap[filterName]){
    case MissingFilterType:
      throw FilterExcept(filterName + " is not a valid filter name.");
      break;
    case ChiSquareFilterType:
      newFilter = new ChiSquareFilter(filterName);
      break;
//    case SeqRepFilterType:
//      newFilter = new SeqRepFilter(filterName);
//      break;
    case MDRFilterType:
      newFilter = new MDRFilter(filterName);
      break;
    case LogRegressFilterType:
      newFilter = new LogRegFilter(filterName);
      break;
//     case TAILSFilterType:
//       newFilter = new TAILSFilter(filterName);
//       break;
//     case ChiAlleleFilterType:
//       newFilter = new ChiAlleleFilter(filterName);
//       break;
//     case FreqCountFilterType:
//       newFilter = new FreqCountFilter(filterName);
//       break;
    case NMIFilterType:
      newFilter = new NMIFilter(filterName);
      break;
    case ContingencyFilterType:
      newFilter = new ContingencyFilter(filterName);
      break;
    case ArmitageFilterType:
      newFilter = new ArmitageFilter(filterName);
      break;
    case OddsRatioFilterType:
      newFilter = new OddsRatioFilter(filterName);
      break;
    case LikelihoodFilterType:
      newFilter = new LikelihoodRatioFilter(filterName);
      break;
    case UncertaintyCoeffFilterType:
      newFilter = new UncertaintyCoeffFilter(filterName);
      break;
    case ConditionalLRType:
      newFilter = new ConditionalLRFilter(filterName);
      break;
    case BioFilterType:
      newFilter = new BioFilter(filterName);
      break;
    case MarkerGenoType:
      newFilter = new MarkerGenoFilter(filterName);
      break;
    case AlleleFreqFilterType:
      newFilter = new AlleleFrequencyFilter(filterName);
      break;
    case MDRPDTFilterType:
      newFilter = new MDRPDTFilter(filterName);
      break;
    case MaxScoreFilterType:
      newFilter = new MaxScoreFilter(filterName);
      break;
    case PermScoreFilterType:
      newFilter = new PermScoreFilter(filterName);
      break;
#ifdef HAVE_R
    case MarsFilterType:
    	newFilter = new MarsFilter(filterName);
    	break;
#endif
    default:
      throw FilterExcept(filterName + " is not a valid filter name.");
  };
  return newFilter;
}


///
/// Establishes the map for use in creating filters
/// @return
///
void FilterFactory::SetFilterMap(){
  FilterMap["CHISQUARE"]=ChiSquareFilterType;
//  FilterMap["SERF"]=SeqRepFilterType;
  FilterMap["MDR"]=MDRFilterType;
  FilterMap["LOGISTICREGRESS"]=LogRegressFilterType;
//  FilterMap["TAILS"]=TAILSFilterType;
  FilterMap["CHIALLELE"] = ChiAlleleFilterType;
//   FilterMap["ALLELEFREQ"] = FreqCountFilterType;
  FilterMap["NMI"] = NMIFilterType;
  FilterMap["CONTINGENCY"] = ContingencyFilterType;
  FilterMap["ARMITAGE"] = ArmitageFilterType;
  FilterMap["ODDSRATIO"] = OddsRatioFilterType;
  FilterMap["LIKELIHOODRATIO"] = LikelihoodFilterType;
  FilterMap["UNCERTAINTYCOEFF"] = UncertaintyCoeffFilterType;
  FilterMap["CONDITIONALLR"] = ConditionalLRType;
  FilterMap["BIOFILTER"] = BioFilterType;
  FilterMap["MARKERGENOEFF"]=MarkerGenoType;
  FilterMap["ALLELEFREQ"] = AlleleFreqFilterType;
  FilterMap["MDRPDT"] = MDRPDTFilterType;
  FilterMap["MAXSCORE"] = MaxScoreFilterType;
  FilterMap["PERMSCORE"] = PermScoreFilterType;
#ifdef HAVE_R
  FilterMap["MARS"] = MarsFilterType;
#endif
}
}


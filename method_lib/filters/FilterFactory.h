//FilterFactory.h
#ifndef __FILTER_FACTORY_H__
#define __FILTER_FACTORY_H__

#include <map>
#include "FilterExcept.h"

namespace Filters{

class Filter;

/// creates filters
class FilterFactory{
  public:

    /// Creates filter and returns.  Throws exception if
    /// no filter matches that name
    static Filter * create_filter(std::string filterName);

    /// Enumeration for filter types with a new type added for each new filter in system
    enum FilterType{
      /// Enum for missing
      MissingFilterType,
      /// Enum for Chi Square filter analysis of genotypes
      ChiSquareFilterType,
      /// Enum for Sequential Replication filter
      SeqRepFilterType,
      /// Enum for MDR filter
      MDRFilterType,
      /// Enum for Logistic Regression filter
      LogRegressFilterType,
      /// Enum for TAILS filter
      TAILSFilterType,
      /// Enum for Chi Square on alleles
      ChiAlleleFilterType,
      /// Enum for Allele frequency calculator
      FreqCountFilterType,
      /// Enum for NMI filter
      NMIFilterType,
      /// Enum for Contingency Table Filter (used for Will's study)
      ContingencyFilterType,
      /// Enum for Armitage Filter
      ArmitageFilterType,
      /// Enum for Odds Ratio Filter
      OddsRatioFilterType,
      /// Enum for Likelihood Ratio Filter
      LikelihoodFilterType,
      /// Enum for Uncertainty Coefficient Filter
      UncertaintyCoeffFilterType,
      /// Enum for Conditional Logistic Regression Filter
      ConditionalLRType,
      /// Enum For BioFilter
      BioFilterType,
      /// Enum For Marker genotype efficiency filter
      MarkerGenoType,
      /// Enum For Allele Frequency filter
      AlleleFreqFilterType,
      /// Enum for MDRPDT Filter
      MDRPDTFilterType,
      /// Enum for Max Score Filter
      MaxScoreFilterType,
      /// Enum for Perm Score Filter
      PermScoreFilterType,
      /// Enum for Earth/Mars Score Filter
      MarsFilterType
    };

  private:

    /// Sets map for use in Filter creation
    static void SetFilterMap();

    static std::map<std::string, FilterType> FilterMap;

};
}

#endif


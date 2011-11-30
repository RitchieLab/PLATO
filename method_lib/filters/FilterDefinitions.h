//FilterDefinitions.h

#ifndef __FILTERDEFINITIONS_H__
#define __FILTERDEFINITIONS_H__

namespace Filters{

typedef map<string, string> PARAMS;

struct FilterParams{
  string filter_name;
  PARAMS filter_params; 
  vector<FilterParams> sub_filters;
};

}

#endif

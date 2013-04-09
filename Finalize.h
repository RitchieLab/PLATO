#ifndef FINALIZE_H
#define FINALIZE_H

#include <vector>
#include "Marker.h"
#include "Sample.h"
#include "Family.h"

class Finalize{
	private:
		
	public:
		Finalize(){
		};
		virtual ~Finalize(){};
		void finish(std::vector<Methods::Marker*>*, std::vector<Methods::Sample*>*, std::vector<Methods::Family*>*);
};
#endif

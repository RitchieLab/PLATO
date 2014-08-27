/*
 * OutputPLINK.cpp
 *
 *  Created on: Dec 3, 2013
 *      Author: jrw32
 */

#include "OutputPLINK.h"

#include "data/Marker.h"
#include "data/Sample.h"

#include "InputManager.h"

using PLATO::Data::Marker;
using PLATO::Data::Sample;

namespace PLATO {
namespace Utility {

void OutputPLINK::printPEDHeader(std::ostream& ped_f, const Sample& s) const {
	ped_f << s.getFID() << "\t" << s.getID() << "\t"
			<< (s.getFather() == 0 ? "0" : s.getFather()->getID()) << "\t"
			<< (s.getMother() == 0 ? "0" : s.getMother()->getID()) << "\t"
			<< (s.isGenderKnown() ? (s.isFemale() + 1) : 0) << "\t"
			<< (s.isAffectedKnown() ? (s.isAffected() + 1) : -9);
}

void OutputPLINK::printMAPInfo(std::ostream& map_f, const Marker& m, bool print_alleles) const {
	map_f << InputManager::chrIntToString(m.getChrom()) << "\t"
		  << m.getID() << "\t" << 0 << "\t"  << m.getLoc();

	if(print_alleles){
		map_f << "\t" << m.getRefAllele() << "\t" << m.getAltAllele();
	}

}

}
}

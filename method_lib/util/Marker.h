#ifndef METHODS_MARKER_H
#define METHODS_MARKER_H

#include <string>
#include <vector>

namespace Methods{

class Sample;

class Marker{

public:
	Marker(const std::string& chrom, unsigned int loc, const std::string& id, unsigned int idx);

	/*!
	 * Adds an allele to the Marker.  Returns the index of the given allele.
	 * (Note: runs in O(n) time)
	 */
	unsigned int addAllele(const std::string& allele);
	/*!
	 * Sets the given allele as the referent allele.  Returns false if the given
	 * allele is not found in the Marker
	 */
	bool setRefAllele(const std::string& allele);

	/*!
	 * Sets the given allele as the alternate allele (in biallelic analysis).
	 * Returns false if only one allele is available or if the given allele is
	 * not found in the allele list.
	 */
	bool setAltAllele(const std::string& allele);

	bool setMAF(float maf){return -1 != (_maf = (maf >= 0 && maf <= 1) ? maf : -1);}

	void setEnabled(bool enabled=true) {_chr = getChrom() | (ENABLED_MASK * enabled); }

	bool issetMAF() const {return _maf >= 0 && _maf <= 1;};
	float getMAF() const {return _maf;}
	/*!
	 * Returns the reference allele, or "0" if no reference allele is available.
	 * Note that the first allele added is by default the reference allele.
	 */
	const std::string& getRefAllele() const{
		return _ref_idx == static_cast<unsigned char>(-1) ? _missing_allele : _alleles[_ref_idx];
	}

	/*!
	 * Returns the alternate allele; the first available allele (as added) which
	 * is not the reference allele.  Returns "0" if no available allele
	 */
	const std::string& getAltAllele() const{
		return _alt_idx == static_cast<unsigned char>(-1) ? _missing_allele : _alleles[_alt_idx];
	}

	const std::string& getAllele(unsigned char idx) const{
		return idx >= _alleles.size() ? _missing_allele : _alleles[idx];
	}

	unsigned char getRefIdx() const{
		return _ref_idx;
	}

	unsigned char getAltIdx() const{
		return _alt_idx;
	}

	/*!
	 * Returns the index into the chromosome array
	 */
	unsigned int getIndex() const {return _idx;}

	/*!
	 * Returns the chromosome
	 */
	unsigned short getChrom() const {return _chr & ~(ENABLED_MASK);}

	/*!
	 * Returns the location
	 */
	unsigned int getLoc() const {return _loc;}

	/*!
	 * Returns the ID (usually an RSID)
	 */
	const std::string& getID() const {return _id;}

	/*!
	 * Sort by chromosome, then position
	 */
	bool operator<(const Marker& o) const {
		return (getChrom() == o.getChrom()) ? (_loc < o._loc) : (_chr < o._chr);
	}

	bool isEnabled() const {return _chr & ENABLED_MASK;	}

	static const std::string& getMissingAllele() {return _missing_allele;}



private:


	// NOTE: the most significant bit of _chr is the enabled bit
	unsigned short _chr;
	unsigned char _ref_idx;
	unsigned char _alt_idx;
	unsigned int _loc;
	// NOTE: this must NEVER change!! It is the index into the genotype memory array!
	const unsigned int _idx;
	float _maf;
	std::string _id;
	std::vector<std::string> _alleles;

	static std::string _missing_allele;
	static const unsigned short ENABLED_MASK = 1 << (sizeof(unsigned short)*8 - 1);
};

}

// define an ordering for Marker pointers
namespace std{
template<>
struct less<Methods::Marker*> {
	bool operator()(const Methods::Marker* x,
			const Methods::Marker* y) const {
		return (y != 0 && x != 0) ? (*x) < (*y) : y < x;
	}
};
}

#endif

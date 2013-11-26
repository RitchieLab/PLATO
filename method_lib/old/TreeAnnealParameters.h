/*
 * TreeAnnealParameters.h
 *
 *  Created on: Nov 3, 2009
 *      Author: gilesjt
 */

#ifndef TREEANNEALPARAMETERS_H_
#define TREEANNEALPARAMETERS_H_

namespace Methods {

class TreeAnnealParameters {
public:
	TreeAnnealParameters();
	TreeAnnealParameters(int t, int o, int mm, int n);
	virtual ~TreeAnnealParameters();

	void setTreeSize(int s){treesize = s;}
	int getTreeSize(){return treesize;}
	void setOpers(int s){opers = s;}
	int getOpers(){return opers;}
	void setMinMass(int s){minmass = s;}
	int getMinMass(){return minmass;}
	void setSampSize(int s){n1 = s;}
	int getSampSize(){return n1;}

private:
	int treesize;  //max number of allowed leaves per logic tree
	int opers; //1 = and + or, 2 = and, 3 = or
	int minmass; //min number of cases for ay tree needs to be 1 and for which any tree needs to be 0 to be considered logic tree in model
	int n1; //sampsize - used to check that minmass is smaller than n1/4
};

}

#endif /* TREEANNEALPARAMETERS_H_ */

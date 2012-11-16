/*
 * AnnealParameters.h
 *
 *  Created on: Nov 3, 2009
 *      Author: gilesjt
 */

#ifndef ANNEALPARAMETERS_H_
#define ANNEALPARAMETERS_H_

namespace Methods {

class AnnealParameters {
public:
	AnnealParameters();
	AnnealParameters(int s, int e, int i, bool e, int u);
	virtual ~AnnealParameters();

	void setStart(int s){start = s;}
	int getStart(){return start;}

	void setEnd(int s){end = s;}
	int getEnd(){return end;}

	void setIter(int s){iter = s;}
	int getIter(){return iter;}

	void setEarlyOut(bool s){earlyout = s;}
	bool getEarlyOut(){return earlyout;}

	void setUpdate(int s){update = s;}
	int getUpdate(){return update;}


private:
	int start;  //upper temp on log10 scale
	int end;    //lower temp on log10 scale
	int iter;   //num iters on annealing chain
	bool earlyout; //stop if hung

	int update; //stdout updating frequency

};

}

#endif /* ANNEALPARAMETERS_H_ */

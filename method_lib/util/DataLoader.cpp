/*
 * DataLoader.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "DataLoader.h"

#include <stdexcept>

using std::string;

void DataLoader::readPed(const string& fn){
	map<string, vector<string> > descinfo;
	vector<string> exclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> fexclude;
	vector<string> finclude;
	vector<string> descheaders;

    FILE* input;
    input = fopen(fn.c_str(), "r");
	if(!input){
		throw std::invalid_argument("Error opening pedfile: " + fn + ".");
	}

	int onind = -1;
    while(!feof(input)){
		onind++;
        Methods::Sample* samp = new Sample();
        int f = 0;
        string temp = "";
        if(readString(input, &temp)){
           samp->setFamID(temp);
              f++;
             temp = "";
         }

		string ftemp = samp->getFamID();
		if(samp->getFamID() == ""){
			delete(samp);
			continue;
		}
        if(ftemp.at(0) == '#'){
			delete(samp);
			while(fgetc(input) != '\n' && !feof(input)){}
            continue;
        }
        string sex = "";
        string pheno = "";
        if(readString(input, &temp)){
            samp->setInd(temp);
            f++;
            temp = "";
        }
		samp->setEnabled(true);

		if(exclude.size() > 0){
			vector<string>::iterator found = find(exclude.begin(), exclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != exclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
        if(readString(input, &temp)){
            samp->setDadID(temp);
            f++;
            temp = "";
        }
        if(readString(input, &temp)){
            samp->setMomID(temp);
            f++;
            temp = "";
        }
        if(readString(input, &sex)) f++;
        if(readString(input, &pheno)) f++;

        if(sex == "1"){
            samp->setSex(true);
        }
        else if(sex == "2"){
            samp->setSex(false);
        }

		if(pheno == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		float mypheno = atof(pheno.c_str());
		samp->setPheno(mypheno);
        if(opts::pedinfo.size() > 0){
            map<string, Methods::Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
            if(sfind != opts::pedinfo.end()){
				Methods::Sample* sfound = sfind->second;
                samp->setDadID(sfound->getDadID());
                samp->setMomID(sfound->getMomID());
                samp->setSex(sfound->getSex());
                samp->setPheno(sfound->getPheno());
            }
        }
		if(samp->getPheno() != 0.0f && samp->getPheno() != 1.0f && samp->getPheno() != 2.0f){
			opts::_BINTRAIT_ = false;
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders.at(i), tokens.at(i));
				}
				else{
					samp->assignDetail(descheaders.at(i), "NA");
				}
			}

		}
        samp->resizeAlleles(markers->size());

        int gn = 0;
        int i = 0;
        bool linedone = false;
        //bool fatal = false;

        string fmsg;
        while(!linedone){
            string one = "";
            string two = "";

            while(1)
            {
                char ch = fgetc(input);

                if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input))
                {
                    if(ch == '\n' || ch == '\r' || feof(input))
                    {
                        linedone = true;
                    }

                    if(one.length() > 0)
                    {
                        gn++;
                        break;
                    }
                    if(ch == '\n' || ch == '\r' || feof(input))
                    {
                        break;
                    }
                }
                else
                {
                    one += ch;
                }
            }
            if(!linedone)
            {
                while(1)
                {
                    char ch = fgetc(input);
                    if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input))
                    {
                        if(ch == '\n' || ch == '\r' || feof(input))
                        {
                            linedone = true;
                        }
                        if(two.length() > 0)
                        {
                            gn++;
                            break;
                        }
                        if(ch == '\n' || ch == '\r' || feof(input))
                        {
                            break;
                        }
                    }
                    else
                    {
                        two += ch;
                    }
                }
                if(linedone && one.length() == 0 && two.length() == 0)
                {
                    break;
                }

				if(i > (int)markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
				Methods::Marker* m = (*markers)[(*marker_map).at(i)];
				if(m->isEnabled()){
					int oldallelecount = m->getNumAlleles();
	                if(one != "0"){
						//new
						if(m->getAlleleLoc(one) < 0){
							m->addAllele(one);
						}

        	        }

            	    if(two != one){
                	    if(two != "0"){
							//new
							if(m->getAlleleLoc(two) < 0){
								m->addAllele(two);
							}

                    	}
	                }


					if(m->getNumAlleles() <= 2){
   		 	            if(one == m->getAllele1() && two == m->getAllele1()){
   	    	     	        samp->addAone(i, false);
   	     		            samp->addAtwo(i, false);
		                }
	    	            else if(one != "0" && two != "0" && one != two){
	        	            samp->addAone(i, false);
 		           	        samp->addAtwo(i, true);
	                	}
		                else if(one == m->getAllele2() && two == m->getAllele2()){
	    	                samp->addAone(i, true);
	        	            samp->addAtwo(i, true);
	            	    }
	                	else if(one == "0" || two == "0"){
		                    samp->addAone(i, true);
	    	                samp->addAtwo(i, false);
	        	        }
					}
					else if(opts::_MICROSATS_){
						samp->addMicroSat(i);
						int loc1 = m->getAlleleLoc(one);
						int loc2 = m->getAlleleLoc(two);

						samp->addAbone(i, loc1);
						samp->addAbtwo(i, loc2);
						if(oldallelecount <= 2){
							remapSamples(samples, markers, marker_map, i);
						}
					}
					else if(m->getNumAlleles() > 2 && !opts::_MICROSATS_){
						opts::printLog("More than 2 unique alleles found for map location: " + getString<int>(i) + ", line: " + getString<int>(onind + 1) + ".  Microsatellites not specified.\n");
						throw MethodException("More than 2 unique alleles found for map location: " + getString<int>(i) + ", line: " + getString<int>(onind + 1) + ".  Microsatellites not specified.\n");
					}
				}
				else{
					samp->addAone(i,true);
					samp->addAtwo(i, false);
				}
    	    	i++;
				if(i > (int)markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
            }/*end !linedone*/
        }/*end while(1)*/
		if(gn != (int)(2* markers->size())){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(((2 * markers->size()) + 6));
					text += " columns but found ";
					text += getString<int>((f + gn));
					text += "\n";
					throw MethodException(text);
		}

        vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Methods::Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				samp->setEnabled(false);
				vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				samp->setEnabled(false);
				vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
    }/*end while(eof)*/

    fclose(input);
}


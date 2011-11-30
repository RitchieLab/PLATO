/*
 * OptionsList.cpp
 *
 *  Created on: Jan 12, 2009
 *      Author: gilesjt
 */
#include <Helper.h>
#include "OptionsList.h"

OptionsList::OptionsList(string f) {
	design_file = f;

	try{
		optionsXml = Gnome::Glade::Xml::create(Vars::PLATO_OPTIONS_GLADE);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	ifstream input;
	input.open("resources/locus-filter.txt", ios::in);//Vars::PLATO_DESIGN_FILES[design_file].c_str(), ios::in);

	if(!input){
		//throw MethodException("Internal Error: Cannot find resource file - " + Vars::PLATO_DESIGN_FILES[design_file] + "\n");
	}

	int count = 1;
	string line = "";
	getline(input, line);
	if(line != design_file){
		input.close();
		//throw MethodException("Internal Error: Design file - " + Vars::PLATO_DESIGN_FILES[design_file] + " doesn't match request for '" + design_file + "'\n");
	}
	while(getline(input, line)){
		count++;
		if(line.size() == 0 || (line.size() > 0 && line.at(0) == '#')){
			continue;
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() == 0){
			continue;
		}

		string option = tokens[0];
		if(option == Vars::DISABLE_TOKEN || option == Vars::REQUIRE_TOKEN){
			//throw MethodException("Internal Error: Design file - " + Vars::PLATO_DESIGN_FILES[design_file] + " incorrect format on line " + getString<int>(count) + "\n");
		}

		OptionData* opt = new OptionData();
		opt->set_name(option);

		option = option.substr(1,option.size() - 1);
		for(int c = 0; c < option.size(); c++){
			if(option[c] == '-'){
				option[c] = '_';
			}
		}
		Gtk::Layout* mylayout = 0;
		string widget_name = option + "_layout";
		optionsXml->get_widget(widget_name, mylayout);
		if(!mylayout){
			throw MethodException("Uhoh, layout is null!! - " + widget_name);
		}
		opt->set_layout(mylayout);

		for(int i = 1; i < tokens.size(); i++){
			if(tokens[i] == Vars::REQUIRE_TOKEN){
				if((i + 1) < tokens.size()){
					string reqs = tokens[++i];
					vector<string> subtokens;
					General::Tokenize(reqs, subtokens, ",");
					for(int st = 0; st < subtokens.size(); st++){
						opt->add_required(subtokens[st]);
					}
				}
			}
			else if(tokens[i] == Vars::DISABLE_TOKEN){
				if((i + 1) < tokens.size()){
					string dis = tokens[++i];
					vector<string> subtokens;
					General::Tokenize(dis, subtokens, ",");
					for(int st = 0; st < subtokens.size(); st++){
						opt->add_disabled(subtokens[st]);
					}
				}
			}
		}

		option_data.push_back(opt);
	}

	input.close();
}

OptionsList::~OptionsList() {
	for(int i = 0; i < option_data.size(); i++){
		delete(option_data[i]);
	}
	option_data.clear();
}


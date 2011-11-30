/*
 * YamlData.h
 *
 *  Created on: Jan 12, 2009
 *      Author: gilesjt
 *
 *Used as parent class for parsing yaml file
 */

#ifndef YAMLDATA_H_
#define YAMLDATA_H_
#include <string>
#include <string.h>
#include <yaml.h>
#include <iostream>
#include <stdio.h>
#include <MethodException.h>

using namespace std;

class YamlData {
public:
	YamlData(string f);
	YamlData(){};

	virtual ~YamlData();

	void parse();

	virtual void yaml_stream_start_event() = 0;
	virtual void yaml_stream_end_event() = 0;
	virtual void yaml_document_start_event() = 0;
	virtual void yaml_document_end_event() = 0;
	virtual void yaml_alias_event() = 0;
	virtual void yaml_scalar_event() = 0;
	virtual void yaml_sequence_start_event() = 0;
	virtual void yaml_sequence_end_event() = 0;
	virtual void yaml_mapping_start_event() = 0;
	virtual void yaml_mapping_end_event() = 0;

protected:
	string ymlfile;


};

#endif /* YAMLDATA_H_ */

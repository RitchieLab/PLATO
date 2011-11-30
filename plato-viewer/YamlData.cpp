/*
 * YamlData.cpp
 *
 *  Created on: Jan 12, 2009
 *      Author: gilesjt
 */

#include "YamlData.h"

YamlData::YamlData(string f) {
	ymlfile = f;
}

YamlData::~YamlData() {
	// TODO Auto-generated destructor stub
}

void YamlData::yaml_stream_start_event(){};
void YamlData::yaml_stream_end_event(){};
void YamlData::yaml_document_start_event(){};
void YamlData::yaml_document_end_event(){};
void YamlData::yaml_alias_event(){};
void YamlData::yaml_scalar_event(){};
void YamlData::yaml_sequence_start_event(){};
void YamlData::yaml_sequence_end_event(){};
void YamlData::yaml_mapping_start_event(){};
void YamlData::yaml_mapping_end_event(){};


void YamlData::parse(){
	yaml_parser_t parser;
	yaml_event_t input_event;
	FILE* file = fopen(ymlfile.c_str(), "rb");
	bool done = false;

	memset(&parser, 0, sizeof(parser));
	memset(&input_event, 0, sizeof(input_event));

	if(!yaml_parser_initialize(&parser)){
		throw MethodException("YAML - Could not initialize parser object\n");
	}

	yaml_parser_set_input_file(&parser, file);

	while(!done){
		int properties, key, value, map, seq;

		//get the next event
		if(!yaml_parser_parse(&parser, &input_event)){
			throw MethodException("YAML - Unknown exception parsing event!\n");
		}

		if(input_event.type == YAML_STREAM_END_EVENT){
			done = true;
		}

		switch(input_event.type){
			case YAML_STREAM_START_EVENT:
				cout << "in yaml_stream_start_event\n";
				yaml_stream_start_event();
				break;
			case YAML_STREAM_END_EVENT:
				cout << "in yaml_stream_end_event\n";
				break;
			case YAML_DOCUMENT_START_EVENT:
				cout << "in yaml_document_start_event\n";
				break;
			case YAML_DOCUMENT_END_EVENT:
				cout << "in yaml_document_end_event\n";
				break;
			case YAML_ALIAS_EVENT:
				cout << "in yaml_alias_event\n";
				break;
			case YAML_SCALAR_EVENT:
				cout << "in yaml_scalar_event\n";
				cout << "myscalar = " << input_event.data.scalar.value << endl;// "::" << input_event.data.scalar.value << endl;
				break;
			case YAML_SEQUENCE_START_EVENT:
				cout << "in yaml_sequence_start_event\n";
				break;
			case YAML_SEQUENCE_END_EVENT:
				cout << "in yaml_sequence_end_event\n";
				break;
			case YAML_MAPPING_START_EVENT:
				cout << "in yaml_mapping_start_event\n";
				break;
			case YAML_MAPPING_END_EVENT:
				cout << "in yaml_mapping_end_event\n";
				break;
			default:
				break;
		}

		yaml_event_delete(&input_event);
	}

	yaml_parser_delete(&parser);
	fclose(file);
}

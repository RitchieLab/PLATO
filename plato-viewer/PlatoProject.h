#ifndef PLATOPROJECT_H
#define PLATOPROJECT_H

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <string.h>
#include <gtkmm.h>
#include "DataSetObject.h"
#include "Process.h"
//#include "SingleModelColumns.h"

using namespace std;

class PlatoProject{

	public:
		PlatoProject();
		~PlatoProject(){};

		void set_project_name(string f){if(f.size() > 0) name = f;}
		string get_project_name(){return name;}
		void add_dataset(DataSetObject* dso){data_sets.push_back(dso);}
		vector<DataSetObject*> get_datasets(){return data_sets;}
		int num_datasets(){return (int)data_sets.size();}
///		void add_pedfile(string f){pedfiles.push_back(f);}
///		void add_mapfile(string f){mapfiles.push_back(f);}
///		void add_tpedfile(string f){tpedfiles.push_back(f);}
///		void add_tfamfile(string f){tfamfiles.push_back(f);}
		void set_pedmap(bool v){pedmap = v;}
		void set_tpedfam(bool v){tpedfam = v;}
		void set_binary(bool v){binary = v;}
		bool has_pedmap(){return pedmap;}
		bool has_tpedfam(){return tpedfam;}
		bool has_binary(){return binary;}
		vector<Process> get_batch(){return batch;}

///		vector<string> get_pedfiles(){return pedfiles;}
///		vector<string> get_mapfiles(){return mapfiles;}
///		vector<string> get_tpedfiles(){return tpedfiles;}
///		vector<string> get_tfamfiles(){return tfamfiles;}
		vector<string> get_outputfiles(){return outputfiles;}
//		Glib::RefPtr<Gtk::TreeStore> get_model(){return m_refTreeStore;}
		void reinitialize();
		void reset(){
//			for(int i = 0; i < data_sets.size(); i++){
//				delete(data_sets[i]);
//			}
			data_sets.clear();
///			pedfiles.clear();
///			mapfiles.clear();
			name = "My Project";
		}

		string get_dataset_summary();

	protected:
		string name;
///		vector<string> pedfiles;
///		vector<string> mapfiles;
///		vector<string> tpedfiles;
///		vector<string> tfamfiles;
		vector<string> outputfiles;
		vector<DataSetObject*> data_sets;
		vector<Process> batch;

		bool pedmap; //has ped/map input files
		bool tpedfam; //has tped/tfam input files
		bool binary; //has binary input files



};


#endif

#include "PlatoProject.h"

PlatoProject::PlatoProject(){
	name = "My Project";
	binary = pedmap = tpedfam = false;
}

void PlatoProject::reinitialize(){
/*	m_refTreeStore->clear();
	m_refTreeStore = Gtk::TreeStore::create(m_Columns);

	Gtk::TreeRow row = *(m_refTreeStore->append());


	row[m_Columns.m_col_value] = this->name;
*/
}

string PlatoProject::get_dataset_summary(){
	string summary = "\n\n";
	for(int i = 0; i < data_sets.size(); i++){
		DataSetObject* ds = data_sets[i];
		summary += ds->summarize();
		summary += "\n";
	}
	return summary;
}

/*
 * ProjectTree.cpp
 *
 *  Created on: Jan 19, 2009
 *      Author: gilesjt
 */

#include "ProjectTree.h"

ProjectTree::ProjectTree() {
	// TODO Auto-generated constructor stub
	TreeView();
	this->set_enable_tree_lines(true);

	Gtk::Menu::MenuList& menulist = m_Menu_Popup.items();
	menulist.push_back(Gtk::Menu_Helpers::MenuElem("_View", sigc::mem_fun(*this, &ProjectTree::on_menu_file_popup_generic)));
	menulist.push_back(Gtk::Menu_Helpers::MenuElem("_Remove", sigc::mem_fun(*this, &ProjectTree::on_menu_file_popup_generic)));
	m_Menu_Popup.accelerate(*this);
}

ProjectTree::~ProjectTree() {
	// TODO Auto-generated destructor stub
}

bool ProjectTree::on_button_press_event(GdkEventButton* event){
	bool return_value = TreeView::on_button_press_event(event);
	if((event->type == GDK_BUTTON_PRESS) && (event->button == 3)){
		m_Menu_Popup.popup(event->button, event->time);
	}
	return return_value;
}

void ProjectTree::on_menu_file_popup_generic(){
	Glib::RefPtr<Gtk::TreeView::Selection> refSelection = get_selection();
	if(refSelection){
		Gtk::TreeModel::iterator iter = refSelection->get_selected();
		if(iter){
			Glib::ustring id = (*iter)[m_columns.m_col_value];
			cout << " selected id: " << id << endl;
		}
	}
}

void ProjectTree::generate(PlatoProject* proj){
	this->remove_all_columns();
	m_refTreeModel = Gtk::TreeStore::create(m_columns);

	Gtk::TreeRow row = *(m_refTreeModel->append());
	row[m_columns.m_col_value] = proj->get_project_name();

	Gtk::TreeRow child_row = *(m_refTreeModel->append(row.children()));
	child_row[m_columns.m_col_value] = "Input Files";

//	vector<string> ped_files = proj->get_pedfiles();
//	vector<string> map_files = proj->get_mapfiles();
//	vector<string> tped_files = proj->get_tpedfiles();
//	vector<string> tfam_files = proj->get_tfamfiles();

	vector<DataSetObject*> data_sets = proj->get_datasets();
	int id_count = 0;

	if(data_sets.size() == 0){//ped_files.size() == 0 && tped_files.size() == 0){
		Gtk::TreeRow input_files = *(m_refTreeModel->append(child_row.children()));
		input_files[m_columns.m_col_value] = "No input files";
		input_files[m_columns.m_col_id] = -1;
	}
	else{
		for(int i = 0; i < data_sets.size(); i++){
			vector<string> files = data_sets[i]->get_files();
			for(int j = 0; j < files.size(); j++){
				Gtk::TreeRow input_files = *(m_refTreeModel->append(child_row.children()));
				input_files[m_columns.m_col_value] = files[j];
				input_files[m_columns.m_col_id] = id_count++;
			}
		}

/*		for(int i = 0; i < ped_files.size(); i++){
			Gtk::TreeRow input_files = *(m_refTreeModel->append(child_row.children()));
			input_files[m_columns.m_col_value] = ped_files[i];
			input_files[m_columns.m_col_id] = id_count++;

			Gtk::TreeRow input_files2 = *(m_refTreeModel->append(child_row.children()));
			input_files2[m_columns.m_col_value] = map_files[i];
			input_files2[m_columns.m_col_id] = id_count++;
		}
		for(int i = 0; i < tped_files.size(); i++){
			Gtk::TreeRow input_files = *(m_refTreeModel->append(child_row.children()));
			input_files[m_columns.m_col_value] = tped_files[i];
			input_files[m_columns.m_col_id] = id_count++;

			Gtk::TreeRow input_files2 = *(m_refTreeModel->append(child_row.children()));
			input_files2[m_columns.m_col_value] = tfam_files[i];
			input_files2[m_columns.m_col_id] = id_count++;
		}
*/
	}

	Gtk::TreeRow child_row2 = *(m_refTreeModel->append(row.children()));
	child_row2[m_columns.m_col_value] = "Output Files";

	vector<string> output_files = proj->get_outputfiles();
	if(output_files.size() == 0){
		Gtk::TreeRow output_files = *(m_refTreeModel->append(child_row2.children()));
		output_files[m_columns.m_col_value] = "No output files";
		output_files[m_columns.m_col_id] = -1;
	}
	else{
		for(int i = 0; i < output_files.size(); i++){
			Gtk::TreeRow outputs = *(m_refTreeModel->append(child_row2.children()));
			outputs[m_columns.m_col_value] = output_files[i];
			outputs[m_columns.m_col_id] = id_count++;
		}
	}

	this->set_model(m_refTreeModel);
	this->append_column("Project", m_columns.m_col_value);
	Gtk::TreeViewColumn* pColumn = get_column(0);
	if(pColumn){

	}

}

/*
 * FileTypeModel.h
 *
 *  Created on: Jan 6, 2009
 *      Author: gilesjt
 */

#ifndef FILETYPEMODEL_H_
#define FILETYPEMODEL_H_
#include <gtkmm.h>
#include "DataSetObject.h"
#include "PlatoProject.h"


using namespace std;

class FileTypeModel {
public:
	FileTypeModel();
	virtual ~FileTypeModel();

	Gtk::Widget* get_widget(){return panes;}
	void assign_results(PlatoProject*);

protected:
	virtual void on_file_type_changed();

	Gtk::VPaned* panes;
	Gtk::Alignment* type_alignment;
	Gtk::Layout* type_layout;
	Gtk::Alignment* form_alignment;
	Gtk::Table* form_table;
	Gtk::Layout* form_layout;
	vector<Gtk::Widget*> form_widgets;
	vector<Gtk::FileChooserButton*> file_choosers;

	Gtk::Label* type_label;
	Gtk::ComboBox* file_type;
	Glib::RefPtr<Gtk::ListStore> m_refTreeModel;

	class ModelColumns : public Gtk::TreeModel::ColumnRecord{
	public:
		ModelColumns(){
			add(m_col_id); add(m_col_name);
		}
		Gtk::TreeModelColumn<int> m_col_id;
		Gtk::TreeModelColumn<Glib::ustring> m_col_name;
	};

	ModelColumns m_Columns;

};


#endif /* FILETYPEMODEL_H_ */

/*
 * FileTypeModel.cpp
 *
 *  Created on: Jan 6, 2009
 *      Author: gilesjt
 */

#include "FileTypeModel.h"

FileTypeModel::FileTypeModel() {
	panes = new Gtk::VPaned();
	panes->set_position(30);
	type_alignment = new Gtk::Alignment();
	form_alignment = new Gtk::Alignment();
	form_layout = new Gtk::Layout();
	type_layout = new Gtk::Layout();
	type_label = new Gtk::Label();
	file_type = new Gtk::ComboBox();

	type_alignment->set_size_request(-1,-1);
	type_alignment->set_size_request(-1, 35);

	form_alignment->set_size_request(-1,-1);
	type_alignment->add(*type_layout);
	form_alignment->add(*form_layout);
	panes->add1(*type_alignment);
	panes->add2(*form_alignment);

	panes->set_size_request(400,400);

	form_table = new Gtk::Table(2,2,true);
	form_table->set_border_width(10);
	form_table->set_row_spacings(5);
	form_layout->add(*form_table);

	type_label->set_text("File Type:");
	type_label->set_size_request(-1,25);
	type_label->set_alignment(0,5);
	type_layout->add(*type_label);

	int count = 0;
	m_refTreeModel = Gtk::ListStore::create(m_Columns);
	file_type->set_model(m_refTreeModel);
	Gtk::TreeModel::Row row = *(m_refTreeModel->append());
	row[m_Columns.m_col_id] = count++;
	row[m_Columns.m_col_name] = "PED/Map Files";

	row = *(m_refTreeModel->append());
	row[m_Columns.m_col_id] = count++;
	row[m_Columns.m_col_name] = "Transposed Ped & Map Files";

	row = *(m_refTreeModel->append());
	row[m_Columns.m_col_id] = count++;
	row[m_Columns.m_col_name] = "Plink Binary Input Files";

	file_type->pack_start(m_Columns.m_col_name);

	file_type->signal_changed().connect(sigc::mem_fun(*this, &FileTypeModel::on_file_type_changed));

	file_type->set_size_request(-1,25);

	type_layout->put(*file_type, 80, 5);
	panes->show_all_children();
}

FileTypeModel::~FileTypeModel() {
	for(int i = 0; i < form_widgets.size(); i++){
		delete form_widgets[i];
	}
	delete(form_table);
	delete(type_label);
	delete(file_type);

	delete(type_alignment);
	delete(type_layout);
	delete(panes);
}

void FileTypeModel::assign_results(PlatoProject* proj){
	Gtk::TreeModel::iterator iter = file_type->get_active();
	if(iter){
		Gtk::TreeModel::Row row = *iter;

		switch(row[m_Columns.m_col_id]){
		case 0:{
			DataSetObject* dso = new DataSetObject();
			dso->set_pedfile(file_choosers[0]->get_filename());
			dso->set_mapfile(file_choosers[1]->get_filename());
			dso->set_type(PEDMAP);
			proj->add_dataset(dso);
			proj->set_pedmap(true);
			break;
		}
		case 1:{
			DataSetObject* dso = new DataSetObject();
			dso->set_tpedfile(file_choosers[0]->get_filename());
			dso->set_tfamfile(file_choosers[1]->get_filename());
			dso->set_type(TPEDFAM);
			proj->add_dataset(dso);
			proj->set_tpedfam(true);
			break;
		}
		case 2:{
			proj->set_binary(true);
			break;
		}
		}
	}
}

void FileTypeModel::on_file_type_changed(){
	for(int c = 0; c < form_widgets.size(); c++){
		form_table->remove(*(form_widgets[c]));
		delete(form_widgets[c]);
	}
	form_widgets.clear();
	file_choosers.clear();

	Gtk::TreeModel::iterator iter = file_type->get_active();
	if(iter){
		Gtk::TreeModel::Row row = *iter;

		switch(row[m_Columns.m_col_id]){
		case 0:{
			Gtk::Label* new_label = new Gtk::Label();
			new_label->set_text("PED File:");
			new_label->set_size_request(-1,25);
			form_widgets.push_back(new_label);
			form_table->attach(*new_label, 0, 1, 0, 1);

			Gtk::FileChooserButton* ped_fcb = new Gtk::FileChooserButton();
			ped_fcb->set_size_request(100,25);
			ped_fcb->set_name("ped_file_chooser");
			form_widgets.push_back(ped_fcb);
			form_table->attach(*ped_fcb, 1, 2, 0, 1);
			file_choosers.push_back(ped_fcb);

			Gtk::Label* new_label2 = new Gtk::Label();
			new_label2->set_text("Map File:");
			new_label2->set_size_request(-1,25);
			form_widgets.push_back(new_label2);
			form_table->attach(*new_label2, 0, 1, 1, 2);

			Gtk::FileChooserButton* map_fcb = new Gtk::FileChooserButton();
			map_fcb->set_size_request(100,25);
			map_fcb->set_name("map_file_chooser");
			form_widgets.push_back(map_fcb);
			form_table->attach(*map_fcb, 1, 2, 1, 2);
			file_choosers.push_back(map_fcb);

			form_layout->show_all_children(true);
			break;
		}
		case 1:
		{
			Gtk::Label* new_label = new Gtk::Label();
			new_label->set_text("TPED File:");
			new_label->set_size_request(-1,25);
			form_widgets.push_back(new_label);
			form_table->attach(*new_label, 0, 1, 0, 1);

			Gtk::FileChooserButton* ped_fcb = new Gtk::FileChooserButton();
			ped_fcb->set_size_request(100,25);
			ped_fcb->set_name("tped_file_chooser");
			form_widgets.push_back(ped_fcb);
			form_table->attach(*ped_fcb, 1, 2, 0, 1);
			file_choosers.push_back(ped_fcb);

			Gtk::Label* new_label2 = new Gtk::Label();
			new_label2->set_text("Family Map File:");
			new_label2->set_size_request(-1,25);
			form_widgets.push_back(new_label2);
			form_table->attach(*new_label2, 0, 1, 1, 2);

			Gtk::FileChooserButton* map_fcb = new Gtk::FileChooserButton();
			map_fcb->set_size_request(100,25);
			map_fcb->set_name("tfam_file_chooser");
			form_widgets.push_back(map_fcb);
			form_table->attach(*map_fcb, 1, 2, 1, 2);
			file_choosers.push_back(map_fcb);

			form_layout->show_all_children(true);

			break;
		}
		case 2:
		{
			Gtk::Label* new_label = new Gtk::Label();
			new_label->set_text("Binary Map File:");
			new_label->set_size_request(-1,25);
			form_widgets.push_back(new_label);
			form_table->attach(*new_label, 0, 1, 0, 1);

			Gtk::FileChooserButton* ped_fcb = new Gtk::FileChooserButton();
			ped_fcb->set_size_request(100,25);
			ped_fcb->set_name("bfile_file_chooser");
			form_widgets.push_back(ped_fcb);
			form_table->attach(*ped_fcb, 1, 2, 0, 1);
			file_choosers.push_back(ped_fcb);

			form_layout->show_all_children(true);
			break;
		}
		case 3:
			break;
		}
	}

}

/*
 * BatchDialog.cpp
 *
 *  Created on: Jan 22, 2009
 *      Author: gilesjt
 */

#include <iostream>
#include <map>
#include "BatchDialog.h"

BatchDialog::BatchDialog() {
	// TODO Auto-generated constructor stub
	Dialog();
}

BatchDialog::BatchDialog(BaseObjectType* cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade)
: Gtk::Dialog(cobject)
{
	summary_textview = 0;
	options_notebook = 0;
	plus_button = 0;
	minus_button = 0;
	ok_button = 0;
	cancel_button = 0;
	process_chooser = 0;
	process_chooser_ok = 0;
	process_chooser_cancel = 0;
	batch_treeview = 0;
	process_chooser_vbox = 0;

	refGlade->get_widget("batch_summary_textview", summary_textview);
	refGlade->get_widget("batch_settings_notebook", options_notebook);
	refGlade->get_widget("batch_plus", plus_button);
	refGlade->get_widget("batch_minus", minus_button);
	refGlade->get_widget("batch_dialog_ok", ok_button);
	refGlade->get_widget("batch_dialog_cancel", cancel_button);
	refGlade->get_widget("process_chooser_dialog", process_chooser);
	refGlade->get_widget("process_chooser_ok", process_chooser_ok);
	refGlade->get_widget("process_chooser_cancel", process_chooser_cancel);
	refGlade->get_widget("process_chooser_vbox", process_chooser_vbox);
	refGlade->get_widget("batch_treeview", batch_treeview);

	color.set_rgb_p(0.863, 0.855, 0.835);

	if(summary_textview){
		summary_textview->modify_base(Gtk::STATE_NORMAL, color);
	}

	if(batch_treeview){
		m_refTreeModel = Gtk::ListStore::create(m_columns);

		batch_treeview->set_model(m_refTreeModel);
		batch_treeview->append_column("Process", m_columns.m_col_value);
		Glib::RefPtr<Gtk::TreeSelection> selection = batch_treeview->get_selection();
		selection->signal_changed().connect(sigc::mem_fun(*this, &BatchDialog::on_batch_treeview_changed));
	}

	if(plus_button){
		plus_button->signal_clicked().connect(sigc::mem_fun(*this, &BatchDialog::on_plus_button_clicked));
	}

	if(minus_button){
		minus_button->signal_clicked().connect(sigc::mem_fun(*this, &BatchDialog::on_minus_button_clicked));
	}

	if(process_chooser_ok){
		process_chooser_ok->signal_clicked().connect(sigc::mem_fun(*this, &BatchDialog::on_process_chooser_ok_clicked));
	}

	if(process_chooser_cancel){
		process_chooser_cancel->signal_clicked().connect(sigc::mem_fun(*this, &BatchDialog::on_process_chooser_cancel_clicked));
	}

	if(ok_button){
		ok_button->signal_clicked().connect(sigc::mem_fun(*this, &BatchDialog::on_ok_button_clicked));
		//put signal to parent here
	}
	if(cancel_button){
		cancel_button->signal_clicked().connect(sigc::mem_fun(*this, &BatchDialog::on_cancel_button_clicked));
	}


	{
		process_chooser_combobox = new Gtk::ComboBoxText();

		std::map<string, int>::iterator iter;
		for(iter = Vars::processes.begin(); iter != Vars::processes.end(); iter++){
			process_chooser_combobox->append_text(iter->first);
		}
		process_chooser_vbox->add(*process_chooser_combobox);
		process_chooser_vbox->show_all();
	}

	initialize();
}


BatchDialog::~BatchDialog() {
	// TODO Auto-generated destructor stub
}

void BatchDialog::initialize(){
	m_refTreeModel->clear();
	while(options_notebook->get_n_pages() > 1){
		options_notebook->remove_page(1);
	}
	for(int i = 0; i < batch_processes.size(); i++){

	}
}


void BatchDialog::on_plus_button_clicked(){

	process_chooser->show();

}

void BatchDialog::on_minus_button_clicked(){
	Glib::RefPtr<Gtk::TreeSelection> refTreeSelection = batch_treeview->get_selection();
//	refTreeSelection->selected_foreach_iter(sigc::mem_fun(*this, &AmdGui::selected_row_remove_callback));

	Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();
	if(iter){
		Gtk::TreeModel::Row row = *iter;
		ostringstream ss;
		ss << row[m_columns.m_col_value];
		refTreeSelection->unselect_all();
		m_refTreeModel->erase(row);
		for(int i = 0; i < batch_processes.size(); i++){
			if(batch_processes[i].get_name() == ss.str()){
				batch_processes.erase(batch_processes.begin() + i);
				cout << "Erased process\n";
			}
		}
	}
	display_process("");

}

void BatchDialog::on_process_chooser_ok_clicked(){
	process_chooser->hide();

	Gtk::TreeModel::Row row = addProcess(process_chooser_combobox->get_active_text());

	Glib::RefPtr<Gtk::TreeSelection> sel = batch_treeview->get_selection();
	sel->select(row);

//	batch_treeview->append_column("Process", m_columns.m_col_value);

}

Gtk::TreeModel::Row BatchDialog::addProcess(string s){
	Gtk::TreeModel::Row row = *(m_refTreeModel->append());
	row[m_columns.m_col_value] = s;

//add code to add process to batch
	switch(	Vars::processes[p]){
	case Processes::p_markerefficiency:
	{
		MarkerEfficiencyProcess* mep = new MarkerEfficiencyProcess();
		batch_processes.push_back(*mep);
		break;
	}
	case Processes::p_sampleefficiency:
	{
		break;
	}
	case -1:
		break;
	default:
		break;
	}

	return row;
}

void BatchDialog::on_batch_treeview_changed(){
	Glib::RefPtr<Gtk::TreeSelection> refTreeSelection = batch_treeview->get_selection();
	Gtk::TreeModel::iterator iter = refTreeSelection->get_selected();
	if(iter){
		Gtk::TreeModel::Row row = *iter;
		ostringstream ss;
		ss << row[m_columns.m_col_value];
		display_process(ss.str());
	}
	else{
		display_process("");
	}
}

void BatchDialog::on_process_chooser_cancel_clicked(){
	process_chooser->hide();
}

void BatchDialog::on_ok_button_clicked(){
	this->hide();
}

void BatchDialog::on_cancel_button_clicked(){
	this->hide();
}

void BatchDialog::display_process(string p){
	Glib::RefPtr<Gtk::TextBuffer> buffer = summary_textview->get_buffer();
	buffer->set_text("\n\n\n\n\n\n\n" + p);
	summary_textview->set_buffer(buffer);
	Gtk::TextBuffer::iterator iter = buffer->get_iter_at_line(buffer->get_line_count());
	iter++;
	summary_textview->scroll_to_iter(iter, 0.0);

	set_notebook(p);
	this->show_all_children(true);
}

void BatchDialog::set_notebook(string p){
	while(options_notebook->get_n_pages() > 1){
		options_notebook->remove_page(1);
	}
	options_notebook->show_all_children(true);
	switch(	Vars::processes[p]){
	case Processes::p_markerefficiency:
	{
		std::map< int, vector<int> > lists = Vars::process_options_map[Processes::p_markerefficiency];
		std::map<int, vector< int > >::iterator oiter;
		for(oiter = lists.begin(); oiter != lists.end(); oiter++){
			int test = oiter->first;
			switch(test){
			case OptionsGroup::og_thresholds:
			{
				ThresholdsModel* tm = new ThresholdsModel();
				vector<int> opts = (vector<int>) oiter->second;
				tm->set_fields(opts);
				options_notebook->append_page(*(tm->get_layout()), "Thresholds");
				options_notebook->show_all_children(true);
			}
			default:
				break;
			}
		}
		break;
	}
	case Processes::p_sampleefficiency:
	{
		break;
	}
	case -1:
		break;
	default:
		break;
	}
}

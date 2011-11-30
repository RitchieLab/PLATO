/*
 * BatchDialog.h
 *
 *  Created on: Jan 22, 2009
 *      Author: gilesjt
 */

#ifndef BATCHDIALOG_H_
#define BATCHDIALOG_H_

#include <stdlib.h>
#include <map>
#include <vector>
#include <algorithm>
#include <gtkmm.h>
#include <libglademm.h>
#include <string>
#include <string.h>
#include <sstream>
#include "Process.h"
#include "ThresholdsModel.h"
#include "Vars.h"
#include "SingleModelColumns.h"

using namespace std;

class BatchDialog : public Gtk::Dialog{
public:
	BatchDialog();
	BatchDialog::BatchDialog(BaseObjectType* cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade);
	virtual ~BatchDialog();

	vector<Process> get_batch(){return batch_processes;}
	void set_batch(vector<Process> b){batch_processes = b;}
	void initialize();

protected:
	Gtk::TextView* summary_textview;
	Gtk::TreeView* batch_treeview;
	Gtk::Notebook* options_notebook;
	Gtk::Button* plus_button;
	Gtk::Button* minus_button;
	Gtk::Button* ok_button;
	Gtk::Button* cancel_button;
	Gtk::Dialog* process_chooser;
	Gtk::Button* process_chooser_ok;
	Gtk::Button* process_chooser_cancel;
	Gtk::ComboBoxText* process_chooser_combobox;
	Gtk::VBox* process_chooser_vbox;

	SingleModelColumns m_columns;
	SingleModelColumns pc_columns;
	Glib::RefPtr<Gtk::ListStore> m_refTreeModel;
	Glib::RefPtr<Gtk::ListStore> pc_store;

	vector<Process> batch_processes;

	Gdk::Color color;

	virtual void on_plus_button_clicked();
	virtual void on_minus_button_clicked();
	virtual void on_process_chooser_ok_clicked();
	virtual void on_process_chooser_cancel_clicked();
	virtual void on_ok_button_clicked();
	virtual void on_cancel_button_clicked();
	virtual void on_batch_treeview_changed();

	void display_process(string);
	void set_notebook(string);
	Gtk::TreeModel::Row addProcess(string);
};

#endif /* BATCHDIALOG_H_ */

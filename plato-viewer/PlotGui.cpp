#include <Helper.h>
#include <map>
#include <iostream>
#include <sstream>
#include <glibmm.h>
#include <glibmm-2.4/glibmm.h>
#include <glibmm-2.4/glibmm/ustring.h>
#include <glibmm/ustring.h>
#include <gtkmm.h>
#include <gtkmm/window.h>
#include <gtkmm/dialog.h>
#include <gtkmm/scale.h>
#include <gtkmm/treeview.h>
//#include <gtkmm/drawingarea.h>
//#include <cairomm-1.0/cairomm/context.h>
#include <libglademm/xml.h>
#include <sigc++/class_slot.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "PlotGui.h"

using namespace std;

PlotGui::PlotGui(string f){

	try{
		refXml = Gnome::Glade::Xml::create(f);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}
//	myproject = new Project();

	plot_assistant = new PlotAssistant(f);

	plot_window = 0;
	plot_new_button = 0;
	refXml->get_widget("plot_window", plot_window);
	refXml->get_widget("plot_new_button", plot_new_button);

	if(plot_new_button){
		plot_new_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_plot_new_button_clicked));
	}

/*	if(save_testing_training_menu_item){
		save_testing_training_menu_item->signal_activate().connect(sigc::mem_fun(*this, &PlotGui::on_save_testing_training_menu_item_activate));
	}
	if(about_menu_item){
		about_menu_item->signal_activate().connect(sigc::mem_fun(*this, &PlotGui::on_about_menu_item_activate));
	}
	if(quit_menu_item){
		quit_menu_item->signal_activate().connect(sigc::mem_fun(*this, &PlotGui::on_quit_menu_item_activate));
	}
	if(load_all_button){
		load_all_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_load_all_button_clicked));
	}
	if(proportion_scale){
		proportion_scale->signal_value_changed().connect(sigc::mem_fun(*this, &PlotGui::on_proportion_scale_value_changed));
		proportion_scale->set_value(50.0);
	}
	if(load_test_files){
		load_test_files->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_load_test_files_clicked));
	}
	if(randomize_button){
		randomize_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_randomize_button_clicked));
	}
	if(view_sets_button){
		view_sets_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_view_sets_button_clicked));
	}
	if(set_stats_button){
		set_stats_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_set_stats_button_clicked));
	}
	if(right_variables_button){
		right_variables_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_right_variables_button_clicked));
	}
	if(left_variables_button){
		left_variables_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_left_variables_button_clicked));
	}
	if(logreg_train_button){
		logreg_train_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_logreg_train_button_clicked));
	}
	if(logreg_results_button){
		logreg_results_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_logreg_results_button_clicked));
	}
	if(mdr_run_button){
		mdr_run_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_mdr_run_button_clicked));
	}
	if(mdr_results_button){
		mdr_results_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_mdr_results_button_clicked));
	}
	if(mdr_save_button){
		mdr_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_mdr_save_button_clicked));
	}
	if(logreg_save_button){
		logreg_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_logreg_save_button_clicked));
	}
	if(folder_chooser){
		folder_chooser->add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
		folder_chooser->add_button("Select", Gtk::RESPONSE_OK);
	}
	if(summary_button){
		summary_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_summary_button_clicked));
	}
	if(summary_save_button){
		summary_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_summary_save_button_clicked));
	}
	if(bayclass_run_button){
		bayclass_run_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_bayclass_run_button_clicked));
	}
	if(bayclass_results_button){
		bayclass_results_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_bayclass_results_button_clicked));
	}
	if(bayclass_save_button){
		bayclass_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_bayclass_save_button_clicked));
	}
	if(set_stats_save_button){
		set_stats_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlotGui::on_set_stats_save_button_clicked));
	}
*/
	//Glib::thread_init();
}

PlotGui::~PlotGui(){
}

void PlotGui::on_plot_new_button_clicked(){
	plot_assistant->show();
}



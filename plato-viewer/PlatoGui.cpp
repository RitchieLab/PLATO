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
#include "PlatoGui.h"

using namespace std;

PlatoGui::PlatoGui(string f){

	try{
		refXml = Gnome::Glade::Xml::create(f);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}
//	myproject = new Project();
	project_config = new ProjectConfig(f);
	plot_gui = new PlotGui(f);

	main_window = 0;
	plot_tool_button = 0;
	new_tool_button = 0;
	edit_tool_button = 0;
	edit_batch_tool_button = 0;
	project_window = 0;
	project_tree_menu = 0;
	display_textview = 0;
	batch_dialog = 0;
	run_tool_button = 0;
	refXml->get_widget("main_app", main_window);
	refXml->get_widget("plot_tool_button", plot_tool_button);
	refXml->get_widget("new_tool_button", new_tool_button);
	refXml->get_widget("edit_tool_button", edit_tool_button);
	refXml->get_widget("edit_batch_tool_button", edit_batch_tool_button);
	refXml->get_widget("project_window", project_window);
	refXml->get_widget("project_tree_menu", project_tree_menu);
	refXml->get_widget("display_textview", display_textview);
	refXml->get_widget_derived("batch_dialog", batch_dialog);
	refXml->get_widget("run_tool_button", run_tool_button);

	eventbox1 = 0;
	textview1 = 0;
//	refXml->get_widget("eventbox1", eventbox1);
//	refXml->get_widget("textview1", textview1);

//	Gdk::Color color;
//	color.set_rgb_p(0.863, 0.855, 0.835);
//	textview1->modify_base(Gtk::STATE_NORMAL, color);




//	Gdk::Color color;
//	color.set_blue(50);
//	color.set_red(0);
//	color.set_green(50);
//	eventbox1->modify_base(Gtk::STATE_NORMAL, color);
//	eventbox1->modify_bg(Gtk::STATE_NORMAL, color);

	if(plot_tool_button){
		plot_tool_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_plot_tool_button_clicked));
		plot_tool_button->set_sensitive(false);
	}
	if(new_tool_button){
		new_tool_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_new_tool_button_clicked));
	}
	if(edit_tool_button){
		edit_tool_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_edit_tool_button_clicked));
		edit_tool_button->set_sensitive(false);
	}
	if(edit_batch_tool_button){
		edit_batch_tool_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_edit_batch_tool_button_clicked));
		edit_batch_tool_button->set_sensitive(true);
	}
	if(project_window){
		project_tree = new ProjectTree();
		project_window->add(*project_tree);
		project_window->show_all();
	}
	if(project_tree_menu){
		cout << "Have project-Tree_menu\n";
	}

/*	if(save_testing_training_menu_item){
		save_testing_training_menu_item->signal_activate().connect(sigc::mem_fun(*this, &PlatoGui::on_save_testing_training_menu_item_activate));
	}
	if(about_menu_item){
		about_menu_item->signal_activate().connect(sigc::mem_fun(*this, &PlatoGui::on_about_menu_item_activate));
	}
	if(quit_menu_item){
		quit_menu_item->signal_activate().connect(sigc::mem_fun(*this, &PlatoGui::on_quit_menu_item_activate));
	}
	if(load_all_button){
		load_all_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_load_all_button_clicked));
	}
	if(proportion_scale){
		proportion_scale->signal_value_changed().connect(sigc::mem_fun(*this, &PlatoGui::on_proportion_scale_value_changed));
		proportion_scale->set_value(50.0);
	}
	if(load_test_files){
		load_test_files->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_load_test_files_clicked));
	}
	if(randomize_button){
		randomize_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_randomize_button_clicked));
	}
	if(view_sets_button){
		view_sets_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_view_sets_button_clicked));
	}
	if(set_stats_button){
		set_stats_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_set_stats_button_clicked));
	}
	if(right_variables_button){
		right_variables_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_right_variables_button_clicked));
	}
	if(left_variables_button){
		left_variables_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_left_variables_button_clicked));
	}
	if(logreg_train_button){
		logreg_train_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_logreg_train_button_clicked));
	}
	if(logreg_results_button){
		logreg_results_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_logreg_results_button_clicked));
	}
	if(mdr_run_button){
		mdr_run_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_mdr_run_button_clicked));
	}
	if(mdr_results_button){
		mdr_results_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_mdr_results_button_clicked));
	}
	if(mdr_save_button){
		mdr_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_mdr_save_button_clicked));
	}
	if(logreg_save_button){
		logreg_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_logreg_save_button_clicked));
	}
	if(folder_chooser){
		folder_chooser->add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
		folder_chooser->add_button("Select", Gtk::RESPONSE_OK);
	}
	if(summary_button){
		summary_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_summary_button_clicked));
	}
	if(summary_save_button){
		summary_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_summary_save_button_clicked));
	}
	if(bayclass_run_button){
		bayclass_run_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_bayclass_run_button_clicked));
	}
	if(bayclass_results_button){
		bayclass_results_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_bayclass_results_button_clicked));
	}
	if(bayclass_save_button){
		bayclass_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_bayclass_save_button_clicked));
	}
	if(set_stats_save_button){
		set_stats_save_button->signal_clicked().connect(sigc::mem_fun(*this, &PlatoGui::on_set_stats_save_button_clicked));
	}
*/
	//Glib::thread_init();
}

PlatoGui::~PlatoGui(){
}

void PlatoGui::on_plot_tool_button_clicked(){
	plot_gui->show();
}

void PlatoGui::on_new_tool_button_clicked(){
//	if(plato_project){
//		throw MethodException("Are you sure you want to create a new project (all data will be lost)?");
//	}
	if(project_assistant){
		delete project_assistant;
	}
	project_assistant = new ProjectAssistant(this, Vars::PLATO_GUI_GLADE);

//	project_assistant->show();
}

void PlatoGui::on_project_assistant_apply(){
	if(project_assistant->apply()){
		plato_project = project_assistant->get_plato_project();
		plato_project.reinitialize();
		edit_tool_button->set_sensitive(true);
		if(plato_project.num_datasets() > 0){
			edit_batch_tool_button->set_sensitive(true);
		}
cout << "getting ready to display tree\n";
		display_project_tree();

cout << "Removing project_assistant\n";
		//project_assistant->hide();
//		if(project_assistant){
//			delete(project_assistant);
//		}
cout << "Loading data\n";
		Controller::load_data(&plato_project);
cout << "done loading assistant data\n";
		display_dataset_summary(plato_project.get_dataset_summary());
	}


}

void PlatoGui::on_edit_tool_button_clicked(){
	project_config->show();
}

void PlatoGui::on_edit_batch_tool_button_clicked(){
	batch_dialog->set_batch(plato_project.get_batch());
	batch_dialog->initialize();
	batch_dialog->show();
}


void PlatoGui::display_project_tree(){
//	Glib::RefPtr<Gtk::TreeStore> m_refTreeStore;
//	SingleModelColumns m_columns;
	project_tree->generate(&plato_project);
}

void PlatoGui::display_dataset_summary(string text){
	Glib::RefPtr<Gtk::TextBuffer> buffer = display_textview->get_buffer();
	buffer->set_text(buffer->property_text() + text);
	display_textview->set_buffer(buffer);
	Gtk::TextBuffer::iterator iter = buffer->get_iter_at_line(buffer->get_line_count());
	iter++;
	display_textview->scroll_to_iter(iter, 0.0);
}


#ifndef PLATOGUI_H
#define PLATOGUI_H

//#include <Helper.h>
#include <iostream>
#include <string>
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <gtkmm/window.h>
#include <glibmm.h>
#include <string.h>
//#include <DataSet.h>
//#include <StepOptions.h>
//#include <MethodException.h>
//#include <AlleleFrequency.h>
#include <vector>
#include <map>
#include "PlotGui.h"
#include "ProjectAssistant.h"
#include "ProjectConfig.h"
#include "PlatoProject.h"
#include "ProjectTree.h"
#include "BatchDialog.h"
#include "SingleModelColumns.h"
#include "Controller.h"


using namespace std;

class ProjectAssistant;

class PlatoGui : public Gtk::Window{

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;
	Gtk::Window* main_window;
	Gtk::ToolButton* plot_tool_button;
	Gtk::ToolButton* new_tool_button;
	Gtk::ToolButton* edit_tool_button;
	Gtk::ToolButton* edit_batch_tool_button;
	Gtk::ScrolledWindow* project_window;
	Gtk::Menu* project_tree_menu;
	Gtk::TextView* display_textview;
	BatchDialog* batch_dialog;
	Gtk::ToolButton* run_tool_button;

	Gtk::TextView* textview1;
	Gtk::EventBox* eventbox1;

	ProjectTree* project_tree;

	PlotGui* plot_gui;
	ProjectAssistant* project_assistant;
	ProjectConfig* project_config;

	PlatoProject plato_project;

	//Project* myproject;
	//
	virtual void on_plot_tool_button_clicked();
	virtual void on_new_tool_button_clicked();
	virtual void on_edit_tool_button_clicked();
	virtual void on_edit_batch_tool_button_clicked();
	void display_project_tree();
	void display_dataset_summary(string);

//	void convertToQuantiles(DataSet*, vector<vector<double> >*, vector<vector<double> >*);
//	int map_value(vector<double>*, double);

	//Glib::RefPtr<Gdk::Pixbuf> generateGraph(int, int);

public:
	PlatoGui(){};
	PlatoGui(string f);
	virtual ~PlatoGui();

	virtual void on_project_assistant_apply();

	Gtk::Window* getMainWindow(){return main_window;};

};

#endif


#ifndef PLOTGUI_H
#define PLOTGUI_H

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
#include "PlotAssistant.h"

using namespace std;

class PlotGui : public Gtk::Window{

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;
	Gtk::Window* plot_window;
	Gtk::ToolButton* plot_new_button;

	PlotAssistant* plot_assistant;

	virtual void on_plot_new_button_clicked();
	//Project* myproject;
	//

//	void convertToQuantiles(DataSet*, vector<vector<double> >*, vector<vector<double> >*);
//	int map_value(vector<double>*, double);

	//Glib::RefPtr<Gdk::Pixbuf> generateGraph(int, int);

public:
	PlotGui(string f);
	virtual ~PlotGui();


	void show(){plot_window->show();}
};

#endif


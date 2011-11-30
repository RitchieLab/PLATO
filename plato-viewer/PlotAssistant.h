#ifndef PLOTASSISTANT_H
#define PLOTASSISTANT_H

#include <iostream>
#include <string>
#include <list>
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string.h>

using namespace std;

class PlotAssistant{

	public:
		PlotAssistant(string);
		~PlotAssistant(){};

		void show(){plot_assistant->show();};

	protected:

		Glib::RefPtr<Gnome::Glade::Xml> refXml;
		Gtk::Assistant* plot_assistant;

		virtual void on_plot_assistant_apply();
		virtual void on_plot_assistant_cancel();

};


#endif

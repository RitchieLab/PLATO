#ifndef PROJECTCONFIG_H
#define PROJECTCONFIG_H

#include <iostream>
#include <string>
#include <list>
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string.h>

using namespace std;

class ProjectConfig{

	public:
		ProjectConfig(string);
		~ProjectConfig(){};

		void show(){project_config_dialog->show();};

	protected:

		Glib::RefPtr<Gnome::Glade::Xml> refXml;
		Gtk::Dialog* project_config_dialog;
//		Gtk::ComboBox* file_type;

};


#endif

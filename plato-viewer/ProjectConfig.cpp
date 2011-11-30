#include <gtkmm.h>
#include <gtkmm/liststore.h>
#include "ProjectConfig.h"

ProjectConfig::ProjectConfig(string f){
	try{
		refXml = Gnome::Glade::Xml::create(f);
	}
	catch(const Gnome::Glade::XmlError& ex){
		cout << "In catch\n";
	}

	project_config_dialog = 0;
//	file_type = 0;
	refXml->get_widget("project_config_dialog", project_config_dialog);
//	refXml->get_widget("file_type", file_type);

//	if(file_type){
}

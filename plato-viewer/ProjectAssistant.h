#ifndef PROJECTASSISTANT_H
#define PROJECTASSISTANT_H

#include <iostream>
#include <string>
#include <list>
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string.h>
//#include "OptionsList.h"
#include "PlatoProject.h"
#include "FileTypeModel.h"
#include "ProjectNameModel.h"
#include "LocusFilterModel.h"
#include "SampleFilterModel.h"
#include "OtherModel.h"
#include "CovTraitModel.h"
#include "PlatoGui.h"

using namespace std;

class PlatoGui;

class ProjectAssistant{

	public:
		ProjectAssistant(){};
		ProjectAssistant(PlatoGui*, string);
		~ProjectAssistant(){
//			if(project_assistant)
//				delete(project_assistant);
			if(ftm)
				delete(ftm);
			if(proj_model)
				delete(proj_model);
			if(locus_filter_model)
				delete(locus_filter_model);
			if(sample_filter_model)
				delete(sample_filter_model);
			if(other_model)
				delete(other_model);
			if(covtrait_model)
				delete(covtrait_model);
		};

		void show(){status = false; reset();project_assistant->show_all();};
		void hide(){project_assistant->hide();}
		PlatoProject get_plato_project(){return plato_project;}
		bool apply(){return status;}

	protected:

		Glib::RefPtr<Gnome::Glade::Xml> refXml;

		Gtk::Assistant* project_assistant;

//		Gtk::Entry* project_assistant_project_name;
//		Gtk::Layout* project_assistant_file_layout;
//		Gtk::VBox* project_assistant_locus_filter_vbox;
//		Gtk::VBox* project_assistant_sample_filter_vbox;
//		Gtk::VBox* project_assistant_covtrait_vbox;


		virtual void on_project_assistant_apply();
		virtual void on_project_assistant_cancel();
		void reset();
		PlatoProject plato_project;
		FileTypeModel* ftm;
		ProjectNameModel* proj_model;
		LocusFilterModel* locus_filter_model;
		SampleFilterModel* sample_filter_model;
		OtherModel* other_model;
		CovTraitModel* covtrait_model;

		PlatoGui* parent;

		bool status;

};


#endif

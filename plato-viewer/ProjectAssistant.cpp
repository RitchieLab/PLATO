#include <gtkmm.h>
#include "ProjectAssistant.h"

ProjectAssistant::ProjectAssistant(PlatoGui* par, string f){
	parent = par;
	try{
		refXml = Gnome::Glade::Xml::create(f);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	project_assistant = 0;
//	project_assistant_project_name = 0;
//	project_assistant_file_layout = 0;
//	project_assistant_locus_filter_vbox = 0;
	refXml->get_widget("project_assistant", project_assistant);
//	refXml->get_widget("project_assistant_project_name", project_assistant_project_name);
//	refXml->get_widget("project_assistant_file_layout", project_assistant_file_layout);
//	refXml->get_widget("project_assistant_locus_filter_vbox", project_assistant_locus_filter_vbox);

//	if(project_assistant_locus_filter_vbox){
/*		OptionsList* opts = new OptionsList(Vars::LOCUS_FILTER_GUI);
		option_lists.push_back(opts);
		Gtk::Table* table = new Gtk::Table(opts->num_opts(), 1, true);
		table->set_border_width(10);
		table->set_row_spacings(5);
		table->set_size_request(-1,-1);
		project_assistant_locus_filter_frame->add(*table);
		for(int o = 0; o < opts->num_opts(); o++){
			Gtk::Widget* mylayout = opts->get_option_data(o)->get_layout();
//			if(mylayout){
				table->attach(*mylayout, 0, 1, o, o + 1);
//				project_assistant_locus_filter_frame->add(*mylayout);
//			}
		}
		cout << opts->num_opts() << endl;
//		Gtk::Widget* mylayout = opts->get_option_data(1)->get_layout();
//		table->attach(*(opts->get_option_data(0)->get_layout()), 0, 1, 0, 1);
//		table->attach(*(opts->get_option_data(1)->get_layout()), 0, 1, 1, 2);
//		project_assistant_locus_filter_frame->add(*mylayout);
*/
//		project_assistant_locus_filter_vbox->show_all_children();
//	}

	ftm = new FileTypeModel();
	proj_model = new ProjectNameModel();
	locus_filter_model = new LocusFilterModel();
	sample_filter_model = new SampleFilterModel();
	other_model = new OtherModel();
	covtrait_model = new CovTraitModel();

	if(project_assistant){
        project_assistant->signal_apply().connect(sigc::mem_fun(*this, &ProjectAssistant::on_project_assistant_apply));
        project_assistant->signal_apply().connect(sigc::mem_fun(*parent, &PlatoGui::on_project_assistant_apply));
        project_assistant->signal_cancel().connect(sigc::mem_fun(*this, &ProjectAssistant::on_project_assistant_cancel));
		project_assistant->signal_close().connect(sigc::mem_fun(*this, &ProjectAssistant::on_project_assistant_cancel));

		project_assistant->append_page(*(proj_model->get_layout()));
		project_assistant->set_page_type(*(proj_model->get_layout()), Gtk::ASSISTANT_PAGE_CONTENT);
		project_assistant->set_page_title(*(proj_model->get_layout()), "Project Information");

		project_assistant->append_page(*(ftm->get_widget()));
		project_assistant->set_page_type(*(ftm->get_widget()), Gtk::ASSISTANT_PAGE_CONTENT);
		project_assistant->set_page_title(*(ftm->get_widget()), "File Input");

		project_assistant->append_page(*(locus_filter_model->get_layout()));
		project_assistant->set_page_type(*(locus_filter_model->get_layout()), Gtk::ASSISTANT_PAGE_CONTENT);
		project_assistant->set_page_title(*(locus_filter_model->get_layout()), "Locus Filters");

		project_assistant->append_page(*(sample_filter_model->get_layout()));
		project_assistant->set_page_type(*(sample_filter_model->get_layout()), Gtk::ASSISTANT_PAGE_CONTENT);
		project_assistant->set_page_title(*(sample_filter_model->get_layout()), "Sample Filters");

		project_assistant->append_page(*(covtrait_model->get_layout()));
		project_assistant->set_page_type(*(covtrait_model->get_layout()), Gtk::ASSISTANT_PAGE_CONTENT);
		project_assistant->set_page_title(*(covtrait_model->get_layout()), "Covariate/Trait Data");

		project_assistant->append_page(*(other_model->get_layout()));
		project_assistant->set_page_type(*(other_model->get_layout()), Gtk::ASSISTANT_PAGE_CONTENT);
		project_assistant->set_page_title(*(other_model->get_layout()), "Other Options");

		Gtk::Label* confirm = new Gtk::Label("Please click \"Apply\" to finalize the project creation.");
		project_assistant->append_page(*confirm);
		project_assistant->set_page_type(*confirm, Gtk::ASSISTANT_PAGE_CONFIRM);
//		Gtk::Widget* page = project_assistant->get_nth_page(0);
//		project_assistant->set_page_complete((*page), true);
//		page = project_assistant->get_nth_page(project_assistant->get_n_pages() - 1);
//		project_assistant->set_page_complete((*page), true);

		Gtk::Widget* page;
		//reset page flags
		for(int i = 0; i < project_assistant->get_n_pages(); i++){
			page = project_assistant->get_nth_page(i);
//			if(i == 0 || i == (project_assistant->get_n_pages() - 1)){
				project_assistant->set_page_complete((*page), true);
//			}
//			else{
//				project_assistant->set_page_complete((*page), false);
//			}
		}
	}

///	project_assistant_file_layout->add(*(ftm->get_widget()));
///	project_assistant_file_layout->show_all_children();


	status = false;
	project_assistant->show_all();
}

//void ProjectAssistant::on_project_assistant_project_name_changed(){
//	Gtk::Widget* page = project_assistant->get_nth_page(project_assistant->get_current_page());
//	if(project_assistant_project_name->get_text_length() > 0){
//		project_assistant->set_page_complete((*page), true);
//	}
//	else{
//		project_assistant->set_page_complete((*page), false);
//	}
//}

void ProjectAssistant::on_project_assistant_apply(){
	status = true;
	plato_project.reset();
	plato_project.set_project_name(proj_model->get_name());
cout << "done assigning name\n";
	ftm->assign_results(&plato_project);


	//plato_project.set_project_name(project_assistant_project_name->get_text());

}

void ProjectAssistant::on_project_assistant_cancel(){
	project_assistant->hide();
}

void ProjectAssistant::reset(){
	Gtk::Widget* page;
	//reset page flags
	for(int i = 0; i < project_assistant->get_n_pages(); i++){
		page = project_assistant->get_nth_page(i);
//		if(i == 0 || i == (project_assistant->get_n_pages() - 1)){
			project_assistant->set_page_complete((*page), true);
//		}
//		else{
//			project_assistant->set_page_complete((*page), false);
//		}
	}

//	project_assistant_project_name->set_text("");
//	delete(ftm);
//	ftm = new FileTypeModel();

//	project_assistant_file_layout->add(*(ftm->get_widget()));
//	project_assistant_file_layout->show_all_children();

	status = false;
}

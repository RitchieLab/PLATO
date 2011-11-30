#include <gtkmm.h>
#include "PlotAssistant.h"

PlotAssistant::PlotAssistant(string f){
	try{
		refXml = Gnome::Glade::Xml::create(f);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	plot_assistant = 0;
	refXml->get_widget("plot_assistant", plot_assistant);
	if(plot_assistant){
        plot_assistant->signal_apply().connect(sigc::mem_fun(*this, &PlotAssistant::on_plot_assistant_apply));
		plot_assistant->signal_cancel().connect(sigc::mem_fun(*this, &PlotAssistant::on_plot_assistant_cancel));
		plot_assistant->signal_close().connect(sigc::mem_fun(*this, &PlotAssistant::on_plot_assistant_cancel));
		Gtk::Widget* page = plot_assistant->get_nth_page(0);
		plot_assistant->set_page_complete((*page), true);
	}

}


void PlotAssistant::on_plot_assistant_apply(){
}

void PlotAssistant::on_plot_assistant_cancel(){
	plot_assistant->hide();
}

/*
 * plato-viewer.cpp
 *
 *  Created on: Dec 12, 2008
 *      Author: gilesjt
 */

#include <gtkmm.h>
#include <gtkmm/main.h>
#include <gdkmm-2.4/gdkmm.h>
#include "PlatoGui.h"
#include "Vars.h"


int main(int argc, char *argv[]){
	Vars::initialize();
	Glib::thread_init();
	gdk_threads_init();
	gdk_threads_enter();
	g_assert(Vars::mutex_Project == NULL);
	Vars::mutex_Project = g_mutex_new();

	Gtk::Main kit(argc, argv);
	PlatoGui* gui = new PlatoGui(Vars::PLATO_GUI_GLADE);
	Gtk::Window* window = gui->getMainWindow();
	kit.run(*window);
	gdk_threads_leave();
	return 0;
}

/*
 * ProjectTree.h
 *
 *  Created on: Jan 19, 2009
 *      Author: gilesjt
 */

#ifndef PROJECTTREE_H_
#define PROJECTTREE_H_
#include <gtkmm.h>
#include "PlatoProject.h"
#include "SingleModelColumns.h"

using namespace std;

class ProjectTree : public Gtk::TreeView{
public:
	ProjectTree();
	virtual ~ProjectTree();
	void generate(PlatoProject*);

protected:
	virtual bool on_button_press_event(GdkEventButton *ev);
	virtual void on_menu_file_popup_generic();

	SingleModelColumns m_columns;
	Glib::RefPtr<Gtk::TreeStore> m_refTreeModel;
	Gtk::Menu m_Menu_Popup;
};

#endif /* PROJECTTREE_H_ */

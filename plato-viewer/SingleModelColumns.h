/*
 * SingleModelColumns.h
 *
 *  Created on: Sep 11, 2008
 *      Author: gilesjt
 */

#ifndef SINGLEMODELCOLUMNS_H_
#define SINGLEMODELCOLUMNS_H_
#include <gtkmm.h>
#include <vector>

using namespace std;

class SingleModelColumns : public Gtk::TreeModel::ColumnRecord{
public:
	SingleModelColumns(){add(m_col_value); add(m_col_id);}
	virtual ~SingleModelColumns();
	Gtk::TreeModelColumn<int> m_col_id;
	Gtk::TreeModelColumn<Glib::ustring> m_col_value;


};

#endif /* SINGLEMODELCOLUMNS_H_ */



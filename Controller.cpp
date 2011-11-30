/*
 * Controller.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: cozartc
 */
#include "Controller.h"
#include "libsqlitewrapped.h"

using namespace Methods;

Controller::Controller()
{
}

Controller::~Controller()
{
}

void Controller::drop_table(Database* db, string table)
{
	string sql = "SELECT NAME FROM SQLITE_MASTER WHERE TYPE = \"table\" AND NAME = \"" + table + "\"";
	Query query(*db);
	if(!query.execute(sql.c_str()))
	{
		//throw PlatoViewerException(query.GetError());
		throw MethodException(query.GetError());
	}
	string tableName = "";
	if(query.fetch_row())
	{
		tableName = query.getval(0);
	}
	query.free_result();
	if(tableName == table)
	{
		sql = "DROP TABLE " + table;
		execute_sql(query, sql);
	}
	else
	{
		cout << "Could not find table: " << table << " " << tableName << endl;
	}
}//end drop_table method

long long Controller::execute_sql(Query& q, string sql)
{
	if (!q.execute(sql))
	{
		cout<<q.GetLastQuery();
		cout<<q.GetError();
		string errorString = q.GetError();
		throw MethodException(sql + "    " + errorString.substr(errorString.length() - 201, 200));
	}
	return (long long)q.insert_id();
}

/*
 * DataLoader.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef METHODS_DATA_LOADER_H
#define METHODS_DATA_LOADER_H

#include <string>

namespace Methods {

class DataLoader {

public:
	DataLoader();

	void read();

private:
	void readPed(const std::string& fn);


};

}

#endif /* DATALOADER_H_ */

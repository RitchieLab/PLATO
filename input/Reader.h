#ifndef READER_H
#define READER_H

/*!
 * This class is an abstract base class for reading files into a dataset from
 * the filesystem
 */
class Reader{

public:

	Reader();

	virtual ~Reader() = 0;

	virtual DataSet& read(DataSet& ds) = 0;

};

#endif

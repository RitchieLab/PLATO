#include <iostream>
#include <string>
#include <vector>
#include "General.h"
namespace Methods{


void General::Tokenize(const string& str, vector<string>& tokens, const string& delimiter){
	string::size_type lastPos = str.find_first_not_of(delimiter, 0);
	string::size_type pos = str.find_first_of(delimiter, lastPos);

	while(string::npos != pos || string::npos != lastPos){
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiter, pos);
        pos = str.find_first_of(delimiter, lastPos);
    }
}

vector<string> General::ParseDelimitedLine(string line){
    vector<string> temp_tokens;
    vector<string> tokens;
    Tokenize(line, temp_tokens, " ");
    for(int i = 0; i < (int)temp_tokens.size(); i++){
        Tokenize(temp_tokens[i], tokens, "\t");
    }
    return tokens;
}

vector<string> General::ParseDelimitedLine(string line, string delim){
    vector<string> tokens;
    Tokenize(line, tokens, delim);
    return tokens;
}

string General::getStringInt(int v){
	stringstream os;
	os << v;
	string temp;
	os >> temp;
	return temp;
}
string General::getStringFloat(float v){
	stringstream os;
	os << v;
	string temp;
	os >> temp;
	return temp;
}

//template <class T> const std::string getString(const T& t){
//	std::stringstream os;
//	os << t;
//	return os.str();
//}
/*
string General::getString(const unsigned int v){
	stringstream os;
	os << v;
	string temp;
	os >> temp;
	return temp;
}
string General::getString(const int v){
	stringstream os;
	os << v;
	string temp;
	os >> temp;
	return temp;
}
string General::getString(const float v){
	stringstream os;
	os << v;
	string temp;
	os >> temp;
	return temp;
}
string General::getString(const double v){
	stringstream os;
	os << v;
	string temp;
	os >> temp;
	return temp;
}
*/
/*
string General::binary(int number) {
    int remainder;

    if(number <= 1) {
        stringstream os;
        os << number;
        string temp;
        os >> temp;
        return temp;
    }

    remainder = number%2;
    string val = binary(number >> 1);
    stringstream os;
    os << remainder;
    string temp;
    os >> temp;
    val += temp;
    return val;
}
*/


/*dynamic_bitset<> General::getBits(int v){
	string val = General::binary(v);
	dynamic_bitset<> b(val);
//	b.resize(val.size());
//	for(int i = 0; i < val.size(); i++){
//		if(val[i] == '1'){
//			b[i] = true;
//		}
//		else{
//			b[i] = false;
//		}
//	}

	return b;
}
*/

}

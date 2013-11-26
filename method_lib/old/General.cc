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
        Tokenize(temp_tokens.at(i), tokens, "\t");
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

}

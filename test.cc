#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "cdflib.h"
//#include "helper.h"
using namespace std;
#define PI 3.141592653589793238462643

float _subuprob(float);
float _subchisqrprob(int, float);

int main (int argc, char* argv[]){
	if(argc < 3){
		exit(1);
	}
/*	string filename = argv[1];
	
	ifstream input;
	input.open(filename.c_str(), ios::in);
	float elems[6 + 496963];
	int read = 0;
	while(!input.eof()){
		input >> elems[0];//famid
		input >> elems[1];//indid
		input >> elems[2];//dad
		input >> elems[3];//mom
		input >> elems[4];//gender
		input >> elems[5];//aff
		cout << elems[0] << endl << elems[1] << endl << elems[2] << endl << elems[3] << endl << elems[4] << endl << elems[5] << endl;
		for(int i = 6; i < 496969; i++){
			input >> elems[i];
			cout << elems[i] << endl;
		}
		read++;
		if(read == 2)
		break;
	}
	cout << "read = " << read;
exit(1);	
*/
	string deg = argv[1];
	string chi = argv[2];

	double x = atof(chi.c_str());
	double df = atof(deg.c_str());
	double p, q, bound;
	int code = 1, status;
//	float result = _subchisqrprob(atoi(deg.c_str()), atof(chi.c_str()));
	//long double result = chiprobP(atof(chi.c_str()), atoi(deg.c_str()));
	long double result = cdfchi(&code, &p, &q, &x, &df, &status, &bound);
	//ofstream out ("testout.out");
	//out.precision(4);
	stringstream s2;
	s2 << result;
	cout << "Result of " << deg << ", " << chi << "\t" << s2.str() << endl;
	//out.close();
}

float _subchisqrprob(int deg, float chi){
	float pval = 0.0;
	if(chi <= 0){
		pval = 1;
	}
	else if(deg > 100){
		pval = _subuprob((pow((chi / (float)deg),((float)1/(float)3)) - (1 - 2/9/deg)) / sqrt((float)2/(float)9/(float)deg));
	}
	else if(chi > 400){
		pval = 0;
	}
	else{
		float a;
		int i1;

		if((deg % 2) != 0){
			pval = 2 * _subuprob(sqrt(chi));
			a = sqrt(2.0/PI) * exp(-chi/2.0) / sqrt(chi);
			i1 = 1;
		}
		else{
			pval = a = exp(-chi/2.0);
			i1 = 2;
		}

		for(int i = i1; i <= (deg - 2); i += 2){
			a *= chi / (float)i;
			pval += a;
		}
	}

	return pval;
}



float _subuprob(float val){
	float p = 0;
	float absx = (float) fabs(val);

	if(absx < 1.9){
		p = pow((1 +
			absx * (.049867347
			+ absx * (.0211410061
			+ absx * (.0032776263
			+ absx * (.0000380036
			+ absx * (.0000488906
			+ absx * .000005383)))))), -16.0/2.0);
	}
	else if(absx <= 100){
		for(int i = 18; i >= 1; i--){
			p = i / (absx + p);
		}
		p = exp(-.5 * absx * absx) / sqrt(2 * PI) / (absx + p);
	}
	if(val < 0){
		p = 1 - p;
	}
	return p;
}


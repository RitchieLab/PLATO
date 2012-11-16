/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates a genotype efficiency for all markers.
*
*
* Files generated:
*	percent_breakdown_by_marker.txt
*	percent_breakdown_by_chrom.txt
*       post_marker_geno_eff_filter_summary.txt
*
*File: LD.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "LD.h"
#include "Options.h"
#include "General.h"
#include "Helper.h"
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"

namespace Methods{
void LD::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void LD::PrintSummary(){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		(*markers)[m]->setFlag(false);
	}

}

void LD::filter(){
}

//courtesy of plink
void LD::calcSetMeanVariance(vector<Marker*> pSNP, vector<double> & mean, vector<vector<double> > & variance){
	int flagged_n = 0;

	int nss = pSNP.size();

	if (nss<1)
	   return;

	mean.resize(nss,0);

	vector<int> cnt(nss,0);
	variance.resize(nss);

	for (int j=0; j<nss; j++)
		variance[j].resize(nss,0);

	////////////////////////////
	//Iterate over SNPs in SET
	//
	int msize = pSNP.size();
	int fsize = families->size();
////	int j = 0;
	for(int m = 0; m < msize; m++){
		Marker* mark = pSNP[m];
		if(mark->isEnabled() && !mark->isFlagged()){
			int mloc = mark->getLoc();
			for(int f = 0; f < fsize; f++){
				Family* fam = (*families)[f];
				if(fam->isEnabled()){
					vector<Sample*>* founders = fam->getFounders();
					for(int s = 0; s < (int)founders->size(); s++){
						Sample* samp = (*founders)[s];
						if(samp->isEnabled()){
							//increase general, flagged sample size (only count once)
							if(m == 0){
								flagged_n++;
							}
							if(!mark->isMicroSat()){
								if(samp->getAone(mloc)){
									if(samp->getAtwo(mloc) && !samp->getAmissing(mloc)){ //11 homozyg
										mean[m]++;
										cnt[m]++;
									}
								}
								else{
									cnt[m]++;
									if(!samp->getAtwo(mloc)){ //00 homozyg
										mean[m]--;
									}
								}
							}
							else{
								if(samp->getAbone(mloc) == samp->getAbtwo(mloc) && samp->getAbone(mloc) != -1){
									mean[m]++;
									cnt[m]++;
								}
								else{
									cnt[m]++;
								}
							}
						}
					}
				}
			}
		}
	}

	  // Having iterated over all individuals, we can now calculate the mean
	  //   // values, perform mean-substitution of missing data, and calculate the
	  //     // second order terms
	  //

	for(int j = 0; j < nss; j++){
		mean[j] /= (double)cnt[j];
	}


	  /////////////////////////////////////
	  //  // Iterate over pairs of SNPs in SET
	  //
	  //    // First SNP
	  //
	for(int m1 = 0; m1 < (int)pSNP.size(); m1++){
		Marker* mark1 = pSNP[m1];
		if(mark1->isEnabled() && !mark1->isFlagged()){
			int mloc1 = mark1->getLoc();
			for(int m2 = m1; m2 < (int)pSNP.size(); m2++){
				Marker* mark2 = pSNP[m2];
				if(mark2->isEnabled() && !mark2->isFlagged()){
					int mloc2 = mark2->getLoc();

					for(int f = 0; f < fsize; f++){
						Family* fam = (*families)[f];
						if(fam->isEnabled()){
							vector<Sample*>* founders = fam->getFounders();
							for(int s = 0; s < (int)founders->size(); s++){
								Sample* samp = (*founders)[s];
								if(samp->isEnabled()){
									double v1 = mean[m1], v2=mean[m2];
									//first snp
									if(!mark1->isMicroSat()){
										if(samp->getAone(mloc1)){
											if(samp->getAtwo(mloc1) && !samp->getAmissing(mloc1)){ //11 homozyg
												v1 = 1;
											}
										}
										else{
											if(!samp->getAtwo(mloc1)){ //00 homozyg
												v1 = -1;
											}
											else{
												v1 = 0; //01 het
											}
										}
									}
									else{
										if(samp->getAbone(mloc1) == samp->getAbtwo(mloc1) && samp->getAbone(mloc1) != -1){
											v1 = 1;
										}
										else if(samp->getAbone(mloc1) != -1){
											v1 = 0;
										}
									}


									//second snp
									if(!mark2->isMicroSat()){
										if(samp->getAone(mloc2)){
											if(samp->getAtwo(mloc2) && !samp->getAmissing(mloc2)){ //11 homozyg
												v2 = 1;
											}
										}
										else{
											if(!samp->getAtwo(mloc2)){ //00 homozyg
												v2 = -1;
											}
											else{
												v2 = 0; //01 het
											}
										}
									}
									else{
										if(samp->getAbone(mloc2) == samp->getAbtwo(mloc2) && samp->getAbone(mloc2) != -1){
											v2 = 1;
										}
										else if(samp->getAbone(mloc2) != -1){
											v2 = 0;
										}
									}
									//contribution to covariance term
									//mean substitution
									variance[m1][m2] += (v1 - mean[m1]) * (v2 - mean[m2]);
								}
							}
						}
					}
				}
			}
		}
	}

	//make symmetric covariance matrix
	for(int i = 0; i < nss; i++){
		for(int j = i; j < nss; j++){
			variance[i][j] /= (double)(flagged_n);
			variance[j][i] = variance[i][j];
		}
	}

	return;
}

/*vector< vector<double> > LD::svd_inverse(vector< vector<double> > & u){
  const double eps = 1e-12;

  if (u.size() == 0){
    cerr << "Internal problem: matrix with no rows (inverse function)\n";
  	exit(1);
  }
  if (u.size() != u[0].size() ){
    cerr << "Internal problem: Cannot invert non-square matrix\n";
	exit(1);
  }
  int n = u.size();

  vector<double> w(n,0);

  vector<vector<double> > v(n);
  for (int i=0; i<n; i++)
    v[i].resize(n,0);

  svdcmp(u,w,v);//////!!!!!!!!!!!!!!!!!!!!!!!

  // Look for singular values
  double wmax = 0;
  for (int i=0; i<n; i++)
    wmax = w[i] > wmax ? w[i] : wmax;
  double wmin = wmax * eps;
  for (int i=0; i<n; i++)
  {
//       cout << w[i] << "\n";
//       //       if ( w[i] < wmin ) cout << "FLAGGIN!\n";
//
    w[i] = w[i] < wmin ? 0 : 1/w[i];
  }

    // u w t(v)
	//
	//   // row U * 1/w
	//

  vector<vector<double> > r(n);
  for (int i=0; i<n; i++)
  {
	r[i].resize(n,0);
    for (int j=0; j<n; j++)
      u[i][j] = u[i][j] * w[j];
  }

 // [nxn].[t(v)]

  for (int i=0; i<n; i++)
	  for (int j=0; j<n; j++)
		  for (int k=0; k<n; k++)
			  r[i][j] += u[i][k] * v[j][k];

  return r;
}

void LD::svdcmp(vector<vector<double> > & a,
		        vector<double> & w,
				        vector<vector<double> > &v){
  bool flag;
  int i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double volatile temp;

  int m=a.size();
  if (m==0){ cerr << "Internal problem in SVD function (no observations left?)\n"; exit(1);}
  int n=a[0].size();

  vector<double> rv1(n);
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k][i]);
      if (scale != 0.0) {
        for (k=i;k<m;k++) {
	      a[k][i] /= scale;
	      s += a[k][i]*a[k][i];
    	}
		f=a[i][i];
		g = -SIGN(sqrt(s),f);
		h=f*g-s;
		a[i][i]=f-g;
		for (j=l-1;j<n;j++) {
		  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
		  f=s/h;
		  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
		}
		for (k=i;k<m;k++) a[k][i] *= scale;
	  }
	}
	w[i]=scale *g;
	g=s=scale=0.0;
	if (i+1 <= m && i+1 != n) {
	  for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
	  if (scale != 0.0) {
	    for (k=l-1;k<n;k++) {
	      a[i][k] /= scale;
	      s += a[i][k]*a[i][k];
	    }
	    f=a[i][l-1];
	    g = -SIGN(sqrt(s),f);////////////!!!!!!!!!!!!!
		h=f*g-s;
		a[i][l-1]=f-g;
		for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
    	for (j=l-1;j<m;j++) {
		  for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
		  for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
		}
	    for (k=l-1;k<n;k++) a[i][k] *= scale;
	  }
	}
	anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));///////!!!!!!!!!!!!!
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
        for (j=l;j<n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<n;j++) {
          for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	      for (k=l;k<n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
        for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
        f=(s/a[i][i])*g;
	    for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	  }
	  for (j=i;j<m;j++) a[j][i] *= g;
	} else for (j=i;j<m;j++) a[j][i]=0.0;
	++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=true;
      for (l=k;l>=0;l--) {
        nm=l-1;
	    temp=fabs(rv1[l])+anorm;
	    if (temp == anorm) {
	      flag=false;
	      break;
	    }
	    temp=fabs(w[nm])+anorm;
	    if (temp == anorm) break;
	  }
	  if (flag) {
    	c=0.0;
	    s=1.0;
		for (i=l;i<k+1;i++) {
		  f=s*rv1[i];
		  rv1[i]=c*rv1[i];
		  temp = fabs(f)+anorm;
		  if (temp == anorm) break;
		  g=w[i];
		  h=pythag(f,g);
		  w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s = -f*h;
		  for (j=0;j<m;j++) {
		    y=a[j][nm];
		    z=a[j][i];
		    a[j][nm]=y*c+z*s;
		    a[j][i]=z*c-y*s;
		  }
		}
	  }
	  z=w[k];
	  if (l == k) {
	    if (z < 0.0) {
	      w[k] = -z;
	      for (j=0;j<n;j++) v[j][k] = -v[j][k];
	    }
        break;
	  }
	  if (its == 29){ cerr << "SVD function cannot converge: multicollinearity issues?\n"; exit(1);}
	  x=w[l];
	  nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
	    y=w[i];
	    h=s*g;
	    g=c*g;
	    z=pythag(f,h);/////////!!!!!!!!!!!!!!!!!
        rv1[j]=z;
	    c=f/z;
	    s=h/z;
	    f=x*c+g*s;
	    g=g*c-x*s;
	    h=y*s;
	    y *= c;
	    for (jj=0;jj<n;jj++) {
	      x=v[jj][j];
	      z=v[jj][i];
	      v[jj][j]=x*c+z*s;
	      v[jj][i]=z*c-x*s;
	    }
	    z=pythag(f,h);
	    w[j]=z;
	    if (z) {
	      z=1.0/z;
	      c=f*z;
	      s=h*z;
	    }
	    f=c*g+s*y;
	    x=c*y-s*g;
	    for (jj=0;jj<m;jj++) {
	      y=a[jj][j];
	      z=a[jj][i];
	      a[jj][j]=y*c+z*s;
	      a[jj][i]=z*c-y*s;
	    }
	  }
	  rv1[l]=0.0;
	  rv1[k]=f;
	  w[k]=x;
	}
  }
}

double LD::pythag(const double a, const double b){
  double absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

double LD::SQR(double a){
	return a*a;
}
*/

vector<bool> LD::vif_prune(vector<vector<double> > m, double threshold){
	//number of variables
	int p = m.size();
	//cout << "cur size = " << p << endl;
	vector<bool> cur(p, true);

	//this on ly is needed if we have 2+ snps
	if(p < 2){
		return cur;
	}

	vector<vector<double> > r = m;

	// Make 'm' a correlation matrix
	for (int i=0; i<p; i++)
		for (int j=0; j<p; j++)
	        r[i][j] = m[i][j] / sqrt(m[i][i] * m[j][j]);

	// Number of excluded items
	int it = 0;

	// Any SNPs with zero variance should be automatically excluded
	for (int i=0; i<p; i++)
	    if ( r[i][i] == 0 || !realnum(r[i][i]) )
	    {
			//cout << r[i][i] << " : setting " << i << " false\n";
	      cur[i] = false;
	      it++;
	    }


	// For any pair of perfectly correlated SNPs, exclude 1
	while(1)
	{
		bool done = true;
		for (int i=0; i<p-1; i++)
		{
			if (cur[i]){
			    for (int j=i+1;j<p; j++)
		        {
			        if (cur[j]) {
				        if ( fabs(r[i][j]) > options.getLDPWthreshold() )
				        {
							//cout << fabs(r[i][j]) << " > " << options.getLDPWthreshold() << " at " << i << endl;
					        cur[i] = false;
				            it++;
				            done = false;
				            break;
				        }
				    }
				}
	        }
    	}
     	if (done) break;
    }

	//skip vif caculation?
	if(options.doLDpairwise()){
		return cur;
	}

  // Calculate r^2 for each element versus all others
  // considering only the current non-pruned elements
  //
    while (1)
    {
      // Build correlation matrix all included items
		vector<vector<double> > u;

		for (int i=0;i<p;i++)
		{
			if ( cur[i] )
			{
				vector<double> mt;
		        for (int j=0;j<p;j++)
		 	       if ( cur[j] )
				       mt.push_back(r[i][j]);
				u.push_back(mt);
			}
		}

//cout << "p = " << p << " u = " << u.size() << endl;
		// Check enough markers left
		if (u.size()<2) break;
//      cout <<  " about to invert\n";
//      cout.precision(12);
//      for (int i=0; i<u.size(); i++)
//      {
//	      for (int j=0; j<u[i].size(); j++)
//		      cout << u[i][j] << " ";
//		      cout << "\n";
//	  }


		// Get inverse
		u = svd_inverse(u);//////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		//Calculate VIFs
        double maxVIF = 0;
        int maxI = 0;
        int cnt=0;
        for (int i=0;i<p;i++)
		    if ( cur[i] )
		    {
        // r^2 = 1 - 1/x where x is diagonal element of inverted
		//         // correlation matrix
		//                 // As VIF = 1 / ( 1 - r^2 ) , implies VIF = x
		//
		        double vif = u[cnt][cnt];
		        if ( dLess(maxVIF, vif) )
                {
		            maxVIF = vif;
		            maxI = i;
		        }
		        cnt++;
            }

		// How are we doing?
		if ( dGreater(maxVIF,threshold) )
		{
			// exclude this item
			cur[maxI] = false;
		}
		else{
			break;
		}

		// Increase count of removed items
		it++;
		// Down to a single item or worse?
		if (it==p-1) break;
	}
	return cur;

}

vector<int> LD::getChromosomeMarkerRange(vector<Marker*>* marks, int chrom){
	vector<int> m(2);
	m[0] = -1;
	m[1] = -1;
	for(int i = 0; i < (int)marks->size(); i++){
		Marker* mark = (*marks)[i];
		if(mark->isEnabled() && mark->getChrom() == chrom){
			if(i <= m[0] || m[0] == -1) m[0] = i;
			if(i >= m[1] || m[1] == -1) m[1] = i;
		}
	}
	return m;
}

void LD::setMarkerRange(){
	//if(run_chr > 0){
		vector<int> m = getChromosomeMarkerRange(markers, run_chr);
		if(m[0] == -1 || m[1] == -1){
			cerr << "blah blah blha error getchrommarkrange\n";
			exit(1);
		}
		run_start = m[0];
		run_end = m[1];
	//}
}

void LD::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	if(options.doLDCalcOnly()){
		calcLDStatistics();
		return;
	}

	map<int, int> chroms;
	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < (int)markers->size(); m++){
		if((*markers)[m]->getChrom() < opts::_CHRX_ && isValidMarker((*markers)[m], &options, prev_base, prev_chrom)){
			if((*markers)[m]->isEnabled() && !(*markers)[m]->isFlagged()){
				chroms[(*markers)[m]->getChrom()] = 1;
			}
		}
	}

	if(options.doLDpairwise() || options.doLDvif()){
		string filename = opts::_OUTPREFIX_ + "ld_keep" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			filename += "." + getString<int>(order);
		}
		string filename2 = opts::_OUTPREFIX_ + "ld_removed_" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			filename2 += "." + getString<int>(order);
		}
		ofstream LDIN(filename.c_str(), ios::out);
		ofstream LDOUT(filename2.c_str(), ios::out);
		if(!LDIN){
			opts::printLog("Error opening " + filename + ". Exiting!\n");
			exit(1);
		}
		if(!LDOUT){
			opts::printLog("Error opening " + filename2 + ".  Exiting!\n");
			exit(1);
		}
////		int win_start = 0;
////		int win_end = win_start + options.getLDWin();

		vector<bool> include(markers->size(), true);
		map<int, int>::iterator chrom;
		opts::printLog("Total number of markers enabled at beginning of LD calculations: " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
		opts::printLog("Looking at ALL markers per chromosome, skipping those that are currently removed or disabled...\n");
		for(chrom = chroms.begin(); chrom != chroms.end(); chrom++){
			run_chr = chrom->first;

			setMarkerRange();
			int s1 = run_start;
			int s2 = run_start + options.getLDWin() - 1;
//cout << s1 << " " << (*markers)[s1]->getRSID() << " " << (*markers)[s1]->getBPLOC() << " - " << s2 << " " << (*markers)[s2]->getRSID() << " " << (*markers)[s2]->getBPLOC() << endl;
			while(s2 <= run_end){
				vector<Marker*> pSNP(0);
				for(int l = s1; l <= s2; l++){
					if(include[l] && (*markers)[l]->isEnabled() && !(*markers)[l]->isFlagged()){
						pSNP.push_back((*markers)[l]);
					}
		//			if(!(*markers)[l]->isEnabled()){
		//				include[l] = false;
		//			}
				}

				if(pSNP.size() < 2){
					if(s2 == run_end){
						break;
					}
					s1 += options.getLDstep();
					s2 += options.getLDstep();

					if(s2 > run_end){
						s2 = run_end;
					}

					if(s2 - s1 < 1){
						break;
					}

					continue;
				}

				vector<double> mean;
				vector<vector<double> > variance;

				cout << "Pruning SNPs " << s1-run_start+1 << " to " << s2-run_start+1 << " of " << run_end - run_start+1 << "       \r";
				cout.flush();

				calcSetMeanVariance(pSNP,mean,variance);

				vector<bool> cur = vif_prune(variance, options.getLDthreshold());
				//exit(1);

				int k=0;
				for(int l=s1; l<=s2; l++){
					//update main list bu do not get back already excluded snps
					if(include[l] && !cur[k++]){
						include[l] = false;
					}
				}

				//advance window
				if(s2 == run_end){
					break;
				}

				s1 += options.getLDstep();
				s2 += options.getLDstep();

				if(s2 > run_end){
					s2 = run_end;
				}
				if(s2 - s1 < 1){
					break;
				}
			}// next window

			cout << "\n";

			// record what is in, what is out
			int cnt_in = 0, cnt_out = 0, cnt_dis = 0, cnt_din = 0, cnt_dout = 0;
			for(int l =run_start; l<= run_end; l++){
				if(!(*markers)[l]->isEnabled()){
					cnt_dis++;
				}
			}
			for(int l = run_start; l<= run_end; l++){
				if(include[l]){
					//output to IN file
					LDIN << (*markers)[l]->toString() << endl;
					cnt_in++;
					if((*markers)[l]->isEnabled() && !(*markers)[l]->isFlagged()){
						cnt_din++;
					}
				}
				else{
					//output to OUT file
					LDOUT << (*markers)[l]->toString() << endl;
					cnt_out++;
					if((*markers)[l]->isEnabled() && !(*markers)[l]->isFlagged()){
						cnt_dout++;
					}
				}
			}

			opts::printLog("Overall results for chrom " + getString<int>(run_chr) + ", " + getString<int>(cnt_out) + " SNPs slated for removal, " + getString<int>(cnt_in) + " remaining\n");
			opts::printLog("For chrom " + getString<int>(run_chr) + ", " + getString<int>(cnt_dis) + " SNPs already removed or disabled, " + getString<int>(cnt_dout) + " more SNPs slated to be chopped out, " + getString<int>(cnt_din) + " would remain\n");
		}// next chrom

		if(options.doLDchop()){
			for(int i = 0; i < (int)include.size(); i++){
				if(!include[i] && (*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
					(*markers)[i]->setEnabled(false);
					orig_num_markers++;
				}
			}
		}

		if(LDIN.is_open()){
			LDIN.close();
		}
		if(LDOUT.is_open()){
			LDOUT.close();
		}
	}

}

/*Borrowed and adapted from Plink*/
vector<double> LD::correlation2SNP(Marker* pri, Marker* sec){
	vector<double> results(2,-1);
	if(pri->isMicroSat() || sec->isMicroSat()){
		return results;
	}

	double X = 0; //sum count of primary minor allele
	double X2 = 0; // sum of pri squares
	double Y = 0; // sum count of secondary minor allele
	double Y2 = 0; // sec. MA squared sum
	double XY = 0; // multiplicative sum of minor allele counts
	double count = 0;

	int ploc = pri->getLoc();
	int sloc = sec->getLoc();

	for(int i = 0; i < (int)samples->size(); i++){
		Sample* samp = (*samples)[i];
		if(!samp->isFounder()){
			continue;
		}

		bool p1 = samp->getAone(ploc);
		bool p2 = samp->getAtwo(ploc);
		bool pM = samp->getAmissing(ploc);
		if(p1 && p2 && pM) continue;

		bool q1 = samp->getAone(sloc);
		bool q2 = samp->getAtwo(sloc);
		bool qM = samp->getAmissing(sloc);
		if(q1 && q2 && qM) continue;

		count++;

		int sp = 0, ss = 0;
		if(pri->getChrom() == opts::_CHRX_ || sec->getChrom() == opts::_CHRX_){
			if(pri->getChrom() == opts::_CHRX_ && samp->getSex()){
				if(!p1){
					sp = 1;
				}
			}
			else{
				if(!p1){
					if(!p2){
						sp = 2;
					}
					else{
						sp = 1;
					}
				}
			}
			if(sec->getChrom() == opts::_CHRX_ && samp->getSex()){
				if(!q1){
					ss = 1;
				}
			}
			else{
				if(!q1){
					if(!q2){
						ss = 2;
					}
					else{
						ss = 1;
					}
				}
			}
		}
		else{
			if(!p1){
				if(!p2){
					sp = 2;
				}
				else{
					sp = 1;
				}
			}
			if(!q1){
				if(!q2){
					ss = 2;
				}
				else{
					ss = 1;
				}
			}
		}
		X += sp;
		Y += ss;
		XY += sp*ss;
		sp *= sp;
		ss *= ss;
		X2 += sp;
		Y2 += ss;
	}

	X /= count;
	X2 /= count;
	Y /= count;
	Y2 /= count;
	XY /= count;

	double var1 = X2 - X*X;
	double var2 = Y2 - Y*Y;
	double D = XY - X*Y;
	results[0] = (D*D)/(var1*var2);

	double dmax1 = (1-X)*Y;
	double dmax2 = (X)*(1-Y);

	double dmax = dmax1 < dmax2 ? dmax1 : dmax2;
	if(dmax == 0){
		results[1] = -1;
	}
	else{
		results[1] = D/dmax;
	}

	return results;
}

void LD::calcLDStatistics(){
	map<int, vector<int> > chroms;
	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < (int)markers->size(); m++){
		if((*markers)[m]->getChrom() <= opts::_CHRX_ && isValidMarker((*markers)[m], &options, prev_base, prev_chrom)){
			if((*markers)[m]->isEnabled() && !(*markers)[m]->isFlagged()){
				chroms[(*markers)[m]->getChrom()].push_back(m);
			}
		}
	}

	ofstream LD;
    string f = opts::_OUTPREFIX_ + "ld_calc" + options.getOut() + ".txt";
	if(!overwrite){
		f += "." + getString<int>(order);
	}
	LD.open(f.c_str(),ios::out);
	LD << "CHR1\tSNP1\tCHR2\tSNP2\tr^2\tD'\n";

	map<int, vector<int> >::iterator citer;
	cout << "Chrom size = " << chroms.size() << endl;
	for(citer = chroms.begin(); citer != chroms.end(); citer++){
		vector<int> snps = citer->second;
		for(int m = 0; m < (int)snps.size(); m++){
			Marker* primary = (*markers)[snps[m]];
			int s1 = m + 1;
			int s2 = m + options.getLDWin() - 1;
			if(s2 >= (int)snps.size()){
				s2 = snps.size() - 1;
			}

			while(s1 <= s2){
				Marker* secondary = (*markers)[snps[s1]];
				if((secondary->getBPLOC() - primary->getBPLOC()) > (options.getLDWinKB() * 1000)){
					break;
				}
				vector<double> results = correlation2SNP(primary, secondary);
				LD << primary->getChrom() << "\t" << primary->getProbeID() << "\t" << secondary->getChrom() << "\t" << secondary->getProbeID() << "\t" << results[0] << "\t" << results[1] << endl;
				s1++;
			}
		}
	}
	if(LD.is_open()){
		LD.close();
	}
}

}

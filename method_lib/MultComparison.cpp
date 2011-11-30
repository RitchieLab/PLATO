//liberally pulled from Plink and adapted to work in Plato

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include "MultComparison.h"
namespace Methods{
// ## QQ plots: concentration bands

// > n <- 5000
// > j <- 1:n
// > x <- (j-0.5)/n
// > d <- sort( runif(n) )
// > plot( -log10(x) , -log10(d) , ylim = range(0,8) , xlim=range(0,8) )
// > lines(-log10(x),-log10(qbeta( 0.05 , j , ( n - j + 1 ) , ncp = d ) ),col="blue")
// > lines(-log10(x),-log10(qbeta( 0.95 , j , ( n - j + 1 ) , ncp = d ) ),col="blue")



/*void pprint(ofstream & MT, double p)
{
  MT << setw(10);

  if ( p == 0 )
    MT << "INF";
  else if( p < 0 )
    MT << "NA";
  else
    {
      if ( par::logscale )
	{
	  double lp = -log10(p);
	  if ( lp <= 0 ) lp = 0;
	  MT << lp;
	}
      else
	MT << p;
    }
  MT << " ";
}*/

void MultComparison::calculate(vector<double> & chi, vector<int> & tcnt)
{

  if (tcnt.size() > 0 && tcnt.size() != chi.size() )
     throw MethodException("Internal problem in multiple comparison routine");

  //printLOG("Computing corrected significance values (FDR, Sidak, etc)\n");

  bool altern_pval = tcnt.size() > 0;

  vector<Pair> sp;
  vector<Pair> schi;
  for (unsigned int l=0; l<chi.size(); l++)
    {
      if ( chi[l]>=0 )
	{
	  double p = altern_pval ? Helpers::pT(sqrt(chi[l]),tcnt[l]) : Helpers::chiprobP(chi[l],1);

	  if (p > -1)
	    {
	      Pair pt;
	      pt.p = p;
	      pt.l = l;
	      sp.push_back(pt);

	      Pair b;
	      b.p = chi[l];
	      b.l = altern_pval ? (int)tcnt[l] : 0;
	      schi.push_back(b);
	    }
	}
    }

  if (schi.size()==0)
    {
      cout << "Zero valid tests computed -- no adjusted values calculated\n";
      return;
    }

  double t = (double)sp.size();
  int ti = sp.size();

  //removed to maintain genomic order
  sort(sp.begin(),sp.end());
  sort(schi.begin(),schi.end());

  for(int i = 0; i < (int)sp.size(); i++){
	  mapping[sp[i].l] = i;
  }



  // Genomic control, based on median chi-square

  double lambda;
  double lambda_mean = 0;
  if (sp.size() % 2 == 0 )
    {
      lambda = ( schi[ ti / 2 - 1 ].p + schi[ ti / 2 ].p ) / 2 ;
    }
  else
    {
      lambda = schi[ (ti-1) / 2 ].p ;
    }

  for (unsigned int i=0; i<schi.size(); i++)
    lambda_mean += schi[i].p;
  lambda_mean /= schi.size();

  lambda /= 0.456;
  if (lambda < 1) lambda = 1.00;

  if ( options.doFixedLambda() )
    {
      lambda = options.get_fixed_lambda();
      cout << "Genomic inflation factor (fixed value) is "+getString<double>(lambda)+"\n";
    }
  else
    cout << "Genomic inflation factor (based on median chi-squared) is "
	     +getString<double>(lambda)+"\n";

  cout << "Mean chi-squared statistic is "+getString<double>(lambda_mean)+"\n";
  cout << "Correcting for "+getString<int>(ti)+" tests\n";

  // Consider each test
  // Bonferroni correction

  pv_GC.resize(sp.size());
  pv_sidakSS.resize(sp.size());
  pv_sidakSD.resize(sp.size());
  pv_holm.resize(sp.size());
  pv_BH.resize(sp.size());
  pv_BY.resize(sp.size());
  pv_bonf.resize(sp.size());

  // Genomic control (reverse order)
  int i2=0;
  for (int i=ti-1;i>=0;i--)
    {
      pv_GC[i2++] = altern_pval ? Helpers::pT(sqrt( schi[i].p / lambda ), schi[i].l ) : Helpers::chiprobP( schi[i].p / lambda, 1 );
    }

  // Base adjust values on GC p-values?
 /// if (par::use_GC)
 ///   {
 ///     printLOG("Using genomic-controlled p-values for adjusted p-values\n");
 ///     for (int i=0; i<ti; i++)
///	sp[i].p = pv_GC[i];
///    }

  // Holm
  pv_holm[0] = sp[0].p*t > 1 ? 1 : sp[0].p*t;
  for (int i=1;i<t;i++)
    {
      double x = (ti-i)*sp[i].p < 1 ? (ti-i)*sp[i].p : 1;
      pv_holm[i] = pv_holm[i-1] > x ? pv_holm[i-1] : x;
    }

  // Sidak SS
  for (int i=0;i<t;i++)
    pv_sidakSS[i] = 1 - pow( 1 - sp[i].p , t );


  // Sidak SD
  pv_sidakSD[0] = 1 - pow( 1 - sp[0].p , t );
  for (int i=1;i<t;i++)
    {
      double x = 1 - pow( 1 - sp[i].p , t - i  );
      pv_sidakSD[i] = pv_sidakSD[i-1] > x ? pv_sidakSD[i-1] : x ;
    }

  // BH
  pv_BH[ti-1] = sp[ti-1].p;
  for (int i=ti-2;i>=0;i--)
    {
      double x = (t/(double)(i+1))*sp[i].p < 1 ? (t/(double)(i+1))*sp[i].p : 1 ;
      pv_BH[i] = pv_BH[i+1] < x ? pv_BH[i+1] : x;
    }

  // BY
  double a = 0;
  for (double i=1; i<=t; i++)
    a += 1/i;

  pv_BY[ti-1] = a * sp[ti-1].p < 1 ? a * sp[ti-1].p : 1 ;

  for (int i=ti-2;i>=0;i--)
    {
      double x = ((t*a)/(double)(i+1))*sp[i].p < 1 ? ((t*a)/(double)(i+1))*sp[i].p : 1 ;
      pv_BY[i] = pv_BY[i+1] < x ? pv_BY[i+1] : x;
    }


  for (int l=0; l<t; l++)
    {
      // Bonferroni, etc
      pv_bonf[l] = sp[l].p*t > 1 ? 1 : sp[l].p*t;
    }

}
}

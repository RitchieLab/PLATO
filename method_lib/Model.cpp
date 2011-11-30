//liberally taken from Plink and adapted to work in the Plato environment

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include "Model.h"
#include "Options.h"
#include "Helpers.h"
//#include "phase.h"
//#include "stats.h"

#include <cmath>
namespace Methods{

	void Model::nullCoefs(){
		for(int i = 0; i < (int)coef.size(); i++){
			coef[i] = 0;
		}
	}

    void Model::display(vector<vector<double> > & m)
		    {
				            cout << "\n";
							              for (int i=0; i< (int)m.size(); i++)
											                        {
																		                                cout << i << ")\t";
																										                                      for (int j=0; j<(int)m[i].size(); j++)
																																				                                                cout << m[i][j] << " ";
																																			                                              cout << "\n";
																																														                                                  }
										                  cout << "\n";
														      }

	    void Model::display(vector<double> & m)
			    {
					            cout << "\n";
								              for (int i=0; i< (int)m.size(); i++)
												                        cout << i << ")\t" << m[i] << "\n";
											                  cout << "\n";
															                    cout << "\n";
																				    }

		    void Model::display(vector<int> & m)
				    {
						            cout << "\n";
									              for (int i=0; i< (int)m.size(); i++)
													                        cout << i << ")\t" << m[i] << "\n";
												                  cout << "\n";
																                    cout << "\n";
																					    }

			Model::Model()
{
				isnull = false;

  np = nind = 0;
  haploid.resize(0);
  xchr.resize(0);
  order.clear();
  sex_effect = false;
  all_valid = true;
  has_snps = true;
  testParameter = 1; // Permutation test parameter

  // Automatically add intercept now
  label.push_back("M"); // Intercept
  type.push_back( INTERCEPT );
  order.push_back(0);


  ///////////////////////////
  // Default additive coding

  mAA = 0;
  mAB = 1;
  mBB = 2;


  ///////////////////////////
  // Set X chromosome coding

  if ( options.get_xchr_model() == 1 )
    {
      mA = 0;
      mB = 1;
    }
  else if ( options.get_xchr_model() == 2 )
    {
      mA = 0;
      mB = 2;
    }
  else if ( options.get_xchr_model() > 2 )
    {
      mA = 0;
      mB = 1;
    }

}

void Model::setDominant()
{
  mAA = 0;
  mAB = 1;
  mBB = 1;
  mA = 0;
  mB = 1;
}

void Model::setRecessive()
{
  mAA = 0;
  mAB = 0;
  mBB = 1;

  // No haploid effect
  mA = mB = 0;
}

void Model::addSexEffect()
{
  sex_effect = true;
  type.push_back( SEX );
  order.push_back(0);
}

bool Model::isSexInModel()
{
  return sex_effect;
}

void Model::hasSNPs(bool b)
{
  has_snps = b;
}

void Model::setMissing()
{

  // Fill in missing data with existing pattern
  // and also optional per-test missingness

	/////need to add PHENOTYPE missing stuff
  miss.clear();
  miss.resize(data_set->num_inds(), false);//P->n,false);
  for (int i=0; i<data_set->num_inds(); i++){//P->n; i++)
		bool pheno_miss = false;
	  if (options.getUsePheno())
	  {
			int index = options.getPhenoLoc();

			if (options.getPhenoName() != "") {
				index = data_set->get_trait_index(options.getPhenoName());
			}
			if (index < 0) {
				throw MethodException(
						"Internal Error: Trait/Phenotype index value < 0 in Linear Regression!");
			}
			if(data_set->get_sample(i)->getPheno(index) == options.getPhenoMissing()){
				pheno_miss = true;
			}
//			Y.push_back(data_set->get_sample(i)->getPheno(index));//data_set->get_trait_index(options.getTraits().at(0))));//P->sample[i]->pperson->phenotype ) ;
		}
	  else
	  {
			if(data_set->get_sample(i)->getPheno() == options.getPhenoMissing()){
				pheno_miss = true;
			}
//			Y.push_back(data_set->get_sample(i)->getPheno());//P->sample[i]->pperson->phenotype ) ;
	  }

	  if(!data_set->get_sample(i)->isEnabled() || pheno_miss ){
		  miss[i] = true;// || data_set->get_sample(i)->isFlagged()) miss[i] = true;
	  }
  }
	  //    if ( P->sample[i]->missing || P->sample[i]->missing2 ) miss[i] = true;

}

vector<bool> Model::getMissing()
{
  return miss;
}

void Model::yokeMissing(Model * m)
{
  //
}

void Model::setMissing(vector<bool> & include)
{

  // Fill in missing data with existing pattern
  if ( (int)include.size() != data_set->num_inds())//P->n )
    throw MethodException("A problem in setMissing()\n");

  miss.resize(data_set->num_inds(), false);//P->n,false);
  for (int i=0; i<data_set->num_inds(); i++){
		bool pheno_miss = false;
	  if (options.getUsePheno()) {
			int index = options.getPhenoLoc();

			if (options.getPhenoName() != "") {
				index = data_set->get_trait_index(options.getPhenoName());
			}
			if (index < 0) {
				throw MethodException(
						"Internal Error: Trait/Phenotype index value < 0 in Linear Regression!");
			}
			if(data_set->get_sample(i)->getPheno(index) == options.getPhenoMissing()){
				pheno_miss = true;
			}
//			Y.push_back(data_set->get_sample(i)->getPheno(index));//data_set->get_trait_index(options.getTraits().at(0))));//P->sample[i]->pperson->phenotype ) ;
		} else {
			if(data_set->get_sample(i)->getPheno() == options.getPhenoMissing()){
				pheno_miss = true;
			}
//			Y.push_back(data_set->get_sample(i)->getPheno());//P->sample[i]->pperson->phenotype ) ;
		}

	    if(!data_set->get_sample(i)->isEnabled() || !include[i] || pheno_miss) miss[i] = true;
	  }
	  //if ( P->sample[i]->missing || ! include[i] ) miss[i] = true;

}


void Model::addAdditiveSNP(Marker* a)
{

  if ( ! has_snps )
    throw MethodException("Cannot add SNP to this MODEL");

  additive.push_back(a);

  if ( opts::_CHRX_ == a->getChrom())//data_set->get_locus(a)->getChrom())//par::chr_sex[P->locus[a]->chr] )
  {  xchr.push_back(true); }//cout << "X CHR FOUND!!!!\n";}
  else
    xchr.push_back(false);

  if ( opts::_CHRXY_ == a->getChrom())//data_set->get_locus(a)->getChrom())//par::chr_haploid[P->locus[a]->chr] )
    haploid.push_back(true);
  else
    haploid.push_back(false);

  type.push_back( ADDITIVE );
  order.push_back( additive.size() - 1 );

}

void Model::addDominanceSNP(Marker* d)
{
  if ( ! has_snps )
    throw MethodException("Cannot add SNP to thie MODEL");

  dominance.push_back(d);

  type.push_back( DOMDEV );
  order.push_back( dominance.size() - 1 );

}

void Model::addCovariate(int c)
{
//	cout << "Adding Covariate: " << getString<int>(c) << "\n";
  covariate.push_back(c);

  type.push_back( COVARIATE );
  order.push_back( covariate.size() - 1 );

}

void Model::addHaplotypeDosage(set<int> & h)
{
  haplotype.push_back(h);

  type.push_back( HAPLOTYPE );
  order.push_back( haplotype.size() - 1 );

}

void Model::addInteraction(int a, int b)
{
  int2 i;
  i.p1 = a;
  i.p2 = b;
  interaction.push_back(i);

  type.push_back( INTERACTION );
  order.push_back( interaction.size() - 1 );

}

void Model::buildDesignMatrix()
{

  // Build X matrix (including intercept)
  // Iterate a person at a time, entering only
  // valid rows into X (non-missing); also build Y
  // at the same time

///  if ( has_snps && par::SNP_major )
///    error("Internal error: must be individual-major to perform this...\n");


  ///////////////////////
  // Number of parameters

  // Standard variables
  // Note: 'additive' really means 'main effect' here,
  //       i.e. which can also be coded recessive or dominant
  //       i.e. the distinction is between the 2df general model

  np = 1
    + additive.size()
    + dominance.size()
    + haplotype.size()
    + covariate.size()
    + interaction.size();

//  cout << "np1: " << getString<int>(np) << "\n";

  // Sex effect?
  if ( sex_effect )
    np++;

//  cout << "np2: " << getString<int>(np) << "\n";


  // QFAM variables
  //  if (par::QFAM_total) np++;
  //  else if (par::QFAM_between || par::QFAM_within1 || par::QFAM_within2) np +=2;

  if (options.getQFAM_total()
      || options.getQFAM_between()
      || options.getQFAM_within1()
      || options.getQFAM_within2())
    {
      np++;

      type.push_back( QFAM );
      order.push_back( 0 );
    }

//  cout << "np3: " << getString<int>(np) << "\n";



  ///////////////////////////
  // Consider each individual
  //cout << "data_set->num_inds() = " << getString<int>(data_set->num_inds()) << "\n";

  for (int i=0; i < data_set->num_inds(); i++)//P->n; i++)
  {

	  Sample * person = data_set->get_sample(i);//P->sample[i];

	  // Ignore if missing phenotype, or the user set this to missing
	  if ( miss[i] )
	  {
		  //cout << "MISSING PHENOTYPE IN buildDesignMatrix()" << endl;
		  continue;
	  }
	  //cout << "NOT missing phenotype in buildDesignMatrix()" << endl;

	  /////////////////////////////
	  // 0) Intercept
	  // 1) Main effects of SNPs
	  // 2) Dominance effects of SNPs
	  // 3) Haplotypes
	  // 4) Covariates
	  // 5) Interactions of the above
	  // 6) QFAM variables

	  // Populate this vector with terms for this
	  // individual

	  skip = false;

	  vector<double> trow(np);

	  //      cout << "np final: " << getString<int>(np) << "\n";
	  //cout << "trow for individual number " << getString<int>(i) << " ";
	  for (int p = 0; p < np; p++)
	  {

		  int pType = type[p];
		  //	  if(isnull && pType == INTERCEPT){
		  //		  trow[p] = buildIntercept();
		  //	  }
		  //	  else{
		  switch ( pType )
		  {
		  case INTERCEPT :
			  trow[p] = buildIntercept();
			  //cout << "Intercept: " << getString<int>(trow[p]) << "    ";
			  break;
		  case ADDITIVE :
			  trow[p] = buildAdditive( person, order[p] );
			  //cout << "Additive: " << getString<int>(trow[p]) << endl;
			  break;
		  case DOMDEV :
			  trow[p] = buildDominance(person, order[p] );
			  //cout << "Domdev: " << getString<int>(trow[p]) << "    ";
			  break;
		  case HAPLOTYPE :
			  trow[p] = buildHaplotype(i, order[p] );
			  //cout << "Haplotype: " << getString<int>(trow[p]) << "    ";
			  break;
		  case SEX :
			  trow[p] = buildSex(person);
			  //cout << "Sex: " << getString<int>(trow[p]) << "    ";
			  break;
		  case COVARIATE :
			  trow[p] = buildCovariate( person, order[p] );
			  //cout << "Covariate: " << getString<int>(trow[p]) << "    ";
			  break;
		  case INTERACTION :
			  trow[p] = buildInteraction( person, order[p], trow );
			  //cout << "Interaction: " << getString<int>(trow[p]) << endl;
			  break;
		  case QFAM :
			  trow[p] = buildQFAM( person );
			  //cout << "QFAM: " << getString<int>(trow[p]) << "    ";
			  break;
		  }
		  //cout << "     " << trow[p] << endl;
		  //	  }
	  }

	  //cout << "\n";

	  if (skip)
	  {
		  miss[i] = true;
		  skip = false;
		  continue;
	  }


	  ////////////////////////////
	  // Add row to design matrix
	  //cout << "Adding row to desing matrix" << endl;
	  X.push_back(trow);

  }


  /////////////////////////////////////////////////
  // Set number of non-missing individuals

  nind = X.size();

//  cout << "nind: " << getString<int>(nind) << "\n";

  //////////////////////////////////
  // Apply a parameter list filter?

///  if ( par::glm_user_parameters )
///    {

      // Intercept always fixed in

///      np = 1;
///      vector<string> label2 = label;
///      label.clear();
///      label.push_back("M");
///      int np2 = label2.size();

///      for (int i=0; i< par::parameter_list.size(); i++)
///	{
///	  if ( par::parameter_list[i] >= 1
///	       && par::parameter_list[i] < np2 )
///	    {
///	      np++;
///	      label.push_back(label2[ par::parameter_list[i] ]);
///	    }
///	}


      // For each individual

///      for ( int i = 0 ; i < X.size() ; i++)
///	{

///	  vector<double> X2(1);

///	  X2[0] = 1;

///	  for ( int j = 0 ; j < par::parameter_list.size() ; j++)
///	    {
///	    if ( par::parameter_list[j] >= 1
///		 && par::parameter_list[j] < np2 )
///	      {
///		X2.push_back(X[i][ par::parameter_list[j] ]);
///	      }
///	    }
///	  X[i] = X2;
///	}
///    }


  /////////////////////////////////////////
  // VIF-based check for multicollinearity

  //cout << "All Valid before checkVIF: " << getString<bool>(all_valid) << "\n";

  all_valid = checkVIF();

 // cout << "All Valid after checkVIF: " << getString<bool>(all_valid) << "\n";


  ///////////////////////
  // Add Y variable also

  setDependent();

  // Now we are ready to perform the analysis

///  if (par::verbose)
///    {
//      cout << "X design matrix\n";
//      display(X);
//      cout << "\n";
///    }

}


vector<bool> Model::validParameters()
{

  // Empty model?
  if (np==0 || nind==0)
    {
      vector<bool> v(np,false);
      all_valid = false;
      return v;
    }


  // Display covariance matrix in verbose mode

///  if ( par::verbose )
///    {
//      cout << "Covariance matrix of estimates\n";
//      display(S);
//     cout << "\n";
///    }

  // Check for multicollinearity

  // For each term, see that estimate is not too strongly (r>0.99)
  // correlated with another, starting at last

  valid.resize(np);

  for (int i = 1; i<np; i++)
    {
      valid[i] = true;
      if ( S[i][i] < 1e-20 ) { valid[i] = all_valid = false; }
      else if ( ! Helpers::realnum(S[i][i]) ) { valid[i] = all_valid = false; }
    }

  if ( all_valid )
    for (int i = np-1; i>0; i--)
      {
	for (int j = i-1; j>=0; j--)
	  {
	    if ( S[i][j] / sqrt( S[i][i] * S[j][j] ) > 0.99999 )
	      {
		valid[i] = false;
		all_valid = false;
		break;
	      }
	  }
      }
  return valid;
}


double Model::getStatistic()
{
  if (all_valid)
    {
      return ( coef[testParameter] * coef[testParameter] )
	/ S[testParameter][testParameter];
    }
  else return 0;
}


// **********************************************
// *** Function to aid testing linear
// *** hypotheses after estimating of a
// *** regression
// **********************************************

double Model::linearHypothesis(vector<vector<double> > & H, vector<double> & h)
{

  // (H v - h)' (H S H')^-1 (H b - h) ~ X^2 with j df
  // where H = constraint matrix (j x (p+1)
  //       h = null coefficient values
  //       S = estimated covariance matrix of coefficients

  //   return ( H*b - h ).transpose()
  //     * ( H*v*H.transpose() ).inverse()
  //     * ( H*b - h ) ;

  int nc = h.size(); // # of constraints

  // 1. Calculate Hb-h

  vector_t outer;
  outer.resize(nc,0);

  for (int r = 0; r < nc; r++)
    for (int c = 0; c < np; c++)
      outer[r] += H[r][c] * coef[c];

  for (int r = 0; r < nc; r++)
    outer[r] -= h[r];

  // 2. Calculate HVH'

  vector<vector<double> > tmp;
  Helpers::sizeMatrix(tmp,nc,np);

  for (int r = 0; r < nc; r++)
    for (int c = 0; c < np; c++)
      for (int k = 0; k < np; k++)
	tmp[r][c] += H[r][k] * S[k][c];

  vector<vector<double> > inner;
  Helpers::sizeMatrix(inner,nc,nc);

  for (int r = 0; r < nc; r++)
    for (int c = 0; c < nc; c++)
      for (int k = 0; k < np; k++)
	inner[r][c] += tmp[r][k] * H[c][k];

  inner = Helpers::svd_inverse(inner);

  vector<double> tmp2;
  tmp2.resize(nc,0);

  for (int c = 0; c < nc; c++)
    for (int k = 0; k < nc; k++)
      tmp2[c] += outer[k] * inner[k][c];

  double result = 0;

  for (int r = 0; r < nc; r++)
    result += tmp2[r] * outer[r];

  return result;

}


bool Model::checkVIF()
{

  // Calculate correlation matrix for X
  // Skip intercept

  int p = X.size();
  if (p<2) return false;

  int q = X[0].size() - 1;
  if ( q < 2 ) return true;

  vector<double> m(q);
  vector<vector<double> > c;
  Helpers::sizeMatrix(c,q,q);

  for (int i=0; i<p; i++)
    for (int j=0; j<q; j++)
      m[j] += X[i][j+1];

  for (int j=0; j<q; j++)
    m[j] /= (double)p;

  for (int i=0; i<p; i++)
    for (int j1=0; j1<q; j1++)
      for (int j2=j1; j2<q; j2++)
	c[j1][j2] += ( X[i][j1+1] - m[j1] ) * ( X[i][j2+1] - m[j2] );

  for (int j1=0; j1<q; j1++)
    for (int j2=j1; j2<q; j2++)
      c[j1][j2] /= (double)(p-1);

  for (int j1=0; j1<q; j1++)
    for (int j2=j1+1; j2<q; j2++)
      {
	c[j1][j2] /= sqrt( c[j1][j1] * c[j2][j2] );
	c[j2][j1] = c[j1][j2];

 	if ( c[j2][j1] > 0.999 )
 	  {
//	    if (par::verbose)
	      //cout << "individual element > 0.999, skipping VIF\n";
	    return false;
	  }
      }


  // Any item with zero variance?

  for (int j=0; j<q; j++)
    {
      if ( c[j][j] == 0 || ! Helpers::realnum( c[j][j] ) )
      {
    	  //cout << "found item with zero variance...VIF returning FALSE\n";
    	  return false;
      }
      c[j][j] = 1;
    }

  // Get inverse
  c = Helpers::svd_inverse(c);

///  if (par::verbose)
///    {
//      cout << "VIF on diagonals\n";
 //     display(c);
 //     cout << "\n";
///    }

  // Calculate VIFs
////  double maxVIF = 0;
  for (int j=0;j<q;j++)
    {

      // r^2 = 1 - 1/x where x is diagonal element of inverted
      // correlation matrix
      // As VIF = 1 / ( 1 - r^2 ) , implies VIF = x

      if ( c[j][j] > options.get_vif_threshold() )
      {
    	  //cout << "found item with vif > threshold, returning false...\n";
    	  return false;
      }
    }


  return true;

}




double Model::buildIntercept()
{
  return 1;
}

double Model::buildAdditive(Sample * person , int snp )
{

  // Additive effects (assuming individual-major mode)
if(isnull)
{
	return 0;
}

  int s = additive[snp]->getLoc();

  bool i1 = person->getAone(s);
  bool i2 = person->getAtwo(s);
  bool i3 = person->getAmissing(s);

  if ( xchr[snp] )
	{

	  /////////////////////////
	  // X chromosome coding

	  if ( person->getSex() ) // male
	  {
			if ( i1 )
			{
				if (  i2 && i3 )
				{
				  skip = true;
				  return 0;
				}
				else if(i2 && !i3)
					return mA;
			}
			else
			{
				if ( i2 )
				{
					// This should not happen...
					skip = true;
					return 0;
				}
				else
					return mB;
			}
		}
	  else // female x-chromosome
	{
		  if ( i1 )
		  {
			  if (  i2 && i3 )
			  {
				  skip = true;
				  return 0;
			  }
			  else if(i2 && !i3)
				  return mAA;
		  }
		  else
		  {
			  if ( i2 )
				  return mAB; // het
			  else
				  return mBB; // hom
		  }
	}

	}
  else if ( haploid[snp] )
  {

      ///////////////////
      // Haploid coding

	  if ( i1 )
	  {
		  if ( i2 && i3 )
		  {
			  skip = true;
			  return 0;
		  }
		  else if(i2 && !i3)
			  return 0;
	  }
	  else
	  {
		  if ( i2 )
		  {
			  // haploid het
			  skip = true;
			  return 0;
		  }
		  else
			  return 1;
	  }

    }
  else
  {

	///////////////////////
	// Autosomal coding

	  if ( i1 )
	  {
		  if ( i2 && i3 )
		  {
			  skip = true;
			  return 0;
		  }
		  else if(i2 && !i3)
			  return mAA;
	  }
      else
      {
    	  if ( i2 )
    		  return mAB; // het
    	  else
    		  return mBB; // hom
      }

  }

  return 0;
}


double Model::buildDominance(Sample * person, int snp)
{

  ////////////////////
  // Dominance effects

  int s = dominance[snp]->getLoc();
  bool i1 = person->getAone(s);
  bool i2 = person->getAtwo(s);
  bool i3 = person->getAmissing(s);

  if ( i1 )
    {
      if ( i2  && i3 )
	{
	  skip = true;
	  return 0;
	}
      else if(i2 && !i3)
	return 0;
    }
  else
    {
      if ( i2 )
	return 1; // het
      else
	return 0; // hom
    }

  return 0;
}

double Model::buildHaplotype(int i, int h )
{

  ////////////////////
  // Haplotype dosage

///  if ( P->haplo->include[i] )
///    return P->haplo->dosage(i,haplotype[h]);

  // No valid haplotypes for this person
  skip = true;
  return 0;

}

double Model::buildSex(Sample * person )
{

  ////////////////////////////////////
  // Sex effect (automatically set for
  // X chromosome models)

  if ( person->getSex() )
    return 1;
  else
    return 0;

}

double Model::buildCovariate(Sample * person, int j)
{

  /////////////
  // Covariates

//	if(isnull){

//	}
	//cout << "index j: " << getString<int>(j) << "\n";
//	cout << "Building Covariate: " << getString<double>(person->getCovariate(covariate[j])) << "\n";
  return person->getCovariate(covariate[j]);//person->clist[covariate[j]];

}

double Model::buildInteraction(Sample * person, int j, vector<double> & trow )
{
  ///////////////
  // Interactions
if(isnull){
	return 0;
}
  return trow[ interaction[j].p1 ] * trow[ interaction[j].p2 ];

}

double Model::buildQFAM(Sample * person)
{

  ///////////////
  // QFAM

  if ( options.getQFAM_total() )
    return person->getT();
  else if ( options.getQFAM_between() )
    return person->getFamily()->getB();
  else if ( options.getQFAM_within1() || options.getQFAM_within2() )
    return person->getW();
  else
    throw MethodException("Internal problem with QFAM model specification");
}
}

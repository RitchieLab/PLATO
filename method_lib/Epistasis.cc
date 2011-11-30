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
*File: Epistasis.cc
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
#include "Epistasis.h"
#include "LogisticRegression.h"
#include "LinearRegression.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void Epistasis::FilterSummary(){
}

void Epistasis::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

void Epistasis::filter(){
}

void Epistasis::evalBioFile(ofstream &EPI)
{
    string fname = opts::_OUTPREFIX_ + "epistasis" + options.getOut() + ".txt";
    if(!overwrite)
    {
        fname += "." + getString<int>(order);
    }

    string logname = fname + ".log";
	ofstream epi_log(logname.c_str());
	if(!epi_log)
	{
		throw MethodException("Cannot open " + logname + " for writing!");
	}
	int nepi = 0;

	// Only output epistatic tests that have p < par::epi_alpha1;
	// Do not even attempt to save any epistatic results -- go straight to STDOUT

	// Also present summary results for all epi1 SNPs
	// (i.e. average / proportion of significant epistatic tests
	//  at a certain alpha level, par::epi_alpha2

	int nl_all = data_set->num_loci();
	vector<bool> sA(nl_all, false);
	vector<bool> sB(nl_all, false);

	vector<int> good_indexes = Helpers::findValidMarkersIndexes(data_set->get_markers(), &options);

	options.readBioTextFile(options.getBioSnpFile());
	map<string, vector<string> > pairs = options.getBioPairs();


	double epi_alpha1 = fabs(Helpers::ltqnorm(options.get_epi_alpha1() / 2));

	vector<int> summary_sig(nl_all, 0);
	vector<int> summary_good(nl_all, 0);
	vector<double> best_score(nl_all, 0);
	vector<int> best_partner(nl_all);

	LinearRegression lr;
	lr.setOptions(options);
	lr.resetDataSet(data_set);

	LogisticRegression logr;
	logr.set_parameters(&options);
	logr.resetDataSet(data_set);

	//begin iterating over bio sets

	int ssize = data_set->num_inds();
	int epcc = 0;
	map<string, vector<string> >::iterator bio_iter;
	for(bio_iter = pairs.begin(); bio_iter != pairs.end(); bio_iter++)
	{
		string snp1 = bio_iter->first;
		vector<string> snp_list = bio_iter->second;

		Marker* s1 = NULL;
		try
		{
			int s1_loc = data_set->get_locus_index(snp1);
			s1 = data_set->get_locus(s1_loc);
		}
		catch(MethodException & ex)
		{
			for(int i = 0; i < (int)snp_list.size(); i++)
			{
				epi_log << snp1 << " " << snp_list[i] << " ---> " << snp1 << " not found!" << endl;
			}
			continue;
		}

		if(s1 != NULL)
		{
			if(!s1->isEnabled())
			{
				for(int i = 0; i < (int)snp_list.size(); i++)
				{
					epi_log << snp1 << " " << snp_list[i] << " ---> " << snp1 << " is disabled!" << endl;
				}
				continue;
			}
			for(int i = 0; i < (int) snp_list.size(); i++)
			{
				Marker* s2 = NULL;
				try
				{
					int s2_loc = data_set->get_locus_index(snp_list[i]);
					if(s2_loc == s1->getLoc())
					{
						epi_log << snp1 << " " << snp_list[i] << " ---> " << "Snps have same index!" << endl;
						continue;
					}
					s2 = data_set->get_locus(s2_loc);
				}
				catch(MethodException & ex)
				{
					epi_log << snp1 << " " << snp_list[i] << " ---> " << snp_list[i] << " not found!" << endl;
				}
				if(s2 != NULL)
				{
					if(!s2->isEnabled())
					{
						epi_log << snp1 << " " << snp_list[i] << " ---> " << snp_list[i] << " is disabled!" << endl;
						continue;
					}

					//perform comparison of bio snps

					if(opts::_CHRX_ == s1->getChrom() || opts::_CHRX_ == s2->getChrom())
						continue;

					// SNPs too close (case-only analysis)
					if (options.doEpiCaseOnly())
						if (s1->getChrom() == s2->getChrom())
							if (fabs((double) (s1->getBPLOC() - s2->getBPLOC()))
									< options.getEpiCaseOnlyKbGap() * 1000)
								continue;

					cout << "Peforming tests of epistasis: group " << ++epcc
							<< "        \r";
					cout.flush();

					//////////////////////////////////
					// Perform test of epistasis here
					if (opts::_BINTRAIT_ && options.doEpiFast())
					{

						double z; // statistic from either method

						// Odds ratio test
						// make two 2x2 tables

						int a11, a12, a21, a22;
						int u11, u12, u21, u22;
						a11 = a12 = a21 = a22 = 0;
						u11 = u12 = u21 = u22 = 0;

						for(int s = 0; s < ssize; s++)
						{
							Sample* samp = data_set->get_sample(s);
							bool b1 = samp->getAone(s2->getLoc());
							bool b2 = samp->getAtwo(s2->getLoc());
							bool b3 = samp->getAmissing(s2->getLoc());
							bool a1 = samp->getAone(s1->getLoc());
							bool a2 = samp->getAtwo(s1->getLoc());
							bool a3 = samp->getAmissing(s1->getLoc());

							if(!samp->isEnabled() || b3 || a3)
							{
								continue;
							}
							if (samp->getAffected())
							{

								if (!b1){
									if (!b2)//   ??x00
									{
										if (!a1){
											if (!a2)
												a11 += 4; // 00 x 00
											else {
												a11 += 2;
												a21 += 2;
											} // 01 x 00
										} else if (a2)
											a21 += 4; // 11 x 00
									} else //   ??x01
									{
										if (!a1){
											if (!a2){
												a11 += 2;
												a12 += 2;
											} // 00 x 01
											else {
												a11++;
												a21++;
												a12++;
												a22++;
											} // 01x01
										} else if (a2){
											a21 += 2;
											a22 += 2;
										} // 11 x 01
									}
								} else if (b2)//*b2) // ?? x 11
								{

									if (!a1){
										if (!a2)
											a12 += 4; // 00 x 01
										else {
											a12 += 2;
											a22 += 2;
										} // 01 x 01
									} else if (a2)
										a22 += 4; // 11 x 01

								}
							}

							// Unaffecteds?
							else if (!options.getEpiCaseOnly())// unaffected
							{
								if (!b1){
									if (!b2)//   ??x00
									{
										if (!a1){
											if (!a2)
												u11 += 4; // 00 x 00
											else {
												u11 += 2;
												u21 += 2;
											} // 01 x 00
										} else if (a2)
											u21 += 4; // 11 x 00
									} else //   ??x01
									{
										if (!a1){
											if (!a2){
												u11 += 2;
												u12 += 2;
											} // 00 x 01
											else {
												u11++;
												u21++;
												u12++;
												u22++;
											} // 01x01
										} else if (a2){
											u21 += 2;
											u22 += 2;
										} // 11 x 01
									}
								} else if (b2)//  ?? x 11
								{

									if (!a1){
										if (!a2)
											u12 += 4; // 00 x 01
										else {
											u12 += 2;
											u22 += 2;
										} // 01 x 01
									} else if (a2)
										u22 += 4; // 11 x 01

								}

							}
						}

						// Calculate log(OR) and SEs

						double or_aff, v_aff, or_unf, v_unf;

						or_aff = log((double) (a11 * a22) / (double) (a12 * a21));
						v_aff = 1 / (double) a11 + 1 / (double) a12 + 1
								/ (double) a21 + 1 / (double) a22;

						// Case-only z-score (if requested)
						if (options.getEpiCaseOnly())
							z = fabs(or_aff / sqrt(v_aff));
						else // Standard case-control analysis
						{
							or_unf = log((double) (u11 * u22)
									/ (double) (u12 * u21));
							v_unf = 1 / (double) u11 + 1 / (double) u12 + 1
									/ (double) u21 + 1 / (double) u22;
							z = fabs((or_aff - or_unf) / sqrt(v_aff + v_unf));
						}

						//////////////////////////////
						// --nop option in effect
						// Just output z score, if valid & above threshold

						if (options.getEpiQuickscan())
						{
							// Is this worth recording?
							if (Helpers::realnum(z))
							{
								nepi++;
							//	cout << "first nepi increased: " << nepi << endl;

								if (z >= epi_alpha1)
									EPI << s1->getChrom() << " "
											<< s1->getRSID() << " "
											<< s2->getChrom() << " "
											<< s2->getRSID()
											<< " " << z * z << "\n";
								EPI.flush();
								continue;
							}
						}

						/////////////////////////////////
						// More full parsing of results

						// Check this is a proper result

						if (options.getEpiFilter() && Helpers::realnum(z))
						{

							// One more test performed
							nepi++;

							// Is this worth recording?

							if (z >= epi_alpha1)
							{
								EPI << s1->getChrom() << " "
								    << s1->getRSID() << " "
										<< s2->getChrom() << " "
										<< s2->getRSID() << " "
										<< " " << z * z << " " << Helpers::normdist(-z) * 2 << " " << "\n";
								EPI.flush();
							}
							else
								continue; // skip to next pair (skip logistic test)

						}
						else if (!options.getEpiFilter())
						{
							// Record all results here, whether NA or otherwise
							EPI << s1->getChrom() << " "
							    << s1->getRSID() << " "
									<< s2->getChrom() << " "
									<< s2->getRSID() << " "
									<< z * z << " "
									<< Helpers::normdist(-z) * 2 << " " << "\n";
							EPI.flush();
						}
						else
							continue; // if bad statistic for this test, do not try logistic

					} // End of binary OR test


					///////////////////////////////////////////////
					// Logistic or linear regression epistasis test

					if (!options.doEpiFast())
					{
						vector<double> b;
						double chisq = 0;
						double F = 0;
						double fstat = 0;
						double r2 = 0;

						if (opts::_BINTRAIT_)
						{
							vector<unsigned int> model;
							model.push_back((unsigned int)s1->getLoc());
							model.push_back((unsigned int)s2->getLoc());
							logr.setFullInteraction(true);
							logr.calculate(model);

							b = logr.getCoefficients();
							vector<double> ses = logr.getCoeffStandardErr();
							if(ses.size() > 0 && b.size() > 0)
							{
								double se = ses[2];
								double Z = b[2] / se;
								chisq = Z*Z;
							}
						}
						else
						{
							//this is not a binary trait, run linear regression...
							lr.reset();
							vector<int> model;
							model.push_back(s1->getLoc());
							model.push_back(s2->getLoc());
							lr.calculate(model);

							b = lr.getCoefs();
							chisq = lr.getStatistic();

							F = lr.findF();
							fstat = lr.getFStat();
							r2 = lr.calculateRSquared();
						}

						bool notvalid = false;
						if(!Helpers::realnum(chisq))
						{
							notvalid = true;
							chisq = 0;
						}

						double pvalue = Helpers::p_from_chi(chisq, 1);
						double z = sqrt(chisq);
						// Is this result worth displaying?

						if (!notvalid)
						{

							// One more valid test performed
							nepi++;
						}
							// Count as a good result

						// Is this result worth displaying?

						if (z >= epi_alpha1)
						{
							EPI << s1->getChrom() << " "
							    << s1->getRSID() << " "
									<< s2->getChrom() << " "
									<< s2->getRSID() << " ";
							if (!notvalid)
							{
								if (opts::_BINTRAIT_)
									EPI << exp(b[2]) << " "
											<< chisq << " " << pvalue
											<< " " << "\n";
								else
									EPI << b[1] << " " << b[2] << " " << b[3] << " "
											<< chisq << " " << pvalue
											<< " " << r2
											<< " " << fstat
											<< " " << F
											<< " " << "\n";
								vector<string>labels = lr.getLabels();

								if(labels.size() > 1 && !opts::_BINTRAIT_)
								{
									vector<double>pvals = lr.getPvalues();
									vector<double>coefs = lr.getCoefs();
									vector<double> vars = lr.getVar();
									vector<double> zs = lr.getZs();

									for(int l = 1; l < (int)labels.size(); l++)
									{
										bool okay = vars[l] < 1e-20 || !Helpers::realnum(vars[l]) ? false : true;
										double se = 0;
										if(okay)
										{
											se = sqrt(vars[l]);
										}
										EPI << s1->getChrom() << "\t" << s1->getRSID() << "\t" << s2->getChrom() << "\t" << s2->getRSID() << "\t";
										EPI << labels[l] << "\t" << "\t" << coefs[l] << "\t" << exp(coefs[l]) << "\t" << se << "\t" << zs[l] << "\t" << pvals[l] << endl;
									}
								}
							}
							else
								EPI << "NA" << " " << "NA"
										<< " " << "NA" << " " << "\n";

							EPI.flush();

						}//end if (z >= epi_alpha1)
					} //end if (!options.doEpiFast())
				}//end if(s2 != NULL)
			}//end for(int i = 0; i < (int) snp_list.size(); i++)
		}//end if(s2 != NULL)
	}// end loop through pairs map

//	if (!par::silent)
		cout << "\n";
}

void Epistasis::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm)
{
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

    string fname = opts::_OUTPREFIX_ + "epistasis" + options.getOut() + ".txt";
    if(!overwrite)
    {
        fname += "." + getString<int>(order);
    }

	ofstream EPI(fname.c_str());
	if(!EPI)
	{
		throw MethodException("Unable to open " + fname + " for writing!\n");
	}
	EPI.precision(4);

	if (!options.doEpiQuickscan())
	{

		EPI << "CHR1" << " " << "SNP1" << " " << "CHR2" << " " << "SNP2" << " ";

		if (!options.doEpiFast())
		{
			if (opts::_BINTRAIT_)
				EPI << "OR_INT" << " ";
			else
				EPI << "BETA_SNP1 BETA_SNP2 BETA_INT" << " ";
		}

		EPI << "STAT_INT" << " " << "INTP" << " ";
		if(!options.doEpiFast())
		{
			if(!opts::_BINTRAIT_)
				EPI << "RSQ" << " " << "F" << " " << "MODELP";
		}
		EPI << "\n";
	}
	else
	{

		EPI << "CHR1" << " " << "SNP1" << " " << "CHR2" << " " << "SNP2" << " "
				<< "CHISQ" << " " << "\n";
	}

	////////////////////////////////////////////////////////////////////
	// epi1 and epi2 thresholds were given in terms of 0.01 (two-sided)
	// calculate appropriate absolute Z scores

	opts::printLog("Threshold for displaying epistatic result (-epi-alpha1) : p <= "
			+ getString<double> (options.get_epi_alpha1()) + "\n");
	opts::printLog("Threshold for counting epistatic result (-epi-alpha2) : p <= "
			+ getString<double> (options.get_epi_alpha2()) + "\n");

	if(options.getBioSnpFile() != "")
	{
		evalBioFile(EPI);
		if(EPI.is_open())
		{
			EPI.close();
		}
		return;
	}
	double epi_alpha1 = fabs(Helpers::ltqnorm(options.get_epi_alpha1() / 2));
	double epi_alpha2 = fabs(Helpers::ltqnorm(options.get_epi_alpha2() / 2));

	// Fast epistasis:  case-only or case/control
	// Regression based test: case/control or quantitative trait

	// Take a list of SNPs, or all SNPs (vector<bool> epi1)
	// Test these against either themselves, or all SNPs (vector<bool> epi2)

	//  A     B
	//  ALL x ALL    skip e1>e2
	//  SET1 x ALL
	//  SET1 x SET1  skip e1>e2
	//  SET1 x SET2

	bool skip_symm = false;

	// Only output epistatic tests that have p < par::epi_alpha1;
	// Do not even attempt to save any epistatic results -- go straight to STDOUT

	// Also present summary results for all epi1 SNPs
	// (i.e. average / proportion of significant epistatic tests
	//  at a certain alpha level, par::epi_alpha2

	int nl_all = data_set->num_loci();
	vector<bool> sA(nl_all, false);
	vector<bool> sB(nl_all, false);

	vector<int> good_indexes = Helpers::findValidMarkersIndexes(data_set->get_markers(), &options);

	// Are we using a test set? If so, construct now
	if (options.doEpiSets())
	{
		if(options.getEpiSetsFilename() != "")
		{
			options.readEpiSets(options.getEpiSetsFilename());
		}
		map<string, vector<string> > snpset = options.getEpiSets();
		if (snpset.size() > 2)
			throw MethodException("Can only specify one or two SETs when testing for epistasis\n");
		if (snpset.size() == 0)
			throw MethodException("There are no valid sets specified\n");

		map<string, vector<string> >::iterator setiter = snpset.begin();
		vector<string> first_set_str = setiter->second;
		vector<string> second_set_str;
		vector<int> first_set;
		for(int e = 0; e < (int)first_set_str.size(); e++)
		{
			int loc = data_set->get_locus_index(first_set_str[e]);
			vector<int>::iterator iter = find(good_indexes.begin(), good_indexes.end(), loc);
			if(iter != good_indexes.end())
			{
				first_set.push_back(loc);
			}
		}
		for (unsigned int e = 0; e < first_set.size(); e++)
		{
			sA[first_set[e]] = true;//snpset[0][e]] = true;
		}
		// Has a second set been specified?

		if (snpset.size() == 2)
		{
			opts::printLog("SET1 x SET2 epistasis mode\n");
			setiter++;
			second_set_str = setiter->second;
			for (unsigned int e = 0; e < second_set_str.size(); e++)
			{
				int loc = data_set->get_locus_index(second_set_str[e]);
				vector<int>::iterator iter = find(good_indexes.begin(), good_indexes.end(), loc);
				if(iter != good_indexes.end())
				{
					sB[loc] = true;//snpset[1][e]] = true;
				}
			}
		}
		else if (options.doEpiSetBySet())// Otherwise, has SET x SET flag been given?
		{
			opts::printLog("SET1 x SET1 epistasis mode\n");
			skip_symm = true;
			for (unsigned int e = 0; e < first_set.size(); e++)
				sB[first_set[e]] = true;
		}
		else // All SNPs in second set
		{
			opts::printLog("SET1 x ALL epistasis mode\n");
			for (unsigned int e = 0; e < good_indexes.size(); e++)
				sB[good_indexes[e]] = true;
		}
	}
	else
	{
		opts::printLog("ALL x ALL epistasis mode\n");
		skip_symm = true;
		for (unsigned int e = 0; e < good_indexes.size(); e++)
		{
			sA[good_indexes[e]] = true;
			sB[good_indexes[e]] = true;
		}
	}

	// Use fast aff coding

	// Count how many items in the SET1

	int epc = 0;
	for (vector<bool>::iterator e1 = sA.begin(); e1 != sA.end(); e1++)
		if (*e1)
			epc++;
	int epcc = 0;
	// Keep track of how many epistatic tests actually performed
	long int nepi = 0;

	vector<int> summary_sig(nl_all, 0);
	vector<int> summary_good(nl_all, 0);
	vector<double> best_score(nl_all, 0);
	vector<int> best_partner(nl_all);


	LinearRegression lr;
	lr.setOptions(options);
	lr.resetDataSet(data_set);

	LogisticRegression logr;
	logr.set_parameters(&options);
	logr.resetDataSet(data_set);

	InputFilter ct_filter;
	vector<string> use_covs = options.getCovars();
	vector<unsigned int> covs;
	vector<unsigned int> traits;
	if(options.doCovarsName())
	{
		ct_filter.add_covariate_list(&use_covs);
		ct_filter.add_covariate_filter(InputFilter::IncludeCovariateFilter);
	}
	for(int c = 0; c < data_set->num_covariates(); c++)
	{
		bool use = true;
		for(int f = 0; f < ct_filter.num_covariate_filters(); f++)
		{
			use = ct_filter.run_covariate_filter(f, data_set->get_covariate_name(c));
		}
		if(use)
		{
			covs.push_back(c);
		}
	}

	//////////////////////////////////////////
	// Begin iterating over pairs : SET x SET

	int ssize = data_set->num_inds();
	for (int e1 = 0; e1 < nl_all; e1++)
	{
		if (sA[e1])
		{
			cout 	<< "Peforming tests of epistasis: group "
					<< getString<int>(++epcc)
					<< " of " << getString<int>(epc)
					<< "        \r";
			cout.flush();

			Marker* mark_e1 = data_set->get_locus(e1);
			if(!mark_e1->isEnabled())
				continue;
			for (int e2 = 0; e2 < nl_all; e2++)
			{
				///////////////////////////////////////////
				// Skip this test under certain conditions

				// The SNP not in the set
				if (!sB[e2])
				{
					continue;
				}

				// We've already performed this test
				if (e1 >= e2 && skip_symm)
				{
					continue;
				}

				// Same SNP
				if (e1 == e2)
				{
					continue;
				}

				Marker* mark_e2 = data_set->get_locus(e2);
				if(!mark_e2->isEnabled())
				{
					continue;
				}
				// Skip X chromosome for now
				if(opts::_CHRX_ == mark_e1->getChrom() || opts::_CHRX_ == mark_e2->getChrom())
				{
					continue;
				}

				// SNPs too close (case-only analysis)
				if (options.doEpiCaseOnly())
					if (mark_e1->getChrom() == mark_e2->getChrom())
						if (fabs((double) (mark_e1->getBPLOC() - mark_e2->getBPLOC()))
								< options.getEpiCaseOnlyKbGap() * 1000)
						{
							continue;
						}

				//////////////////////////////////
				// Perform test of epistasis here

				if (opts::_BINTRAIT_ && options.doEpiFast())
				{

					double z; // statistic from either method

					// Odds ratio test
					// make two 2x2 tables

					int a11, a12, a21, a22;
					int u11, u12, u21, u22;
					a11 = a12 = a21 = a22 = 0;
					u11 = u12 = u21 = u22 = 0;


					for(int s = 0; s < ssize; s++)
					{
						Sample* samp = data_set->get_sample(s);
						bool b1 = samp->getAone(mark_e2->getLoc());
						bool b2 = samp->getAtwo(mark_e2->getLoc());
						bool b3 = samp->getAmissing(mark_e2->getLoc());
						bool a1 = samp->getAone(mark_e1->getLoc());
						bool a2 = samp->getAtwo(mark_e1->getLoc());
						bool a3 = samp->getAmissing(mark_e1->getLoc());



						if(!samp->isEnabled() || b3 || a3)
						{
							continue;
						}
						if (samp->getAffected())
						{

							if (!b1)
							{
								if (!b2)
								{
									if (!a1)
									{
										if (!a2)
										{
											a11 += 4;
										}
										else
										{
											a11 += 2;
											a21 += 2;
										}
									}
									else if (a2)
									{
										a21 += 4;
									}
								}
								else
								{
									if (!a1)
									{
										if (!a2)
										{
											a11 += 2;
											a12 += 2;
										}
										else
										{
											a11++;
											a21++;
											a12++;
											a22++;
										}
									}
									else if (a2)
									{
										a21 += 2;
										a22 += 2;
									}
								}
							}
							else if (b2)
							{
								if (!a1)
								{
									if (!a2)
									{
										a12 += 4;
									}
									else
									{
										a12 += 2;
										a22 += 2;
									}
								}
								else if (a2)
								{
									a22 += 4;
								}

							}
						}

						// Unaffecteds?
						else if (!options.getEpiCaseOnly())
						{
							if (!b1)
							{
								if (!b2)
								{
									if (!a1)
									{
										if (!a2)
										{
											u11 += 4;
										}
										else
										{
											u11 += 2;
											u21 += 2;
										}
									}
									else if (a2)
									{
										u21 += 4;
									}
								}
								else
								{
									if (!a1)
									{
										if (!a2)
										{
											u11 += 2;
											u12 += 2;
										}
										else
										{
											u11++;
											u21++;
											u12++;
											u22++;
										}
									}
									else if (a2)
									{
										u21 += 2;
										u22 += 2;
									}
								}
							}
							else if (b2)
							{
								if (!a1)
								{
									if (!a2)
									{
										u12 += 4;
									}
									else
									{
										u12 += 2;
										u22 += 2;
									}
								}
								else if (a2)
								{
									u22 += 4;
								}

							}

						}

					}

					// Calculate log(OR) and SEs

					double or_aff, v_aff, or_unf, v_unf;

					or_aff = log((double) (a11 * a22) / (double) (a12 * a21));
					v_aff = 1 / (double) a11 + 1 / (double) a12 + 1
							/ (double) a21 + 1 / (double) a22;

					// Case-only z-score (if requested)
					if (options.getEpiCaseOnly())
					{
						z = fabs(or_aff / sqrt(v_aff));
					}
					else // Standard case-control analysis
					{
						or_unf = log((double) (u11 * u22)
								/ (double) (u12 * u21));
						v_unf = 1 / (double) u11 + 1 / (double) u12 + 1
								/ (double) u21 + 1 / (double) u22;
						z = fabs((or_aff - or_unf) / sqrt(v_aff + v_unf));
					}

					//////////////////////////////
					// --nop option in effect
					// Just output z score, if valid & above threshold

					if (options.getEpiQuickscan())
					{
						// Is this worth recording?
						if (Helpers::realnum(z))
						{
							nepi++;

							if (z >= epi_alpha1)
							{
								EPI << mark_e1->getChrom() << " "
										<< mark_e1->getRSID() << " "
										<< mark_e2->getChrom() << " "
										<< mark_e2->getRSID()
										<< " " << z * z << "\n";
							}
							EPI.flush();
							continue;
						}
					}

					/////////////////////////////////
					// More full parsing of results

					// Check this is a proper result

					if (options.getEpiFilter() && Helpers::realnum(z))
					{
						// One more test performed
						nepi++;
						// Count as a good result
						summary_good[e1]++;
						if (sA[e2])
						{
							summary_good[e2]++;
						}

						// Do we want to record this as part of the summary for the first set?
						if (z >= epi_alpha2)
						{
							// first variable will always be in A set
							summary_sig[e1]++;
							// but the second may also be in A set
							if (sA[e2])
							{
								summary_sig[e2]++;
							}
						}

						// Is this result the best scrore yet for marker in set A?
						if (z > best_score[e1])
						{
							best_score[e1] = z;
							best_partner[e1] = e2;
						}

						// The second marker might also be in set A
						if (sA[e2])
						{
							if (z > best_score[e2])
							{
								best_score[e2] = z;
								best_partner[e2] = e1;
							}
						}

						// Is this worth recording?

						if (z >= epi_alpha1)
						{
							EPI << mark_e1->getChrom() << " "
							    << mark_e1->getRSID() << " "
									<< mark_e2->getChrom() << " "
									<< mark_e2->getRSID() << " "
									<< " " << z * z << " " << Helpers::normdist(-z) * 2 << " " << "\n";
							EPI.flush();
						}
						else
						{
							continue; // skip to next pair (skip logistic test)
						}

					}
					else if (!options.getEpiFilter())
					{
						// Record all results here, whether NA or otherwise
						EPI << mark_e1->getChrom() << " "
						    << mark_e1->getRSID() << " "
								<< mark_e2->getChrom() << " "
								<< mark_e2->getRSID() << " "
								<< z * z << " "
								<< Helpers::normdist(-z) * 2 << " " << "\n";
						EPI.flush();
					}
					else
					{
						continue; // if bad statistic for this test, do not try logistic
					}

				}
				// End of binary OR test


				///////////////////////////////////////////////
				// Logistic or linear regression epistasis test

				if (!options.doEpiFast())
				{
					vector<double> b;
					double chisq = 0;
					double F = 0;
					double fstat = 0;
					double r2 = 0;

					if (opts::_BINTRAIT_)
					{
						//do Logistic Regression because binary trait...
						vector<unsigned int> model;
						model.push_back((unsigned int)e1);
						model.push_back((unsigned int)e2);

						logr.setFullInteraction(true);
						if(covs.size() > 0)
						{
							logr.calculate(model, covs, traits);
						}
						else
						{
							logr.calculate(model);
						}

						b = logr.getCoefficients();
						vector<double> ses = logr.getCoeffStandardErr();
						if(ses.size() > 0 && b.size() > 0)
						{
							double se = ses[2];
							double Z = b[2] / se;
							chisq = Z*Z;
						}
					}
					else
					{
						lr.reset();
						vector<int> model;
						model.push_back(e1);
						model.push_back(e2);
						//have the model...need to be able to give model, ?covariates?, and ?traits? to lr.calculate()...
						lr.calculate(model);

						b = lr.getCoefs();
						chisq = lr.getStatistic();

						F = lr.findF();
						fstat = lr.getFStat();
						r2 = lr.calculateRSquared();

					}
					bool notvalid = false;
					if(!Helpers::realnum(chisq))
					{
						notvalid = true;
						chisq = 0;
					}
					double pvalue = Helpers::p_from_chi(chisq, 1);
					double z = sqrt(chisq);

					// Is this result worth displaying?

					if (!notvalid)
					{
						// One more valid test performed
						nepi++;

						// Count as a good result

						summary_good[e1]++;
						if (sA[e2])
						{
							summary_good[e2]++;
						}

						// Do we want to record this as part of the summary for the first set?
						if (z >= epi_alpha2)
						{
							// first variable will always be in A set
							summary_sig[e1]++;

							// but the second may also be in A set
							if (sA[e2])
							{
								summary_sig[e2]++;
							}
						}

						// Is this result the best scrore yet for marker in set A?

						if (z > best_score[e1])
						{
							best_score[e1] = z;
							best_partner[e1] = e2;
						}

						// The second marker might also be in set A

						if (sA[e2])
						{
							if (z > best_score[e2])
							{
								best_score[e2] = z;
								best_partner[e2] = e1;
							}
						}
					}

					// Is this result worth displaying?
					// TODO:  figure out if this needs to be changed or if my data just doesn't give good enough results....
					//if (z >= epi_alpha1)
					if (true)
					{
						EPI << mark_e1->getChrom() << " "
						    << mark_e1->getRSID() << " "
							<< mark_e2->getChrom() << " "
							<< mark_e2->getRSID() << " ";
						if (!notvalid)
						{
							if (opts::_BINTRAIT_)
								EPI << exp(b[2]) << " "
									<< chisq << " " << pvalue
									<< " " << "\n";
							else
								EPI << b[1] << " " << b[2] << " " << b[3] << " "
									<< chisq << " " << pvalue
									<< " " << r2
									<< " " << fstat
									<< " " << F
									<< " " << "\n";
							vector<string>labels = lr.getLabels();

							if(labels.size() > 1 && !opts::_BINTRAIT_)
							{
								vector<double>pvals = lr.getPvalues();
								vector<double>coefs = lr.getCoefs();
								vector<double> vars = lr.getVar();
								vector<double> zs = lr.getZs();
								for(int l = 1; l < (int)labels.size(); l++)
								{
									bool okay = vars[l] < 1e-20 || !Helpers::realnum(vars[l]) ? false : true;
									double se = 0;
									if(okay)
									{
										se = sqrt(vars[l]);
									}
									EPI << mark_e1->getChrom() << "\t" << mark_e1->getRSID() << "\t" << mark_e2->getChrom() << "\t" << mark_e2->getRSID() << "\t";
									EPI << labels[l] << "\t";
									EPI << coefs[l] << "\t";
									EPI << exp(coefs[l]) << "\t";
									EPI << se << endl;
								}
							}
						}
						else
							EPI << "NA" << " " << "NA"
								<< " " << "NA" << " " << "\n";

						EPI.flush();

					}
				}
			} // Next pair of SNPs
		}// end [if (sA[e1])

	}//end main loop through all loci

		cout << "\n";
	EPI.close();

	//////////////////////
	// Summary of results

	// Skip this for now

	if (true)
	{
		fname += ".summary";
		EPI.open(fname.c_str(), ios::out);
		EPI.clear();

		opts::printLog("Performed a total of " + getString<int>(nepi)
				+ " valid SNPxSNP tests\n");

		opts::printLog("Writing epistasis summary results to [ " + fname + " ] \n");

		EPI.precision(4);
		EPI << "CHR" << " " << "SNP" << " "
				<< "N_SIG" << " " << "N_TOT" << " "
				<< "PROP" << " " << "BEST_CHISQ" << " "
				<< "BEST_CHR" << " "
				<< "BEST_SNP" << " " << "\n";


		for (int e1 = 0; e1 < nl_all; e1++)
		{
			if (sA[e1])
			{
				EPI << data_set->get_locus(e1)->getChrom() << " "
						<< data_set->get_locus(e1)->getRSID() << " "
						<< summary_sig[e1] << " "
						<< summary_good[e1] << " "
						<< (double) summary_sig[e1] / (double) summary_good[e1]
						<< " " << best_score[e1] * best_score[e1]
						<< " " << data_set->get_locus(best_partner[e1])->getChrom()
						<< " "
						<< data_set->get_locus(best_partner[e1])->getRSID() << " " << "\n";

			}
		}

		EPI.close();
	}
}// end Epistasis::process()


}

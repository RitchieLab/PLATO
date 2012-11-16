#include <unistd.h>

#include <math.h>
#include "AlleleInfo.h"
#define PI 3.141592653589793238462643
#define DEBUG 0

void AlleleInfo::merge(allele_struct_data a){
	allele1.assign(a.allele1);
	int i = allele1.find(" ", 0);
	if(i != string::npos){
		allele1.erase(i);
	}
	allele2.assign(a.allele2);
	i = allele2.find(" ", 0);
	if(i != string::npos){
		allele2.erase(i);
	}
//	major.assign(a.major);
//	i = major.find(" ", 0);
//	if(i != string::npos){
//		major.erase(i);
//	}
	parent_major_count = a.parent_major_count;
	parent_minor_count = a.parent_minor_count;
	father_major_count = a.father_major_count;
	father_minor_count = a.father_minor_count;
	mother_major_count = a.mother_major_count;
	mother_minor_count = a.mother_minor_count;
	child_major_count = a.child_major_count;
	child_minor_count = a.child_minor_count;
	child_major_count_f = a.child_major_count_f;
	child_minor_count_f = a.child_minor_count_f;
	child_major_count_m = a.child_major_count_m;
	child_minor_count_m = a.child_minor_count_m;
	//casecontrol
	case_major_count_f = a.case_major_count_f;
	case_minor_count_f = a.case_minor_count_f;
	case_major_count_m = a.case_major_count_m;
	case_minor_count_m = a.case_minor_count_m;
	control_major_count_f = a.control_major_count_f;
	control_minor_count_f = a.control_minor_count_f;
	control_major_count_m = a.control_major_count_m;
	control_minor_count_m = a.control_minor_count_m;
	
	parent_het_count = a.parent_het_count;
	father_het_count = a.father_het_count;
	mother_het_count = a.mother_het_count;
	child_het_count = a.child_het_count;
	child_het_count_f = a.child_het_count_f;
	child_het_count_m = a.child_het_count_m;
	//casecontrol
	case_het_count_f = a.case_het_count_f;
	case_het_count_m = a.case_het_count_m;
	control_het_count_f = a.control_het_count_f;
	control_het_count_m = a.control_het_count_m;
	
	parent_homo1_count = a.parent_homo1_count;
	parent_homo2_count = a.parent_homo2_count;
	father_homo1_count = a.father_homo1_count;
	father_homo2_count = a.father_homo2_count;
	mother_homo1_count = a.mother_homo1_count;
	mother_homo2_count = a.mother_homo2_count;
	child_homo1_count = a.child_homo1_count;
	child_homo2_count = a.child_homo2_count;
	child_homo1_count_f = a.child_homo1_count_f;
	child_homo2_count_f = a.child_homo2_count_f;
	child_homo1_count_m = a.child_homo1_count_m;
	child_homo2_count_m = a.child_homo2_count_m;
	//casecontrol
	case_homo1_count_f = a.case_homo1_count_f;
	case_homo2_count_f = a.case_homo2_count_f;
	case_homo1_count_m = a.case_homo1_count_m;
	case_homo2_count_m = a.case_homo2_count_m;
	control_homo1_count_f = a.control_homo1_count_f;
	control_homo2_count_f = a.control_homo2_count_f;
	control_homo1_count_m = a.control_homo1_count_m;
	control_homo2_count_m = a.control_homo2_count_m;
	
	parents_used = a.parents_used;
	fathers_used = a.fathers_used;
	mothers_used = a.mothers_used;
	children_used = a.children_used;
	children_used_f = a.children_used_f;
	children_used_m = a.children_used_m;
	//casecontrol
	cases_used_f = a.cases_used_f;
	cases_used_m = a.cases_used_m;
	controls_used_f = a.controls_used_f;
	controls_used_m = a.controls_used_m;
	
	parent_hw = a.parent_hw;
	father_hw = a.father_hw;
	mother_hw = a.mother_hw;
	child_hw = a.child_hw;
	child_hw_f = a.child_hw_f;
	child_hw_m = a.child_hw_m;
	//casecontrol
	case_hw_f = a.case_hw_f;
	case_hw_m = a.case_hw_m;
	control_hw_f = a.control_hw_f;
	control_hw_m = a.control_hw_m;
	
/*cout << "Alleles merged!" << endl;
cout << "--";
cout << "allele1 = " << allele1 << endl;
cout << "--";
cout << "allele2 = " << allele2 << endl;
cout << "--";
cout << "pmc = " << parent_major_count << endl;
cout << "--";
cout << "pminc = " << parent_minor_count << endl;
cout << "--";
cout << "fmc = " << father_major_count << endl;
cout << "--";
cout << "fminc = " << father_minor_count << endl;
cout << "--";
cout << "mmc = " << mother_major_count << endl;
cout << "--";
cout << "mminc = " << mother_minor_count << endl;
cout << "--";
cout << "cmc = " << child_major_count << endl;
cout << "--";
cout << "cminc = " << child_minor_count << endl;
cout << "--";
cout << "phc = " << parent_het_count << endl;
cout << "--";
cout << "fhc = " << father_het_count << endl;
cout << "--";
cout << "mhc = " << mother_het_count << endl;
cout << "--";
cout << "chc = " << child_het_count << endl;
cout << "--";
cout << "ph1c = " << parent_homo1_count << endl;
cout << "--";
cout << "ph2c = " << parent_homo2_count << endl;
cout << "--";
cout << "fh1c = " << father_homo1_count << endl;
cout << "--";
cout << "fh2c = " << father_homo2_count << endl;
cout << "--";
cout << "mh1c = " << mother_homo1_count << endl;
cout << "--";
cout << "mh2c = " << mother_homo2_count << endl;
cout << "--";
cout << "ch1c = " << child_homo1_count << endl;
cout << "--";
cout << "ch2c = " << child_homo2_count << endl;
cout << "--";
cout << "puc = " << parents_used << endl;
cout << "--";
cout << "fuc = " << fathers_used << endl;
cout << "--";
cout << "muc = " << mothers_used << endl;
cout << "--";
cout << "cuc = " << children_used << endl;
cout << "--";
cout << "phw = " << parent_hw << endl;
cout << "--";
cout << "fhw = " << father_hw << endl;
cout << "--";
cout << "mhw = " << mother_hw << endl;
cout << "--";
cout << "chw = " << child_hw << endl;
*/
	calcFreqs();
}

float AlleleInfo::getFatherHomo1Exp(){
	    float p2 = (float)pow(father_major_freq, 2);
		    int population = father_homo1_count + father_homo2_count + father_het_count;
			    float exp = (float) p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getFatherHomo2Exp(){
	    float p2 = (float)pow(father_minor_freq, 2);
		    int population = father_homo1_count + father_homo2_count + father_het_count;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}

float AlleleInfo::getMotherHomo1Exp(){
	    float p2 = (float)pow(mother_major_freq, 2);
		    int population = mother_homo1_count + mother_homo2_count + mother_het_count;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getMotherHomo2Exp(){
	    float p2 = (float)pow(mother_minor_freq, 2);
		    int population = mother_homo1_count + mother_homo2_count + mother_het_count;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getParentHomo1Exp(){
	    float p2 = (float)pow(parent_major_freq, 2);
		    int population = parent_homo1_count + parent_homo2_count + parent_het_count;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}

float AlleleInfo::getParentHomo2Exp(){
	    float p2 = (float)pow(parent_minor_freq, 2);
		    int population = parent_homo1_count + parent_homo2_count + parent_het_count;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getChildHomo1Exp(){
	    float p2 = (float)pow(child_major_freq, 2);
		    int population = child_homo1_count + child_homo2_count + child_het_count;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getChildHomo2Exp(){
	    float p2 = (float)pow(child_minor_freq, 2);
		    int population = child_homo1_count + child_homo2_count + child_het_count;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getChildHomo1Exp_F(){
	    float p2 = (float)pow(child_major_freq_f, 2);
		    int population = child_homo1_count_f + child_homo2_count_f + child_het_count_f;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getChildHomo2Exp_F(){
	    float p2 = (float)pow(child_minor_freq_f, 2);
		    int population = child_homo1_count_f + child_homo2_count_f + child_het_count_f;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getChildHomo1Exp_M(){
	    float p2 = (float)pow(child_major_freq_m, 2);
		    int population = child_homo1_count_m + child_homo2_count_m + child_het_count_m;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getChildHomo2Exp_M(){
	    float p2 = (float)pow(child_minor_freq_m, 2);
		    int population = child_homo1_count_m + child_homo2_count_m + child_het_count_m;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
//casecontrol
float AlleleInfo::getCaseHomo1Exp_F(){
	    float p2 = (float)pow(case_major_freq_f, 2);
		    int population = case_homo1_count_f + case_homo2_count_f + case_het_count_f;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getCaseHomo2Exp_F(){
	    float p2 = (float)pow(case_minor_freq_f, 2);
		    int population = case_homo1_count_f + case_homo2_count_f + case_het_count_f;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getCaseHomo1Exp_M(){
	    float p2 = (float)pow(case_major_freq_m, 2);
		    int population = case_homo1_count_m + case_homo2_count_m + case_het_count_m;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getCaseHomo2Exp_M(){
	    float p2 = (float)pow(case_minor_freq_m, 2);
		    int population = case_homo1_count_m + case_homo2_count_m + case_het_count_m;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getControlHomo1Exp_F(){
	    float p2 = (float)pow(control_major_freq_f, 2);
		    int population = control_homo1_count_f + control_homo2_count_f + control_het_count_f;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getControlHomo2Exp_F(){
	    float p2 = (float)pow(control_minor_freq_f, 2);
		    int population = control_homo1_count_f + control_homo2_count_f + control_het_count_f;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getControlHomo1Exp_M(){
	    float p2 = (float)pow(control_major_freq_m, 2);
		    int population = control_homo1_count_m + control_homo2_count_m + control_het_count_m;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}
float AlleleInfo::getControlHomo2Exp_M(){
	    float p2 = (float)pow(control_minor_freq_m, 2);
		    int population = control_homo1_count_m + control_homo2_count_m + control_het_count_m;
			    float exp = (float)p2 * (float) population;
				
				    return exp;
}


float AlleleInfo::getFatherHetExp(){
	    float pq = 2 * father_major_freq * father_minor_freq;
		    int population = father_homo1_count + father_homo2_count + father_het_count;
			    float exp = (float)pq * (float) population;
				
				    return exp;
}

float AlleleInfo::getMotherHetExp(){
	    float pq = 2 * mother_major_freq * mother_minor_freq;
		    int population = mother_homo1_count + mother_homo2_count + mother_het_count;
			    float exp = (float)pq * (float) population;
				
				    return exp;
}

float AlleleInfo::getParentHetExp(){
    float pq = 2 * parent_major_freq * parent_minor_freq;
    int population = parent_homo1_count + parent_homo2_count + parent_het_count;
    float exp = (float)pq * (float) population;

    return exp;
}

float AlleleInfo::getChildHetExp(){
    float pq = 2 * child_major_freq * child_minor_freq;
    int population = child_homo1_count + child_homo2_count + child_het_count;
    float exp = (float)pq * (float) population;

    return exp;
}
float AlleleInfo::getChildHetExp_F(){
    float pq = 2 * child_major_freq_f * child_minor_freq_f;
    int population = child_homo1_count_f + child_homo2_count_f + child_het_count_f;
    float exp = (float)pq * (float) population;

    return exp;
}
float AlleleInfo::getChildHetExp_M(){
    float pq = 2 * child_major_freq_m * child_minor_freq_m;
    int population = child_homo1_count_m + child_homo2_count_m + child_het_count_m;
    float exp = (float)pq * (float) population;

    return exp;
}
//casecontrol
float AlleleInfo::getCaseHetExp_F(){
    float pq = 2 * case_major_freq_f * case_minor_freq_f;
    int population = case_homo1_count_f + case_homo2_count_f + case_het_count_f;
    float exp = (float)pq * (float) population;

    return exp;
}
float AlleleInfo::getCaseHetExp_M(){
    float pq = 2 * case_major_freq_m * case_minor_freq_m;
    int population = case_homo1_count_m + case_homo2_count_m + case_het_count_m;
    float exp = (float)pq * (float) population;

    return exp;
}
float AlleleInfo::getControlHetExp_F(){
    float pq = 2 * control_major_freq_f * control_minor_freq_f;
    int population = control_homo1_count_f + control_homo2_count_f + control_het_count_f;
    float exp = (float)pq * (float) population;

    return exp;
}
float AlleleInfo::getControlHetExp_M(){
    float pq = 2 * control_major_freq_m * control_minor_freq_m;
    int population = control_homo1_count_m + control_homo2_count_m + control_het_count_m;
    float exp = (float)pq * (float) population;

    return exp;
}


float AlleleInfo::calcHW_exact(int obs_hets, int obs_hom1, int obs_hom2){
/*
 // This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
	// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
	// // Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000  
	// //
	// // Written by Jan Wigginton
	// */
	if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) 
	{
		printf("FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a"
		" negative count", obs_hets, obs_hom1, obs_hom2);
		exit(EXIT_FAILURE);
	}
	
	int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;
	
	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes   = obs_hets + obs_homc + obs_homr;
	if(genotypes <= 0){
		return -1;
	}
	
	double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	if (het_probs == NULL) 
	{
		printf("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities" );
		exit(EXIT_FAILURE);
	}
	                         
	int i;
	for (i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;
	   /* start at midpoint */
	    int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
	    /* check to ensure that midpoint and rare alleles have same parity */
	    if ((rare_copies & 1) ^ (mid & 1))
		    mid++;
	
		int curr_hets = mid;
	    int curr_homr = (rare_copies - mid) / 2;
	    int curr_homc = genotypes - curr_hets - curr_homr;
	
	    het_probs[mid] = 1.0;
	    double sum = het_probs[mid];
	    for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
	    {
		    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
	        sum += het_probs[curr_hets - 2];
	
	       /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
	        curr_homr++;
	        curr_homc++;
	    }
	
	    curr_hets = mid;
	    curr_homr = (rare_copies - mid) / 2;
	    curr_homc = genotypes - curr_hets - curr_homr;
	    for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
	    {
		    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc /((curr_hets + 2.0) * (curr_hets + 1.0));
	        sum += het_probs[curr_hets + 2];
	
	       /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
	        curr_homr--;
	        curr_homc--;
	    }
	
	    for (i = 0; i <= rare_copies; i++)
		    het_probs[i] /= sum;
	
	 	/* alternate p-value calculation for p_hi/p_lo
	    double p_hi = het_probs[obs_hets];
	    for (i = obs_hets + 1; i <= rare_copies; i++)
		    p_hi += het_probs[i];
	
		double p_lo = het_probs[obs_hets];
	    for (i = obs_hets - 1; i >= 0; i--)
		    p_lo += het_probs[i];
	
	    double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
	    */
	
	    double p_hwe = 0.0;
	    /*  p-value calculation for p_hwe  */
	    for (i = 0; i <= rare_copies; i++)
	    {
		    if (het_probs[i] > het_probs[obs_hets])
	 	       continue;
	        p_hwe += het_probs[i];
	    }
	    p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
	    free(het_probs);
	 return (float) p_hwe;
}

float AlleleInfo::calcHW(float a1freq, float a2freq, int obsa1, int obsa2, int obshet, int population){
	if(DEBUG){
		cout << "Initial:\t" << a1freq << "\t" << a2freq << "\t" << obsa1 << "\t" << obsa2 << "\t" << obshet << "\t" << population << endl;
	}
	float p2 = pow(a1freq, 2);
	float q2 = pow(a2freq,2);
	float pq = 2 * a1freq * a2freq;
	if(DEBUG){
		cout << "P2,Q2,PQ:\t" << p2 << "\t" << q2 << "\t" << pq << endl;
	}
	float expp2 = p2 * (float)population;
	float expq2 = q2 * (float)population;
	float exppq = pq * (float)population;
	if(expq2 < 5){
		return calcHW_exact(obshet, obsa1, obsa2);
	}
	if(DEBUG){
		cout << "Exp:\t" << expp2 << "\t" << expq2 << "\t" << exppq << endl;
	}
	float chi1 = 0.0;
	float chi2 = 0.0;
	float chi3 = 0.0;

	if(expp2 > 0){
		chi1 = (pow(((float)obsa1 - expp2), 2)) / expp2;
	}
	if(expq2 > 0){
		chi2 = (pow(((float)obsa2 - expq2),2)) / expq2;
	}
	if(exppq > 0){
		chi3 = (pow(((float)obshet - exppq),2)) / exppq;
	}

	
	float chi = chi1 + chi2 + chi3;

	if(DEBUG){
		cout << "ChiTot:\t" << chi1 << "\t" << chi2 << "\t" << chi3 << "\t=" << chi << endl;
	}
	//float results = _subchisqrprob(1,chi);
	float results = ChiSquare::pfromchi(chi, 1);
	if(DEBUG){
		cout << "Results: " << results << endl;
	}
	return results;
}

float AlleleInfo::_subchisqrprob(int deg, float chi){
	float pval = 0.0;
	if(chi <= 0){
		pval = 1;
	}
	else if(deg > 100){
		pval = _subuprob((pow((chi / (float)deg),((float)1/(float)3)) - (float)(1.0 - 2.0/9.0/(float)deg)) / sqrt(2.0/9.0/(float)deg));
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

float AlleleInfo::_subuprob(float val){
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
			p = (float)i / (absx + p);
		}
		p = exp(-.5 * absx * absx) / sqrt(2 * PI) / (absx + p);
	}
	if(val < 0){
		p = 1 - p;
	}
	return p;
}

void AlleleInfo::process(bool a1, bool a2, string type){
	bool found = false;
    if(type == "C"){
        if(!a1){
            child_major_count++;
            found = true;
        }
	    if(!a2){
	        child_major_count++;
	        found = true;
	    }
	    if(a1){
	        child_minor_count++;
	        found = true;
	    }
	    if(a2){
	        child_minor_count++;
	        found = true;
	    }
	    if(found){
	        if(!a1 && !a2){
		        child_homo1_count++;
		    }
		    if(!a1 && a2){
		        child_het_count++;
		    }
            if(a1 && a2){
                child_homo2_count++;
            }
            children_used++;
        }
    }
	else if(type == "CF"){
        if(!a1){
            child_major_count_f++;
            found = true;
        }
	    if(!a2){
	        child_major_count_f++;
	        found = true;
	    }
	    if(a1){
	        child_minor_count_f++;
	        found = true;
	    }
	    if(a2){
	        child_minor_count_f++;
	        found = true;
	    }
	    if(found){
	        if(!a1 && !a2){
		        child_homo1_count_f++;
		    }
		    if(!a1 && a2){
		        child_het_count_f++;
		    }
            if(a1 && a2){
                child_homo2_count_f++;
            }
            children_used_f++;
        }
    }
	else if(type == "CM"){
        if(!a1){
            child_major_count_m++;
            found = true;
        }
	    if(!a2){
	        child_major_count_m++;
	        found = true;
	    }
	    if(a1){
	        child_minor_count_m++;
	        found = true;
	    }
	    if(a2){
	        child_minor_count_m++;
	        found = true;
	    }
	    if(found){
	        if(!a1 && !a2){
		        child_homo1_count_m++;
		    }
		    if(!a1 && a2){
		        child_het_count_m++;
		    }
            if(a1 && a2){
                child_homo2_count_m++;
            }
            children_used_m++;
        }
    }
    else if(type == "P"){
		if(!a1){
            parent_major_count++;
            found = true;
        }
        if(!a2){
            parent_major_count++;
            found = true;
        }
        if(a1){
            parent_minor_count++;
            found = true;
        }
        if(a2){
            parent_minor_count++;
            found = true;
        }
        if(found){
	    if(!a1 && !a2){
                parent_homo1_count++;
            }
            if(!a1 && a2){
                parent_het_count++;
            }
            if(a1 && a2){
                parent_homo2_count++;
            }
            parents_used++;
        }
    }
    else if(type == "F"){
        if(!a1){
            father_major_count++;
            found = true;
        }
        if(!a2){
            father_major_count++;
            found = true;
        }
        if(a1){
            father_minor_count++;
            found = true;
        }
        if(a2){
            father_minor_count++;
            found = true;
        }
        if(found){
            if(!a1 && !a2){
                father_homo1_count++;
            }
            if(!a1 && a2){
                father_het_count++;
            }
            if(a1 && a2){
                father_homo2_count++;
            }
            fathers_used++;
        }
    }
    else if(type == "M"){
        if(!a1){
            mother_major_count++;
            found = true;
        }
        if(!a2){
            mother_major_count++;
            found = true;
        }
        if(a1){
            mother_minor_count++;
            found = true;
        }
        if(a2){
            mother_minor_count++;
            found = true;
        }
        if(found){
            if(!a1 && !a2){
                mother_homo1_count++;
            }
            if(!a1 && a2){
                mother_het_count++;
            }
            if(a1 && a2){
                mother_homo2_count++;
            }
            mothers_used++;
        }
    }
    else if(type == "CASEF"){
        if(!a1){
            case_major_count_f++;
            found = true;
        }
        if(!a2){
            case_major_count_f++;
            found = true;
        }
        if(a1){
            case_minor_count_f++;
            found = true;
        }
        if(a2){
            case_minor_count_f++;
            found = true;
        }
        if(found){
            if(!a1 && !a2){
                case_homo1_count_f++;
            }
            if(!a1 && a2){
                case_het_count_f++;
            }
            if(a1 && a2){
                case_homo2_count_f++;
            }
            cases_used_f++;
        }
    }
    else if(type == "CASEM"){
        if(!a1){
            case_major_count_m++;
            found = true;
        }
        if(!a2){
            case_major_count_m++;
            found = true;
        }
        if(a1){
            case_minor_count_m++;
            found = true;
        }
        if(a2){
            case_minor_count_m++;
            found = true;
        }
        if(found){
            if(!a1 && !a2){
                case_homo1_count_m++;
            }
            if(!a1 && a2){
                case_het_count_m++;
            }
            if(a1 && a2){
                case_homo2_count_m++;
            }
            cases_used_m++;
        }
    }
    else if(type == "CONTROLM"){
        if(!a1){
            control_major_count_m++;
            found = true;
        }
        if(!a2){
            control_major_count_m++;
            found = true;
        }
        if(a1){
            control_minor_count_m++;
            found = true;
        }
        if(a2){
            control_minor_count_m++;
            found = true;
        }
        if(found){
            if(!a1 && !a2){
                control_homo1_count_m++;
            }
            if(!a1 && a2){
                control_het_count_m++;
            }
            if(a1 && a2){
                control_homo2_count_m++;
            }
            controls_used_m++;
        }
    }
    else if(type == "CONTROLF"){
        if(!a1){
            control_major_count_f++;
            found = true;
        }
        if(!a2){
            control_major_count_f++;
            found = true;
        }
        if(a1){
            control_minor_count_f++;
            found = true;
        }
        if(a2){
            control_minor_count_f++;
            found = true;
        }
        if(found){
            if(!a1 && !a2){
                control_homo1_count_f++;
            }
            if(!a1 && a2){
                control_het_count_f++;
            }
            if(a1 && a2){
                control_homo2_count_f++;
            }
            controls_used_f++;
        }
    }
}

void AlleleInfo::process(string a1, string a2, string type){
	bool found = false;
    if(type == "C"){
        if(a1 == allele1){
            child_major_count++;
            found = true;
        }
	    if(a2 == allele1){
	        child_major_count++;
	        found = true;
	    }
	    if(a1 == allele2){
	        child_minor_count++;
	        found = true;
	    }
	    if(a2 == allele2){
	        child_minor_count++;
	        found = true;
	    }
	    if(found){
	        if(a1 == allele1 && a2 == allele1){
		        child_homo1_count++;
		    }
		    if(a1 == allele1 && a2 == allele2 && a1 != a2){
		        child_het_count++;
		    }
            if(a1 == allele2 && a2 == allele2){
                child_homo2_count++;
            }
            children_used++;
        }
    }
	else if(type == "CF"){
        if(a1 == allele1){
            child_major_count_f++;
            found = true;
        }
	    if(a2 == allele1){
	        child_major_count_f++;
	        found = true;
	    }
	    if(a1 == allele2){
	        child_minor_count_f++;
	        found = true;
	    }
	    if(a2 == allele2){
	        child_minor_count_f++;
	        found = true;
	    }
	    if(found){
	        if(a1 == allele1 && a2 == allele1){
		        child_homo1_count_f++;
		    }
		    if(a1 == allele1 && a2 == allele2 && a1 != a2){
		        child_het_count_f++;
		    }
            if(a1 == allele2 && a2 == allele2){
                child_homo2_count_f++;
            }
            children_used_f++;
        }
    }
	else if(type == "CM"){
        if(a1 == allele1){
            child_major_count_m++;
            found = true;
        }
	    if(a2 == allele1){
	        child_major_count_m++;
	        found = true;
	    }
	    if(a1 == allele2){
	        child_minor_count_m++;
	        found = true;
	    }
	    if(a2 == allele2){
	        child_minor_count_m++;
	        found = true;
	    }
	    if(found){
	        if(a1 == allele1 && a2 == allele1){
		        child_homo1_count_m++;
		    }
		    if(a1 == allele1 && a2 == allele2 && a1 != a2){
		        child_het_count_m++;
		    }
            if(a1 == allele2 && a2 == allele2){
                child_homo2_count_m++;
            }
            children_used_m++;
        }
    }
    else if(type == "P"){
		if(a1 == allele1){
            parent_major_count++;
            found = true;
        }
        if(a2 == allele1){
            parent_major_count++;
            found = true;
        }
        if(a1 == allele2){
            parent_minor_count++;
            found = true;
        }
        if(a2 == allele2){
            parent_minor_count++;
            found = true;
        }
        if(found){
	    if(a1 == allele1 && a2 == allele1){
                parent_homo1_count++;
            }
            if(a1 == allele1 && a2 == allele2 && a1 != a2){
                parent_het_count++;
            }
            if(a1 == allele2 && a2 == allele2){
                parent_homo2_count++;
            }
            parents_used++;
        }
    }
    else if(type == "F"){
        if(a1 == allele1){
            father_major_count++;
            found = true;
        }
        if(a2 == allele1){
            father_major_count++;
            found = true;
        }
        if(a1 == allele2){
            father_minor_count++;
            found = true;
        }
        if(a2 == allele2){
            father_minor_count++;
            found = true;
        }
        if(found){
            if(a1 == allele1 && a2 == allele1){
                father_homo1_count++;
            }
            if(a1 == allele1 && a2 == allele2 && a1 != a2){
                father_het_count++;
            }
            if(a1 == allele2 && a2 == allele2){
                father_homo2_count++;
            }
            fathers_used++;
        }
    }
    else if(type == "M"){
        if(a1 == allele1){
            mother_major_count++;
            found = true;
        }
        if(a2 == allele1){
            mother_major_count++;
            found = true;
        }
        if(a1 == allele2){
            mother_minor_count++;
            found = true;
        }
        if(a2 == allele2){
            mother_minor_count++;
            found = true;
        }
        if(found){
            if(a1 == allele1 && a2 == allele1){
                mother_homo1_count++;
            }
            if(a1 == allele1 && a2 == allele2 && a1 != a2){
                mother_het_count++;
            }
            if(a1 == allele2 && a2 == allele2){
                mother_homo2_count++;
            }
            mothers_used++;
        }
    }
    else if(type == "CASEF"){
        if(a1 == allele1){
            case_major_count_f++;
            found = true;
        }
        if(a2 == allele1){
            case_major_count_f++;
            found = true;
        }
        if(a1 == allele2){
            case_minor_count_f++;
            found = true;
        }
        if(a2 == allele2){
            case_minor_count_f++;
            found = true;
        }
        if(found){
            if(a1 == allele1 && a2 == allele1){
                case_homo1_count_f++;
            }
            if(a1 == allele1 && a2 == allele2 && a1 != a2){
                case_het_count_f++;
            }
            if(a1 == allele2 && a2 == allele2){
                case_homo2_count_f++;
            }
            cases_used_f++;
        }
    }
    else if(type == "CASEM"){
        if(a1 == allele1){
            case_major_count_m++;
            found = true;
        }
        if(a2 == allele1){
            case_major_count_m++;
            found = true;
        }
        if(a1 == allele2){
            case_minor_count_m++;
            found = true;
        }
        if(a2 == allele2){
            case_minor_count_m++;
            found = true;
        }
        if(found){
            if(a1 == allele1 && a2 == allele1){
                case_homo1_count_m++;
            }
            if(a1 == allele1 && a2 == allele2 && a1 != a2){
                case_het_count_m++;
            }
            if(a1 == allele2 && a2 == allele2){
                case_homo2_count_m++;
            }
            cases_used_m++;
        }
    }
    else if(type == "CONTROLM"){
        if(a1 == allele1){
            control_major_count_m++;
            found = true;
        }
        if(a2 == allele1){
            control_major_count_m++;
            found = true;
        }
        if(a1 == allele2){
            control_minor_count_m++;
            found = true;
        }
        if(a2 == allele2){
            control_minor_count_m++;
            found = true;
        }
        if(found){
            if(a1 == allele1 && a2 == allele1){
                control_homo1_count_m++;
            }
            if(a1 == allele1 && a2 == allele2 && a1 != a2){
                control_het_count_m++;
            }
            if(a1 == allele2 && a2 == allele2){
                control_homo2_count_m++;
            }
            controls_used_m++;
        }
    }
    else if(type == "CONTROLF"){
        if(a1 == allele1){
            control_major_count_f++;
            found = true;
        }
        if(a2 == allele1){
            control_major_count_f++;
            found = true;
        }
        if(a1 == allele2){
            control_minor_count_f++;
            found = true;
        }
        if(a2 == allele2){
            control_minor_count_f++;
            found = true;
        }
        if(found){
            if(a1 == allele1 && a2 == allele1){
                control_homo1_count_f++;
            }
            if(a1 == allele1 && a2 == allele2 && a1 != a2){
                control_het_count_f++;
            }
            if(a1 == allele2 && a2 == allele2){
                control_homo2_count_f++;
            }
            controls_used_f++;
        }
    }
 //   calcFreqs();
}

void AlleleInfo::calcFreqs(){
	int parent_total = parent_major_count + parent_minor_count;
	int child_total = child_major_count + child_minor_count;
	int child_total_f = child_major_count_f + child_minor_count_f;
	int child_total_m = child_major_count_m + child_minor_count_m;
	//casecontrol
	int case_total_f = case_major_count_f + case_minor_count_f;
	int case_total_m = case_major_count_m + case_minor_count_m;
	int control_total_f = control_major_count_f + control_minor_count_f;
	int control_total_m = control_major_count_m + control_minor_count_m;
	
    int father_total = father_major_count + father_minor_count;
	int mother_total = mother_major_count + mother_minor_count;
	if(parent_major_count < parent_minor_count){
		int temp = parent_major_count;
		parent_major_count = parent_minor_count;
		parent_minor_count = temp;
	}
	if(parent_homo1_count < parent_homo2_count){
		int temp = parent_homo1_count;
		parent_homo1_count = parent_homo2_count;
		parent_homo2_count = temp;
	}
	if(parent_total > 0){
		parent_major_freq = (float)((float)parent_major_count / (float)parent_total);
		parent_minor_freq = (float)((float)parent_minor_count / (float)parent_total);
		if(parent_major_freq < parent_minor_freq){
			float temp = parent_major_freq;
			parent_major_freq = parent_minor_freq;
			parent_minor_freq = temp;
		}
	}
	else{
		parent_major_freq = 0.0;
		parent_minor_freq = 0.0;
	}

	if(child_major_count < child_minor_count){
		int temp = child_major_count;
		child_major_count = child_minor_count;
		child_minor_count = temp;
	}
	if(child_homo1_count < child_homo2_count){
		int temp = child_homo1_count;
		child_homo1_count = child_homo2_count;
		child_homo2_count = temp;
	}
	if(child_total > 0){
		child_major_freq = (float)child_major_count / (float)child_total;
		child_minor_freq = (float)child_minor_count / (float)child_total;
		if(child_major_freq < child_minor_freq){
			float temp = child_major_freq;
			child_major_freq = child_minor_freq;
			child_minor_freq = temp;
		}
	}
	else{
		child_major_freq = 0.0;
		child_minor_freq = 0.0;
	}
	
	if(child_major_count_f < child_minor_count_f){
		int temp = child_major_count_f;
		child_major_count_f = child_minor_count_f;
		child_minor_count_f = temp;
	}
	if(child_homo1_count_f < child_homo2_count_f){
		int temp = child_homo1_count_f;
		child_homo1_count_f = child_homo2_count_f;
		child_homo2_count_f = temp;
	}
	if(child_total_f > 0){
		child_major_freq_f = (float)child_major_count_f / (float)child_total_f;
		child_minor_freq_f = (float)child_minor_count_f / (float)child_total_f;
		if(child_major_freq_f < child_minor_freq_f){
			float temp = child_major_freq_f;
			child_major_freq_f = child_minor_freq_f;
			child_minor_freq_f = temp;
		}
	}
	else{
		child_major_freq_f = 0.0;
		child_minor_freq_f = 0.0;
	}
	
	if(child_major_count_m < child_minor_count_m){
		int temp = child_major_count_m;
		child_major_count_m = child_minor_count_m;
		child_minor_count_m = temp;
	}
	if(child_homo1_count_m < child_homo2_count_m){
		int temp = child_homo1_count_m;
		child_homo1_count_m = child_homo2_count_m;
		child_homo2_count_m = temp;
	}
	if(child_total_m > 0){
		child_major_freq_m = (float)child_major_count_m / (float)child_total_m;
		child_minor_freq_m = (float)child_minor_count_m / (float)child_total_m;
		if(child_major_freq_m < child_minor_freq_m){
			float temp = child_major_freq_m;
			child_major_freq_m = child_minor_freq_m;
			child_minor_freq_m = temp;
		}
	}
	else{
		child_major_freq_m = 0.0;
		child_minor_freq_m = 0.0;
	}

	//casecontrol
	if(case_major_count_f < case_minor_count_f){
		int temp = case_major_count_f;
		case_major_count_f = case_minor_count_f;
		case_minor_count_f = temp;
	}
	if(case_homo1_count_f < case_homo2_count_f){
		int temp = case_homo1_count_f;
		case_homo1_count_f = case_homo2_count_f;
		case_homo2_count_f = temp;
	}
	if(case_total_f > 0){
		case_major_freq_f = (float)case_major_count_f / (float)case_total_f;
		case_minor_freq_f = (float)case_minor_count_f / (float)case_total_f;
		if(case_major_freq_f < case_minor_freq_f){
			float temp = case_major_freq_f;
			case_major_freq_f = case_minor_freq_f;
			case_minor_freq_f = temp;
		}
	}
	else{
		case_major_freq_f = 0.0;
		case_minor_freq_f = 0.0;
	}
	
	if(case_major_count_m < case_minor_count_m){
		int temp = case_major_count_m;
		case_major_count_m = case_minor_count_m;
		case_minor_count_m = temp;
	}
	if(case_homo1_count_m < case_homo2_count_m){
		int temp = case_homo1_count_m;
		case_homo1_count_m = case_homo2_count_m;
		case_homo2_count_m = temp;
	}
	if(case_total_m > 0){
		case_major_freq_m = (float)case_major_count_m / (float)case_total_m;
		case_minor_freq_m = (float)case_minor_count_m / (float)case_total_m;
		if(case_major_freq_m < case_minor_freq_m){
			float temp = case_major_freq_m;
			case_major_freq_m = case_minor_freq_m;
			case_minor_freq_m = temp;
		}
	}
	else{
		case_major_freq_m = 0.0;
		case_minor_freq_m = 0.0;
	}
	if(control_major_count_f < control_minor_count_f){
		int temp = control_major_count_f;
		control_major_count_f = control_minor_count_f;
		control_minor_count_f = temp;
	}
	if(control_homo1_count_f < control_homo2_count_f){
		int temp = control_homo1_count_f;
		control_homo1_count_f = control_homo2_count_f;
		control_homo2_count_f = temp;
	}
	if(control_total_f > 0){
		control_major_freq_f = (float)control_major_count_f / (float)control_total_f;
		control_minor_freq_f = (float)control_minor_count_f / (float)control_total_f;
		if(control_major_freq_f < control_minor_freq_f){
			float temp = control_major_freq_f;
			control_major_freq_f = control_minor_freq_f;
			control_minor_freq_f = temp;
		}
	}
	else{
		control_major_freq_f = 0.0;
		control_minor_freq_f = 0.0;
	}
	
	if(control_major_count_m < control_minor_count_m){
		int temp = control_major_count_m;
		control_major_count_m = control_minor_count_m;
		control_minor_count_m = temp;
	}
	if(control_homo1_count_m < control_homo2_count_m){
		int temp = control_homo1_count_m;
		control_homo1_count_m = control_homo2_count_m;
		control_homo2_count_m = temp;
	}
	if(control_total_m > 0){
		control_major_freq_m = (float)control_major_count_m / (float)control_total_m;
		control_minor_freq_m = (float)control_minor_count_m / (float)control_total_m;
		if(control_major_freq_m < control_minor_freq_m){
			float temp = control_major_freq_m;
			control_major_freq_m = control_minor_freq_m;
			control_minor_freq_m = temp;
		}
	}
	else{
		control_major_freq_m = 0.0;
		control_minor_freq_m = 0.0;
	}
	
	
	if(father_major_count < father_minor_count){
		int temp = father_major_count;
		father_major_count = father_minor_count;
		father_minor_count = temp;
	}
	if(father_homo1_count < father_homo2_count){
		int temp = father_homo1_count;
		father_homo1_count = father_homo2_count;
		father_homo2_count = temp;
	}
    if(father_total > 0){
        father_major_freq = (float)father_major_count / (float)father_total;
        father_minor_freq = (float)father_minor_count / (float)father_total;
		if(father_major_freq < father_minor_freq){
			float temp = father_major_freq;
			father_major_freq = father_minor_freq;
			father_minor_freq = temp;
		}
    }
    else{
        father_major_freq = 0.0;
        father_minor_freq = 0.0;
    }

	if(mother_major_count < mother_minor_count){
		int temp = mother_major_count;
		mother_major_count = mother_minor_count;
		mother_minor_count = temp;
	}
	if(mother_homo1_count < mother_homo2_count){
		int temp = mother_homo1_count;
		mother_homo1_count = mother_homo2_count;
		mother_homo2_count = temp;
	}
    if(mother_total > 0){
        mother_major_freq = (float)mother_major_count / (float)mother_total;
        mother_minor_freq = (float)mother_minor_count / (float)mother_total;
		if(mother_major_freq < mother_minor_freq){
			float temp = mother_major_freq;
			mother_major_freq = mother_minor_freq;
			mother_minor_freq = temp;
		}
    }
    else{
        mother_major_freq = 0.0;
        mother_minor_freq = 0.0;
    }
}

float AlleleInfo::getABSParentMinorFreq(){
	calcFreqs();
	//cout << "Parent minorfreq: " << parent_minor_freq << " " << parent_major_freq << endl;
	if(parent_minor_freq > parent_major_freq){
		return parent_major_freq;
	}
	return parent_minor_freq;
}

float AlleleInfo::getABSFatherMinorFreq(){
	calcFreqs();
    if(father_minor_freq > father_major_freq){
        return father_major_freq;
    }
    return father_minor_freq;
}

float AlleleInfo::getABSMotherMinorFreq(){
	calcFreqs();
    if(mother_minor_freq > mother_major_freq){
        return mother_major_freq;
    }
    return mother_minor_freq;
}


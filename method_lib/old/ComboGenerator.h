//ComboGenerator.h

#ifndef __COMBOGENERATOR_H__
#define __COMBOGENERATOR_H__

#include <vector>
namespace Methods{
///
/// Generates combinations of numbers <br>
/// For example, will generate all 2 number combinations from
/// 1 to 100: <br> <i> 1-2, 1-3, 1-4, 1-5, 1-6, 1-7, etc.</i> <br>
/// Adapted from Will Bush's sMDR code
///

/// Generates numerical combinations
class ComboGenerator{
  
  public:
    ComboGenerator();
    bool GenerateCombinations();  
    bool GenerateCombinations(int new_interval);
    
    bool AdvanceParameters();
    
    void SetComboInterval(int combInterval);
    int GetComboInterval(){return ComboInterval;}
    
    void ComboEnds(int combStart, int combEnd);
    void SetLoci(int nLoci);
    std::vector < std::vector<unsigned int> > ComboList;

    int get_combo_max(){return ComboEnd;}
    int get_combo_min(){return ComboStart;}
    
    int get_count_generated(){return counter+1;}
    
    void initialize_state();
    
    int param_j(){return j;}
    void param_j(int new_j){j=new_j;}
    
    int param_x(){return x;}
    void param_x(int new_x){x=new_x;}
    
    int param_kdec(){return kdec;}
    void param_kdec(int new_kdec){kdec=new_kdec;}
    
    bool param_AlreadyStarted(){return AlreadyStarted;}
    void param_AlreadyStarted(bool new_AlreadyStarted){AlreadyStarted=new_AlreadyStarted;}
    
    int param_c(int index){return c[index];}
    void set_param_c(int index, int value){c[index] = value;}
    
    double calc_combinations(double n, unsigned int k);
  
  private:
    void initialize(); 
    int ComboInterval, ComboStart, ComboEnd, NumLoci;
    int j, x, kdec, counter;
    int * c;
    bool AlreadyStarted;
    
};
};

#endif

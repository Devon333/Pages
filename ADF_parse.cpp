#include "ADF_File.h"
//#include "TAPE_READ.h"
#include <iostream> 
using namespace std;
//      if(orbitalSymmetry[i] == orbitalSymmetry[i+1])
//      {
//          AllPercent += abs(stod(percents[i]));
//          cout <<"if 1 " <<  this->orbitalSymmetry[i] << "  "<< this->orbitalSymmetry[i+1] << "  "<<this->orbitalType[i]<< "  " << percentSumm << endl; 
//          cout << orbitalSymmetry[i] << "  " << percents[i] << endl;
//      }
//
//      else if(orbitalSymmetry[i] != orbitalSymmetry[i+1])
//      {
//          AllPercent += abs(stod(percents[i]));
//          cout <<orbitalSymmetry[i] << "  " << percents[i] << "  All percents " << AllPercent << endl;
//          this->AllOrbSum.push_back(AllPercent);
//          allSto=AllPercent;
//          cout <<"if 1 " <<  this->orbitalSymmetry[i] << "  "<< this->orbitalSymmetry[i+1] << "  "<<this->orbitalType[i]<< "  " << percentSumm << endl; 
//          cout <<"new sao "<< percentSto / allSto << endl;
//          AllPercent = 0;  
//      }

  






/* Using a timer example
    clock_t start, end;
    start = clock();
    ios_base::sync_with_stdio(true);

        end = clock();
         double time_taken = double(end-start)/double(CLOCKS_PER_SEC);

         cout << "Time taken to calculate coupings between "<< numMOs << " molecular orbitals " << fixed << time_taken<< setprecision(16);
         cout << " sec" << endl;
//END Using a timer example 
*/

int main(int argc , char** argv)
{
  string filename = argv[1];//ADF TDDFT calculation output file
  string param = argv[2];//TAPE FILE included or not
  if(param == "TAPE")
  {
      string filename2 = argv[3];//TAPE file name
      string TAPESYMM = argv[4];//Symmetry of excited states in TAPE File
      int i_state = stoi(argv[5]);//initial state number
      int f_state = stoi(argv[6]);//final state number
//      READ_TAPE test(filename, filename2);
      vector<string> supOrb;
      for(int i=7; i< argc; i++)
      {
          supOrb.push_back(argv[i]);
      }
      cout <<argc << "  " << supOrb.size() << "  " << supOrb[0] << endl;
      ADF_File nap(filename, filename2, supOrb, TAPESYMM, i_state, f_state);
  }  
  else
  {
      vector<string> supOrb;
      for(int i=2; i< argc; i++)
      {
          supOrb.push_back(argv[i]);
      }
      cout <<argc << "  " << supOrb.size() << "  " << supOrb[0] << endl;
      ADF_File nap(filename, supOrb);
  }



}//END MAIN BRACKET


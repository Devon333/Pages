#include<iostream>
#include<fstream>
#include<string> 
#include<sstream>
#include<vector>
#include<algorithm>
#include<bits/stdc++.h>//insert function
#include"TAPE_READ.h"

using namespace std;
///put transition data into vectors
class ADF_File:public READ_TAPE{

  private: 

  public:
//  struct MO_Transitions
//  {
//      vector<string> dipoleAllowedIrreps;
//      vector<string> excNum;
//      vector<string> occOrb;
//      vector<string> virOrb;
//      vector<string> weight;
//      vector<string> dipx;
//      vector<string> dipy;
//      vector<string> dipz;
//  };

 

  


  void ADF_READ(string filename);
  void Find_Irreps(vector<string> FileLines);
  //void Find_Irrep_Exc(vector<string> FileLines);
  void Find_Irrep_Exc(vector<string> FileLines, string TAPESYMM);
  void Find_MOTransition(vector<string> FileLines, string TAPESYMM);
  void Find_MOTransition(vector<string> FileLines);
  void showTransitions();
  void Find_MOS(vector<string> FileLines);
  void show_MOs(); //shows orbital symmetry, percentage orbital contribution, orbital type
  void calcSuperAtomic(vector<string> supOrb); // calculates the sum of superatomic contributions to an orbital
  void show_SuperAtomic(); //shows the orbital symmetry and percentage of superatomic character
  void match_Orbitals(string filename); //matches orbital symmetry in excited state with superatomic character
  void newSet(string filename);
  void newSet(string filename, string TAPESYMM, vector<double> TAPEdipadd);
  void calcCollect();
  void calcDipAdd();
  void showDipAdd();
  void check_MOs();
  void readTAPEadditivities(int i_state, int f_state);


  string dotColFind(string lineElement, string SorP);
  string replaceSubstring(string element, string old , string newstring );
  vector<string> split_Line(string Line);  
  

  //class variables/containers
  vector<string> FileLines;
  vector<string> Irreps;
  vector<vector<string> > IrrepExc;
  vector<string> dipoleAllowedIrreps;
  int IrrepWithExc;
  vector<string> orbitalSymmetry;
  vector<string> percents;
  vector<string> orbitalType;

  //holds transition information
  vector<int> positionCounts;//shows positions in vector where excitations of a different symmetry begin
  vector<string> pointGroup;
  vector<string> excNum;
  vector<string> occOrb;
  vector<string> virOrb;
  vector<string> weight;
  vector<string> dipx;
  vector<string> dipy;
  vector<string> dipz;
  vector<double> dipadd;

  //holds energy and oscillator information about excitations
  vector<int> energyCounts;//shows positions in vector where excitations of a different symmetry begin
  vector<string> energies;
  vector<string> energyNumbers;
  vector<string> oscStr;
  vector<double> collectivity;
  vector<double> SAOCoeffSum;
  vector<double> TAPEADDI;
  //holds orbital symmetry and percent of superatomic orbital
  vector<string> SAOSymmetry;//###
  vector<double> SAOSum;//########
  vector<double> AllOrbSum;//##### 

  //class constructors
  ADF_File()
  {
    cout << "constructing ADF File " << endl;
  }

  ADF_File(string filename, vector<string> supOrb){
      //MO_Transitions Trans[5];
      ADF_READ(filename);
      Find_Irreps(FileLines);
      cout << "made it here " << endl;
      //Find_Irrep_Exc(FileLines);
      cout << "made it here " << endl;
      //MO_Transitions Trans[IrrepWithExc+1];
      Find_MOTransition(FileLines);
      showTransitions();
      Find_MOS(FileLines);
      //show_MOs();
      //vector<string> orbs=[1P:z];
      //vector<string> supOrb={"1P:z","2P:z"};
      check_MOs();
      calcSuperAtomic(supOrb);
      show_SuperAtomic();
      calcDipAdd();
      showDipAdd();
      calcCollect();
      newSet(filename);
      //match_Orbitals(filename);
  }

  ADF_File(string filename, string filename2, vector<string> supOrb,string TAPESYMM, int i_state, int f_state){

      //READ_TAPE(filename, filename2, this->TAPEADDI, i_state, f_state);
      ADF_READ(filename);
      Find_Irreps(FileLines);
      cout << "made it here " << endl;
      Find_Irrep_Exc(FileLines,TAPESYMM);
      cout << "made it here " << endl;
      Find_MOTransition(FileLines,TAPESYMM);
      showTransitions();
      Find_MOS(FileLines);
      //show_MOs();
      //vector<string> orbs=[1P:z];
      //vector<string> supOrb={"1P:z","2P:z"};
      calcSuperAtomic(supOrb);
      show_SuperAtomic();
      check_MOs();
      //exit(0);
      //calcDipAdd();
      //showDipAdd();
      calcCollect();
      //readTAPEadditivities(i_state, f_state);
      //cout << "additivity size " << TAPEADDI.size() << endl;
      newSet(filename,TAPESYMM,TAPEADDI);
      //match_Orbitals(filename);
  }


};

void ADF_File::readTAPEadditivities(int i_state, int f_state)
{
  string GrabLine;
  string runner;
  ifstream TAPEout;
  TAPEout.open("Tapeadditivities"+to_string(i_state)+"_"+to_string(f_state)+".out", ios::in);

  if(!TAPEout){cout << "Unable to open " << "Tapeadditivities.out" << " !! \n"; exit(1) ;} /*check if the file is open*/
  cout << "Tapeadditivities"+to_string(i_state)+"_"+to_string(f_state)+".out" << " is open" << endl;
  int counter = 0;
  while(getline(TAPEout,runner))
  {  
      GrabLine = runner;
      TAPEADDI.push_back(stod(GrabLine));
      counter++;
  }
  TAPEout.close();
}

vector<string> ADF_File::split_Line(string Line)
{
    istringstream LineSplit(Line);
    string LineElement;
    vector<string> splitLine; 
    while(LineSplit >> LineElement)
    {
        splitLine.push_back(LineElement); 
    }
    return splitLine;
}



void ADF_File::newSet(string Filename)
{
  ofstream outFile;
  string name = Filename.erase(Filename.size()-4,Filename.size());
  outFile.open(name+".dat");
  int j=0; int betweenIrreps=0;
  double SAOSUMM=0; int orbcount=0;
  int state_count=0;

  outFile<< "vecInt \t ExcNum \t Energy \t OscStr \t dipAdd \t SAOChar \t collectivity \t dipAllowIrrep \t IrrepTag" << endl;
  for(int i=0;i<IrrepWithExc;i++)
  {
      betweenIrreps+=positionCounts[i]; 
      outFile<< i << " " << dipoleAllowedIrreps[i] << endl;
      outFile <<"position Counts " << positionCounts[i] <<  "  sizes  irreps with Exc " << IrrepWithExc << "  exc Num " << excNum.size() << " SAO SYmm " << SAOSymmetry.size() << "SAO summ " << SAOSum.size()<<  "  Occ orb " << occOrb.size()  << " dipadd " << dipadd.size() << endl;
//      outFile << betweenIrreps << endl;
      for(j;j<betweenIrreps;j++)
      {
      //cout << "occOrb[j] " << occOrb[j] << " weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << endl;
      string stateNumber=excNum[j];
//          if(stateNumber != excNum[j+1])
//          {
//              for(int k=0; k<SAOSymmetry.size() ; k++)
//              {
//                  //cout << "neq occOrb[j] " << occOrb[j] << " SAO sum " << SAOSUMM << endl;
//                  if(SAOSymmetry[k] == occOrb[j])
//                  {
//                     cout << "I'm here occOrb[j] " << occOrb[j] << " SAO sum " << SAOSUMM << endl;
//                     // outFile << "made it here" << endl;
//                      SAOSUMM += abs(stod(weight[j]));//sqrt(stod(weight[j]));//SAOSum[k]*stod(weight[j]);
//                      orbcount++;
//                  }
//              }
//          
//          }

          while(stateNumber == excNum[j+1])
          {
              cout << " in new set occOrb[j] " << occOrb[j] << endl;
              for(int k=0; k< SAOSymmetry.size() ; k++)
              {
                  cout << "all occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << endl;
                  if(SAOSymmetry[k] == occOrb[j] )//&& SAOSum[k]>0.5)
                  {
                      cout << "eq occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << endl;
                      SAOSUMM += abs(stod(weight[j]));//SAOSum[k];//stod(weight[j]);
                      orbcount++;
                  }
              }
              j++;
          }
              for(int k=0; k< SAOSymmetry.size() ; k++)
              {
                  cout << "all occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << endl;
                  if(SAOSymmetry[k] == occOrb[j] )//&& SAOSum[k]>0.5)
                  {
                      cout << "eq occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << endl;
                      SAOSUMM += abs(stod(weight[j]));//SAOSum[k];//stod(weight[j]);
                      orbcount++;
                  }
              }
          cout << "sum after set " << SAOSUMM << endl;

          if(SAOSUMM > 0.000){
          outFile << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
               <<"  "<< dipadd[state_count] << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i << endl;   
          SAOSUMM=0; orbcount=0; state_count++;
          }
          else{
          outFile << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
               <<"  "<< dipadd[state_count] << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i << endl;   
          SAOSUMM=0; orbcount=0; state_count++;
          }
       
      }
      outFile << "\n\n";
  }
  outFile.close();
}

void ADF_File::newSet(string Filename, string TAPESYMM, vector<double> TAPEdipadd)
{
  cout << "using newSet() " << endl;
  ofstream outFile;
  string name = Filename.erase(Filename.size()-4,Filename.size());
  outFile.open(name+TAPESYMM+".dat");
  int j=0; int betweenIrreps=0;
  double SAOSUMM=0; int orbcount=0;
  int state_count=0;
  cout << "SUPER ATOMIC ORBITAL SYMMETRIES " << SAOSymmetry.size() << endl;
  for(int symm=0;symm<SAOSymmetry.size();symm++)
  {
    cout << SAOSymmetry[symm] << endl;
  }  
  cout << "SUPER ATOMIC CHARACTER SUM " << SAOSum.size() << endl;
  for(int cuts=0;cuts<SAOSum.size();cuts++)
  {
    cout << SAOSum[cuts] << endl;
  }

  outFile<< "vecInt \t ExcNum \t Energy \t OscStr \t dipAdd \t SAOChar \t collectivity \t dipAllowIrrep \t IrrepTag" << endl;
  cout << "additivities size " << TAPEADDI.size() << endl;
  cout << "exc Num size " << excNum.size() << endl;
  cout << "collectivity size " << collectivity.size() << endl;
  cout << "IrrepWithExc num " << IrrepWithExc << endl;
  cout << "occOrb size " << occOrb.size() << endl;
  cout << "weight size " << weight.size() << endl;
  cout << "dipoleAllowedIrreps size " << dipoleAllowedIrreps.size() << endl;
  int i=0;
  for(i ;i<IrrepWithExc;i++)
  {
      if(TAPESYMM == dipoleAllowedIrreps[i])
      {
          
          if(i !=0)
          {
            j += positionCounts[i];
              for(int k=0;k<=i;k++)
              {
                  betweenIrreps+=positionCounts[i];
              }
          }
          else if(i==0)
          { 
              j=0;
              betweenIrreps +=positionCounts[0];
          }
          cout << "position count " << j << endl;
          cout << "betweenIrreps " << betweenIrreps << endl;
      }
 
          outFile<< i << " " << dipoleAllowedIrreps[i] << endl;
          outFile <<"position Counts " << positionCounts[i] <<  "  sizes  irreps with Exc " << IrrepWithExc << "  exc Num " << excNum.size() << " SAO SYmm " << SAOSymmetry.size() << "SAO summ " << SAOSum.size()<<  "  Occ orb " << occOrb.size()  << " dipadd " << dipadd.size() << endl;
  
    //      outFile << betweenIrreps << endl;
          for(j;j<=betweenIrreps-1;j++)
          {
          cout << "occOrb[j] " << occOrb[j] <<" occOrb[j+1] "<< occOrb[j+1]<< "  weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << "  " << j << "  " << state_count << endl;
          string stateNumber=excNum[j];
    //          if(stateNumber != excNum[j+1])
    //          {
    //              for(int k=0; k<SAOSymmetry.size() ; k++)
    //              {
    //                  //cout << "neq occOrb[j] " << occOrb[j] << " SAO sum " << SAOSUMM << endl;
    //                  if(SAOSymmetry[k] == occOrb[j])
    //                  {
    //                     cout << "I'm here occOrb[j] " << occOrb[j] << " SAO sum " << SAOSUMM << endl;
    //                     // outFile << "made it here" << endl;
    //                      SAOSUMM += abs(stod(weight[j]));//sqrt(stod(weight[j]));//SAOSum[k]*stod(weight[j]);
    //                      orbcount++;
    //                  }
    //              }
    //          
    //          }
           if(stoi(excNum[j]) < 300 && j+1 < excNum.size() && j < weight.size())
           {   
              while(stateNumber == excNum[j+1])
              {
                  //if(state_count < 10){break;}
                  cout << " in new set occOrb[j] " << occOrb[j] << endl;
                  for(int k=0; k< SAOSymmetry.size() ; k++)
                  {
                      cout << "all occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << endl; //<< "  "<<excNum[j+1] << endl;
                      string s1=SAOSymmetry[k]; string s2=occOrb[j];
                      int orbCompare = s1.compare(s2);
                      //if(SAOSymmetry[k] == occOrb[j] && SAOSum[k]>0.5)
                      if(orbCompare == 0 && SAOSum[k]>0.0)
                      {
                          cout << "eq occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << endl;// << "  "<<excNum[j+1] << endl;
                          SAOSUMM += abs(stod(weight[j]));//SAOSum[k];//stod(weight[j]);
                          orbcount++;
                      }
                  }

                  j++;
              }
              for(int k=0; k< SAOSymmetry.size() ; k++)
              {
                  cout << "all occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << endl;
                  string s1=SAOSymmetry[k]; string s2=occOrb[j];
                  int orbCompare = s1.compare(s2);
                  if(orbCompare == 0 && SAOSum[k]>=0.0)
                  {
                      cout << "final contribution eq occOrb[j] " << occOrb[j] <<" SAOSymmetry " << SAOSymmetry[k]<< " weight " << weight[j]<< " exc num "<< excNum[j] << "  "<<excNum[j+1] << endl;
                      SAOSUMM += abs(stod(weight[j]));//SAOSum[k];//stod(weight[j]);
                      orbcount++;
                  }
                  else{continue;}
              }
              cout << "Final sum before printing " << SAOSUMM << endl;
              //if(SAOSUMM > 0.000){
              //testing cout << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
              //testing     <<"  "<< /*TAPEADDI[state_count]*/0.0 << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i <<  "  " << j << endl;   
              //testing outFile << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
              //testing     <<"  "<< /*TAPEADDI[state_count]*/0.0 << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i << endl;   
              cout << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
                   <<"  "<< /*TAPEADDI[state_count]*/0.0 << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i <<  "  " << j << endl;   
              outFile << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
                   <<"  "<< /*TAPEADDI[state_count]*/ 0.0 << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i << endl;   
              SAOSUMM=0; orbcount=0; state_count++;

              //}
              if(stod(energies[state_count]) > 15.00)
              {
                  break;
              }
              cout << "sum after set " << SAOSUMM << endl;
            
//              else{
//              //testing cout << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
//              //testing     <<"  "<< /*TAPEADDI[state_count]*/0.0 << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i <<  "  " << j << endl;   
//              //testing outFile << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
//              //testing     <<"  "<< /*TAPEADDI[state_count]*/0.0 << "  " << SAOSUMM <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i << endl;   
//              cout << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
//                   <<"  "<< TAPEADDI[state_count] << "  " << 0.00 <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i << "  " << j << endl;   
//              outFile << state_count << "  " << excNum[j] <<"  "<<stod(energies[state_count]) <<"  "<< stod(oscStr[state_count])
//                   <<"  "<< TAPEADDI[state_count] << "  " << 0.00 <<"  "<< collectivity[state_count] << "  " << dipoleAllowedIrreps[i] << "  " << i << endl;   
//              SAOSUMM=0; orbcount=0; state_count++;
//              if(stod(energies[state_count]) > 15.00)
//              {
//                  break;
//              }
//
//              }
           
          
          }
          }
          outFile << "\n\n";
  }
outFile.close();
}





void ADF_File::match_Orbitals(string Filename)
{
  cout << "sizes " << SAOSum.size() << "  " << SAOSymmetry.size() << endl;
  string name = Filename.erase(Filename.size()-4,Filename.size());
  ofstream outFile;
  outFile.open(name+".dat");
  int j=0; int betweenIrreps=0; int state_count=0;
  outFile << "excNum  energy  osc.str.  dipoleSum  SAOSum  Coll.  totalweight  Symmetry  SymmetryID" << "\n";
  for(int i=0;i<IrrepWithExc;i++)
  {
      double sumx=0;double sumy=0;double sumz=0;
      double SAOSumm=0; double coll=0; double total=0;
      double sumdipx=0; double absumdipx=0;
      double sumdipy=0; double absumdipy=0;
      double sumdipz=0; double absumdipz=0;
      betweenIrreps+=positionCounts[i];
//      cout << dipoleAllowedIrreps[i] <<endl << endl;
      for(j;j<betweenIrreps;j++)
      {
          string statenumber = excNum[j];
          total += stod(weight[j]);
          coll += stod(weight[j]) * stod(weight[j]);
          sumdipx += (stod(weight[j])) * stod(dipx[j]) ; absumdipx +=(stod(weight[j])) * abs(stod(dipx[j]));
          sumdipy += (stod(weight[j])) * stod(dipy[j]) ; absumdipy +=(stod(weight[j])) * abs(stod(dipy[j]));
          sumdipz += (stod(weight[j])) * stod(dipz[j]) ; absumdipz +=(stod(weight[j])) * abs(stod(dipz[j]));
         
         //cout << j << "  " << excNum[j] <<"  "<< occOrb[j] 
         //      <<"  "<< virOrb[j] <<"  "<<sqrt(stod(weight[j])) * stod(dipx[j])
         //      <<"  "<<sqrt(stod(weight[j])) * stod(dipy[j])<<"  "<<sqrt(stod(weight[j])) * stod(dipz[j]) 
         //      << endl;   
         //cout << sumdipx << "  " <<sumdipy << "  " << sumdipz<<endl;
         // j++;
          int orbcount=0;
          while(statenumber==excNum[j+1])
          {
              j++;
             // cout << j << "  " << excNum[j] <<"  "<< occOrb[j] 
             //       <<"  "<< virOrb[j] <<"  "<<sqrt(stod(weight[j])) * stod(dipx[j]) 
             //       <<"  "<<sqrt(stod(weight[j])) * stod(dipy[j])<<"  "<<sqrt(stod(weight[j])) * stod(dipz[j]) 
             //       << endl;   
              total += stod(weight[j]);
              coll += stod(weight[j])* stod(weight[j]);
              sumdipx += (stod(weight[j])) * stod(dipx[j]) ; absumdipx +=(stod(weight[j])) * abs(stod(dipx[j]));
              sumdipy += (stod(weight[j])) * stod(dipy[j]) ; absumdipy +=(stod(weight[j])) * abs(stod(dipy[j]));
              sumdipz += (stod(weight[j])) * stod(dipz[j]) ; absumdipz +=(stod(weight[j])) * abs(stod(dipz[j]));

             // cout << sumdipx << "  " <<sumdipy << "  " << sumdipz<<endl;
//              if(pointGroup[pointGroup.size()-1] ==" Symmetry: D(7D)")
//              {
//                   if(occOrb[j].find("s+.u") != std::string::npos)
//                   {
//                       SAOSumm += stod(weight[j]);
//                   }
//
//              }
//              else if(pointGroup[pointGroup.size()-1] !=" Symmetry: D(7D)")
//              {
                  for(int k=0; k<SAOSymmetry.size();k++)
                  {
                      if(occOrb[j] == SAOSymmetry[k]) //&& SAOSum[k] > 0.50)
                      {
                          cout << excNum[j] << "  " << occOrb[j] << "  " << SAOSymmetry[k] <<endl;
                          //cout << "Look here  " << SAOSum[state_count] << endl;
                          //if(SAOSum[k] > 0.50)
                          //{
                          //cout << "Look here  " <<k << "   "<< SAOSum[k] << "  " << weight[j] << endl;
                          SAOSumm += SAOSum[k];
                          orbcount++;
                          //cout <<"part " << SAOSumm << endl;
                          //cout << SAOSymmetry[k] << "  ";
                          //}
                         // else
                         // { 
                         //     SAOSumm+=0.000;             
                         // }
                      }
                      
                  }
//              }
          }
          if(SAOSumm != 0)
          {
          SAOSumm *= 1/orbcount;
          cout << "SAO SUMM " << SAOSumm << endl;
          orbcount =0;
          }
          //cout <<"total " << SAOSumm << endl;
          outFile <<std::fixed << std::setprecision(6)<< excNum[j] <<"\t"<< stod(energies[state_count]) << "\t"<< stod(oscStr[state_count])<<"\t"; 
//          cout <<"excitation number " << excNum[j] << 
//                 " energy " << energies[state_count] << 
//                 " osc.str. " << oscStr[state_count] << 
//                 " dipole x:"; if(absumdipx==0){sumx = 0.0; cout<< sumx;}else{sumx =(sumdipx/absumdipx)*100; cout << sumx;} 
//                 cout << " y: " ; if(absumdipy==0){sumy = 0.0; cout << 0.00;}else{sumy = (sumdipy/absumdipy)*100; cout << sumy;} 
//                 cout << " z: " ; if(absumdipz==0){sumz = 0.0; cout << 0.00;}else{sumz = (sumdipz/absumdipz)*100; cout << sumz;} 
//                 cout << " dipole sum " <<abs( (sumx) + (sumy) + (sumz) );
//                 cout << " SAO SUM " << SAOSumm<<endl;
                 if(absumdipx==0){sumx = 0.0;coll = 1/coll;}else{sumx =(sumdipx/absumdipx)*100; coll = 1/coll;} 
                 if(absumdipy==0){sumy = 0.0;}else{sumy = (sumdipy/absumdipy)*100;} 
                 if(absumdipz==0){sumz = 0.0;}else{sumz = (sumdipz/absumdipz)*100;} 
                 //cout << "final sum percent " << abs( (sumx) + (sumy) + (sumz) ) << endl;
                 //outFile << abs( (sumx) + (sumy) + (sumz) ); 
                 outFile << this->dipadd[state_count];
                 outFile << "\t"<< SAOSumm*100 << "\t" << coll << "\t" << total << "\t" << dipoleAllowedIrreps[i] <<  "  " << i << endl;
          total = 0;
          coll = 0;
          SAOSumm = 0;
          sumdipx = 0;
          sumdipy = 0;
          sumdipz = 0;
          sumx=0; sumy=0; sumz=0;
          state_count++;
      }
//      cout << "\n\n";
  }
  outFile.close();
}


void ADF_File::ADF_READ(string name)
{
  vector<string> lines;
  ifstream InFile;
  string runner, GrabLine; 
  InFile.open(name, ios::in);

  if(!InFile){cout << "Unable to open " << name << " !! \n"; exit(1) ;} /*check if the file is open*/
  cout << name << " is open" << endl;

  while(getline(InFile,runner))
  {
    GrabLine = runner; 
    lines.push_back(GrabLine);
    
  }
  FileLines = lines;
 
}

void ADF_File::Find_Irreps(vector<string> name)
{
  cout << "looking for Irreps" << endl;
  vector<int> IrrepLineNumbers;
  string IrrepID = "Irreducible Representations, including subspecies";
  for(int i=0;i < FileLines.size(); i++)
  {
      vector<string> temp;
      string line_ele;

      if(FileLines[i].find(IrrepID) != string::npos)
      {
          //cout << "Found Irrep " << endl;
          pointGroup.push_back(FileLines[i-2]);
          //cout << FileLines[i-2] << endl; 
          IrrepLineNumbers.push_back(i);
      }
  } 

  //for(int i=0; i<IrrepLineNumbers.size();i++){cout << IrrepLineNumbers[i] << endl;}
  //cout << IrrepLineNumbers[IrrepLineNumbers.size()-1];

      vector<string> Irrep;
      int FinalIrrepLine = IrrepLineNumbers[IrrepLineNumbers.size()-1] + 2;  
      string IrrepsDone = "Configuration of Valence Electrons";
      while(FileLines[FinalIrrepLine].find(IrrepsDone) == string::npos)
      {
          istringstream line(FileLines[FinalIrrepLine]);
          string temp;
          while(line >> temp)
          {
              Irrep.push_back(temp);
          }
          FinalIrrepLine++; 
      }

     // for(int j=0; j<Irrep.size(); j++)
     // {
     //     cout << Irrep[j] << endl;
     // }
  this->Irreps = Irrep;
  //return Irreps;
}


void ADF_File::Find_Irrep_Exc(vector<string> FileLines, string TAPESYMM)
{
  cout << "Function looking for Irreps Excited States \nNumber of lines in file: " << FileLines.size()<< "  Number of Irreps: " << Irreps.size() << endl;
  vector<string> IrrepExcitations;
  vector<vector<string> > ExcHolder(200000,vector<string>(Irreps.size()));//MAY RUNINTO ISSUE With the number of rows
  string ExcDoneID="Transition dipole moments mu (x,y,z) in a.u.";
  int IrrepCount=0;
  if(pointGroup.size() == 0){Irreps = {"A"};}
  for(int i=0;i<pointGroup.size();i++)
  {
      cout << "Point groups " << pointGroup[i] << endl;
      Irreps.push_back(TAPESYMM);

      if(pointGroup[i] ==" Symmetry: D(7D)")
      {
          cout << "IRREP WILL BE " << TAPESYMM << endl;
          Irreps.push_back(TAPESYMM);
          //dipoleAllowedIrreps = {TAPESYMM};
          //Irreps = {"S+.u","Pi.u"};
          //vector<vector<string> > ExcHolder(200000,vector<string>(Irreps.size()));//MAY RUNINTO ISSUE With the number of rows
      }
  }

  if(pointGroup[pointGroup.size()-1] == " Symmetry: C(S)")
  {
      Irreps = {"A'", "A''"};
      //vector<vector<string> > ExcHolder(200000,vector<string>(Irreps.size()));//MAY RUNINTO ISSUE With the number of rows
  }
 
  
   
  for(int irr=0; irr< 1;irr++)//Irreps.size(); irr++) //WORKING HERE
  {
      string IrrepStringID="Symmetry "+ TAPESYMM + " ";//Irreps[irr]+" ";
      cout << "These are the Irreps " << Irreps[irr] << endl;
      ExcHolder[0][IrrepCount] = Irreps[irr];
      cout << ExcHolder[0][IrrepCount] << endl;
      for(int i=0; i<FileLines.size() ; i++)
      {
          if(FileLines[i].find(IrrepStringID) != string::npos )
          {
              //cout << FileLines[i] << endl;
              i = i+7;
              int eleCounter=1;
              while(FileLines[i].find(ExcDoneID) == string::npos)
              {
                  if(FileLines[i] ==" "){break;}
                  //cout << FileLines[i] << endl;
                  istringstream splitLine(FileLines[i]);
                  string line_ele;
                  vector<string> line_holder;
                  while(splitLine >> line_ele)
                  {
                      line_holder.push_back(line_ele);
                  }
                  ExcHolder[eleCounter][IrrepCount] = line_holder[0];//no.
//                  cout << ExcHolder[eleCounter][IrrepCount] << "  ";
                  eleCounter++;
                  ExcHolder[eleCounter][IrrepCount] = line_holder[2];//E(eV)
//                  cout << ExcHolder[eleCounter][IrrepCount]<< "  ";
                  eleCounter++;
                  ExcHolder[eleCounter][IrrepCount] = line_holder[3];//f
//                  cout << ExcHolder[eleCounter][IrrepCount]<< "  ";
                  eleCounter++;
                  ExcHolder[eleCounter][IrrepCount] = line_holder[4];//dE(a.u.)
//                  cout << ExcHolder[eleCounter][IrrepCount] << endl;
                  eleCounter++;
                  //cout << eleCounter << "  " << irr << endl;
                  i++;
              }
              IrrepCount++;
              break;
          }
      }
  }
  cout << IrrepCount;
  
  this-> IrrepWithExc = IrrepCount;
  this-> IrrepExc = ExcHolder; 
}

// 40 to 50% cut off for superatomic excitaions
//
void ADF_File::Find_MOTransition(vector<string> FileLines)
{
  //MO_Transitions Trans[IrrepWithExc+1];
  cout <<"Looking for MO Transitions" << endl;
  cout <<"number of irreps found " << IrrepWithExc+1 << endl;
  string TransitionID ="Major MO -> MO transitions for the above excitations";
  string EnergyID ="Excitation energies E in a.u. and eV, dE wrt prev. cycle,";
  cout << FileLines.size()<< endl;
  cout << Irreps.size() << endl;
 // vector<string> excNum1;
 // vector<string> occOrb1;
 // vector<string> virOrb1;
 // vector<string> weight1;
 // vector<string> dipx1;
 // vector<string> dipy1;
 // vector<string> dipz1;
  //vector<vector<string> > excNum1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > occOrb1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > virOrb1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > weight1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > dipx1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > dipy1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > dipz1(10000,vector<string>(IrrepWithExc+1));


  int transition_count=0;
  int energy_counter=0;
  int searchDone = 0;
  for(int i=0; i<Irreps.size();i++)//WORKING HERE
  {
      int posCount=0;
      string IrrepID = "Symmetry "+Irreps[i]; 
      for(int l=0; l<FileLines.size();l++)
      {
//          cout << FileLines[l] << endl;
          if(FileLines[l].find(IrrepID+" ") != string::npos)
          {
              cout << "Irrep ID Found " << FileLines[l] << endl;
              this->dipoleAllowedIrreps.push_back(Irreps[i]);
              l +=7;
              while(FileLines[l] != " ")
              {
                  cout << FileLines[l] << endl;

                  vector<string> line = split_Line(FileLines[l]); //split file line by white space and puts it into a vector 

                  energyNumbers.push_back(line[0]);//excited state number
                  energies.push_back(line[2]);//energy in eV
                  oscStr.push_back(line[3]);//oscillator strength
                  l++;
                  energy_counter++;
              }
              energyCounts.push_back(energy_counter);
              cout <<"energy counter " << energy_counter<< endl;
              energy_counter=0;

              while(FileLines[l].find(TransitionID) == string::npos)
              {
                  l++;
              }
              //cout << FileLines[l] << endl;
                  if(FileLines[l].find(TransitionID) != string::npos)
                  {
                      posCount = 0;
                      cout << "found Transition ID" << endl << "Transition set found " << transition_count<< endl;
                      transition_count++;
                      l = l+9;
                  }
                      //cout << FileLines[l] << endl;
                      while(l <FileLines.size())
                      {
                          if(FileLines[l].find("Eigenvalues of small (approximate) problem") != string::npos){break;}
                          if(FileLines[l].find(" All SINGLET-SINGLET excitation energies") != string::npos){searchDone=1 ; break;}
                          if(FileLines[l].find("Transition dipole moments mu (x,y,z) in a.u.") != string::npos){searchDone=1;break;}
                          else//if(FileLines[l] != "" )//|| FileLines[l] !="\n" )
                            {
                                //cout << FileLines[l]<< endl;
                                istringstream splitstream(FileLines[l]);
                                string lineEle; 
                                vector<string> splitLine=split_Line(FileLines[l]);
                               // while(splitstream >> lineEle)
                               // {
                               //     splitLine.push_back(lineEle);
                                  //  cout << lineEle << "  " ;
                               // }
                              //  cout << endl;
                         //       cout << splitLine[0] << endl;
                                //if(splitLine.size()<8 ){break;}
                                if(splitLine.empty() != true)
                                {
                                   splitLine[0].pop_back();
                                   string excNum1=splitLine[0],occOrb1=splitLine[1], virOrb1=splitLine[3];
                                   cout << "occupied orbitals " << occOrb1 << endl;
                                    //cout << FileLines[l] << endl;
                                   this->excNum.push_back(excNum1);
                                   this->occOrb.push_back(occOrb1);
                                   this->virOrb.push_back(virOrb1);
                                   this->weight.push_back(splitLine[4]);
				   this->dipx.push_back(splitLine[5]);
				   this->dipy.push_back(splitLine[6]);
                                   this->dipz.push_back(splitLine[7]);
                                   cout << "look " << excNum1 << " "<<IrrepID << endl;
                                   //excNum1[eleCount][transition_count]=(excNum);
                                   //occOrb1[eleCount][transition_count]=(occOrb);
                                   //virOrb1[eleCount][transition_count]=(virOrb);
                                   //weight1[eleCount][transition_count]=(splitLine[4]);
				   //dipx1[eleCount][transition_count]=(splitLine[5]);
				   //dipy1[eleCount][transition_count]=(splitLine[6]);
                                   //dipz1[eleCount][transition_count]=(splitLine[7]);
                                   posCount++;
                                   cout << posCount << " position Count" << endl; 
                                }
                            }
                          //cout << eleCount << endl;
                          //eleCount++;
                          l++;
                      }
               if(searchDone == 1){break;}
          }
      }
      cout << "Position Counter " << posCount << endl;
      if(posCount != 0){
      positionCounts.push_back(posCount);}
  }

 // for(int j=0;j<transition_count;j++)
 // {
 //     for(int i=0;i<10000;i++)
 //     {
 //         if(!excNum1[i][j].empty())
 //         {
 //             this->Trans[j].excNum.push_back(excNum1[i][j]);
 //             this->Trans[j].occOrb.push_back(occOrb1[i][j]);
 //             this->Trans[j].virOrb.push_back(virOrb1[i][j]);
 //             this->Trans[j].weight.push_back(weight1[i][j]);
 //             this->Trans[j].dipx.push_back(dipx1[i][j]);
 //             this->Trans[j].dipy.push_back(dipy1[i][j]);
 //             this->Trans[j].dipz.push_back(dipz1[i][j]);
 //         } 
 //     }
 // }
}

void ADF_File::Find_MOTransition(vector<string> FileLines, string TAPESYMM)
{
  //MO_Transitions Trans[IrrepWithExc+1];
  cout <<"Looking for MO Transitions" << endl;
  cout <<"number of irreps found " << IrrepWithExc+1 << endl;
  string TransitionID ="Major MO -> MO transitions for the above excitations";
  string EnergyID ="Excitation energies E in a.u. and eV, dE wrt prev. cycle,";
  cout << FileLines.size()<< endl;
  cout << Irreps.size() << endl;
 // vector<string> excNum1;
 // vector<string> occOrb1;
 // vector<string> virOrb1;
 // vector<string> weight1;
 // vector<string> dipx1;
 // vector<string> dipy1;
 // vector<string> dipz1;
  //vector<vector<string> > excNum1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > occOrb1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > virOrb1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > weight1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > dipx1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > dipy1(10000,vector<string>(IrrepWithExc+1));
  //vector<vector<string> > dipz1(10000,vector<string>(IrrepWithExc+1));


  int transition_count=0;
  int energy_counter=0;
  int searchDone = 0;
  for(int i=0; i<1;i++)//Irreps.size();i++)//WORKING HERE
  {
      int posCount=0;
      string IrrepID = "Symmetry "+TAPESYMM+" ";//Irreps[i]; 
      
      for(int l=0; l<FileLines.size();l++)
      {
//          cout << FileLines[l] << endl;
          if(FileLines[l].find(IrrepID) != string::npos)
          {
              cout << "Irrep ID Found " << FileLines[l] << endl;
              this->dipoleAllowedIrreps.push_back(TAPESYMM);
              l +=7;
              while(FileLines[l] != " ")
              {
                  cout << FileLines[l] << endl;

                  vector<string> line = split_Line(FileLines[l]); //split file line by white space and puts it into a vector 

                  energyNumbers.push_back(line[0]);//excited state number
                  energies.push_back(line[2]);//energy in eV
                  oscStr.push_back(line[3]);//oscillator strength
                  l++;
                  energy_counter++;
              }
              energyCounts.push_back(energy_counter);
              cout <<"energy counter " << energy_counter<< endl;
              energy_counter=0;

              while(FileLines[l].find(TransitionID) == string::npos)
              {
                  l++;
              }
              //cout << FileLines[l] << endl;
                  if(FileLines[l].find(TransitionID) != string::npos)
                  {
                      posCount = 0;
                      cout << "found Transition ID" << endl << "Transition set found " << transition_count<< endl;
                      transition_count++;
                      l = l+9;
                  }
                      //cout << FileLines[l] << endl;
                      while(l <FileLines.size())
                      {
                          if(FileLines[l].find("Eigenvalues of small (approximate) problem") != string::npos){break;}
                          if(FileLines[l].find(" All SINGLET-SINGLET excitation energies") != string::npos){searchDone=1 ; break;}
                          if(FileLines[l].find("Transition dipole moments mu (x,y,z) in a.u.") != string::npos){searchDone=1;break;}
                          else//if(FileLines[l] != "" )//|| FileLines[l] !="\n" )
                            {
                                //cout << FileLines[l]<< endl;
                                istringstream splitstream(FileLines[l]);
                                string lineEle; 
                                vector<string> splitLine=split_Line(FileLines[l]);
                               // while(splitstream >> lineEle)
                               // {
                               //     splitLine.push_back(lineEle);
                                  //  cout << lineEle << "  " ;
                               // }
                              //  cout << endl;
                         //       cout << splitLine[0] << endl;
                                //if(splitLine.size()<8 ){break;}
                                if(splitLine.empty() != true)
                                {
                                   splitLine[0].pop_back();
                                   string excNum1=splitLine[0],occOrb1=splitLine[1], virOrb1=splitLine[3];
                                   cout << "occupied orbitals " << occOrb1 << endl;
                                    //cout << FileLines[l] << endl;
                                   this->excNum.push_back(excNum1);
                                   this->occOrb.push_back(occOrb1);
                                   this->virOrb.push_back(virOrb1);
                                   this->weight.push_back(splitLine[4]);
				   this->dipx.push_back(splitLine[5]);
				   this->dipy.push_back(splitLine[6]);
                                   this->dipz.push_back(splitLine[7]);
                                   cout << "look " << excNum1 << " "<<IrrepID << endl;
                                   //excNum1[eleCount][transition_count]=(excNum);
                                   //occOrb1[eleCount][transition_count]=(occOrb);
                                   //virOrb1[eleCount][transition_count]=(virOrb);
                                   //weight1[eleCount][transition_count]=(splitLine[4]);
				   //dipx1[eleCount][transition_count]=(splitLine[5]);
				   //dipy1[eleCount][transition_count]=(splitLine[6]);
                                   //dipz1[eleCount][transition_count]=(splitLine[7]);
                                   posCount++;
                                   cout << posCount << " position Count" << endl; 
                                }
                            }
                          //cout << eleCount << endl;
                          //eleCount++;
                          l++;
                      }
               if(searchDone == 1){break;}
          }
      }
      cout << "Position Counter " << posCount << endl;
      if(posCount != 0){
      positionCounts.push_back(posCount);}
  }

 // for(int j=0;j<transition_count;j++)
 // {
 //     for(int i=0;i<10000;i++)
 //     {
 //         if(!excNum1[i][j].empty())
 //         {
 //             this->Trans[j].excNum.push_back(excNum1[i][j]);
 //             this->Trans[j].occOrb.push_back(occOrb1[i][j]);
 //             this->Trans[j].virOrb.push_back(virOrb1[i][j]);
 //             this->Trans[j].weight.push_back(weight1[i][j]);
 //             this->Trans[j].dipx.push_back(dipx1[i][j]);
 //             this->Trans[j].dipy.push_back(dipy1[i][j]);
 //             this->Trans[j].dipz.push_back(dipz1[i][j]);
 //         } 
 //     }
 // }
}



void ADF_File::calcCollect()
{
  cout << "using calcCollect() " << endl;
  ofstream outFile;
  outFile.open("collectivity.txt");
  for(int i=0; i<positionCounts.size();i++){outFile<<"position counts size " << positionCounts[i] << endl;}
  cout << "position counts size " << positionCounts.size()<< endl;
  int j=0; int betweenIrreps=0;
  double collect=0.000;
  double summing=0.000;
  for(int i=0;i<IrrepWithExc;i++)
  {
      betweenIrreps+=positionCounts[i];
      cout << dipoleAllowedIrreps[i] <<endl;
      outFile<< i << " " << dipoleAllowedIrreps[i] << endl;
      
      for(j;j<betweenIrreps;j++)
      {
          for(int k=0;k<SAOSymmetry.size();k++)
          {
              cout << SAOSymmetry[k] << "  " << occOrb[j] << endl;
              if(SAOSymmetry[k] == occOrb[j])
              {
                  summing += stod(weight[j]); 
              } 
          }
          collect += stod(weight[j]) * stod(weight[j]);
          while(excNum[j] == excNum[j+1])
          {
              j++;
              collect += stod(weight[j]) * stod(weight[j]);
              for(int k=0;k>SAOSymmetry.size();k++)
              {
                  if(SAOSymmetry[k] == occOrb[j])
                  {
                      summing += stod(weight[j]); 
                  } 
              }
          }          
          collectivity.push_back(1/collect);
          SAOCoeffSum.push_back(summing);
          outFile << j << "  " << excNum[j]  <<"  "<< 1/collect << "  " << summing << endl;   
          summing = 0.000;
          collect = 0.000;
         // cout << j << "  " << excNum[j] <<"  "<< occOrb[j] 
         //      <<"  "<< virOrb[j] <<"  "<<weight[j] 
         //      <<"  "<<dipx[j]<<"  "<<dipy[j] 
         //      <<"  "<<dipz[j] << endl;   
      }
      outFile << "\n\n";
      //cout << "\n\n";
  }

  outFile.close();


}



void ADF_File::showDipAdd()
{
  cout<< "using showDipAdd() " << endl;
  ofstream outFile;
  outFile.open("dipoleAdditivity.txt");
  for(int i=0; i<dipadd.size() ; i++)
  {
    outFile<< i << "  " << dipadd[i] << endl;
    //cout << i << "  " << dipadd[i] << endl;
  }
  outFile.close();
}

void ADF_File::calcDipAdd()
{
  cout << "using calcDipAdd() " << endl;
  int j=0; int betweenIrreps=0;
  vector<double> tempdip;

  for(int i=0; i<IrrepWithExc; i++)
  {
      betweenIrreps += positionCounts[i];
      cout << dipoleAllowedIrreps[i] << endl;
      double absdipsumx=0; double absdipsumy=0; double absdipsumz=0;
      double dipsumx=0; double dipsumy=0; double dipsumz=0;
      double sum; double abssum;
      for(j; j<betweenIrreps; j++)
      {
          dipsumx=0; dipsumy=0; dipsumz=0;
          absdipsumx=0; absdipsumy=0; absdipsumz=0;
          absdipsumx+=abs(stod(weight[j]) * stod(dipx[j])); absdipsumy+=abs(stod(weight[j]) * stod(dipy[j])); absdipsumz+=abs(stod(weight[j]) * stod(dipz[j]));
          dipsumx+=stod(weight[j]) * stod(dipx[j]); dipsumy+=stod(weight[j]) * stod(dipy[j]); dipsumz+=stod(weight[j]) * stod(dipz[j]);
          //cout << excNum[j] << "  " << dipx[j] << "  " << dipy[j] << "  " << dipz[j] << endl;
          string Num = excNum[j]; string Num2 = excNum[j+1];
          if(j+1 < betweenIrreps)
          {
          //cout <<excNum[j] << "  " << excNum[j+1] << endl;
          if((excNum[j]) != (excNum[j+1]))
          {
              dipsumx+=stod(weight[j]) * stod(dipx[j]); dipsumy+=stod(weight[j]) * stod(dipy[j]); dipsumz+=stod(weight[j]) * stod(dipz[j]);
              absdipsumx+=abs(stod(weight[j]) * stod(dipx[j])); absdipsumy+=abs(stod(weight[j]) * stod(dipy[j])); absdipsumz+=abs(stod(weight[j]) * stod(dipz[j]));
              cout <<"Final contribution and sum " <<  excNum[j] << "  " << dipx[j] << "  " << dipy[j] << "  " << dipz[j] << "  " << dipsumx <<"  "<< dipsumy << "  "<< dipsumz<<  endl;
          //    dipsumx=0; dipsumy=0; dipsumz=0;
          }
          while((excNum[j]) == (excNum[j+1]))
          {
              j++;
              dipsumx+=stod(weight[j]) * stod(dipx[j]); dipsumy+=stod(weight[j]) * stod(dipy[j]); dipsumz+=stod(weight[j]) * stod(dipz[j]);
              absdipsumx+=abs(stod(weight[j]) * stod(dipx[j])); absdipsumy+=abs(stod(weight[j]) * stod(dipy[j])); absdipsumz+=abs(stod(weight[j]) * stod(dipz[j]));
              cout <<"It is equal " << excNum[j] << "  " << dipx[j] << "  " << dipy[j] << "  " << dipz[j] << endl;
          }
          sum = (dipsumx + dipsumy + dipsumz); 
          abssum = (absdipsumx + absdipsumy + absdipsumz);
          tempdip.push_back(abs(sum/abssum)*100);
          //cout <<"Final contribution and sum " <<  excNum[j] << "  " << dipx[j] << "  " << dipy[j] << "  " << dipz[j] << "  " << 
          //       (dipsumx + dipsumy + dipsumz)/(absdipsumx + absdipsumy + absdipsumz) << endl;
          }
      }

  }   
  cout << "done with dipole additivity " << endl;
  for(int i=0; i<tempdip.size();i++)
  {
      this->dipadd.push_back(tempdip[i]);
  }
  //this->dipadd = tempdip;

}

void ADF_File::showTransitions()
{
  cout << "using showTransitions() " << endl;
  ofstream outFile;
  outFile.open("transitions.txt");
  for(int i=0; i<positionCounts.size();i++){outFile<<"position counts size " << positionCounts[i] << endl;}
  //cout << "position counts size " << positionCounts.size()<< endl;
  int j=0; int betweenIrreps=0;
  for(int i=0;i<IrrepWithExc;i++)
  {
      betweenIrreps+=positionCounts[i];
      //cout << dipoleAllowedIrreps[i] <<endl;
      outFile<< i << " " << dipoleAllowedIrreps[i] << endl;
      for(j;j<betweenIrreps;j++)
      {
          
          outFile << j << "  " << excNum[j] <<"  "<< occOrb[j] 
               <<"  "<< virOrb[j] <<"  "<<weight[j] 
               <<"  "<<dipx[j]<<"  "<<dipy[j] 
               <<"  "<<dipz[j] << endl;   
         // cout << j << "  " << excNum[j] <<"  "<< occOrb[j] 
         //      <<"  "<< virOrb[j] <<"  "<<weight[j] 
         //      <<"  "<<dipx[j]<<"  "<<dipy[j] 
         //      <<"  "<<dipz[j] << endl;   
      }
      outFile << "\n\n";
      //cout << "\n\n";
  }

  outFile.close();
}


void ADF_File::Find_MOS(vector<string> FileLines)
{
  cout << "using Find_MOS " << endl;
  cout << "Looking for Molecular orbitals" << endl;
  string MOsID="List of all MOs, ordered by energy, with the most significant SFO gross populations";
  string MOsEnd=" ";
  string orbitalSymm;
  string goru,finSym;
  vector<string> symmorb, percent,  orbtype;
  string temp;
  int count=0;
  double SSum, PSum;  
  for(int i=0; i<FileLines.size(); i++)
  {
    if(FileLines[i].find(MOsID) != string::npos)
        {
            cout << "Found Molecular orbitals " << endl;
            //cout << FileLines[i] << endl;
            i = i+13;
            while(FileLines[i].find(MOsEnd) != string::npos)
            {
                vector<string> MOLine;
                string element;
                istringstream MOInfo(FileLines[i]);
                while(MOInfo >> element)
                {
                    MOLine.push_back(element);
                }
                cout << FileLines[i] << endl;
                if(MOLine.size() > 9)
                {
                    for(int elem=0; elem<MOLine.size(); elem++)
                    {
                        cout << MOLine[elem] << "  ";
                    }
                    cout << endl;
                    cout << "dipole allowed irreps size " << dipoleAllowedIrreps.size() << endl;
                    if(dipoleAllowedIrreps[0] != "A'")
                    {
                        transform(MOLine[3].begin(), MOLine[3].end(), MOLine[3].begin(), ::tolower); //May cause problems with matching orbitals
                        cout <<" Made it here " << endl;
                    }
                    orbitalSymm = MOLine[2]+MOLine[3];
                    cout << " Made it here " << endl; 
//                    if(Irreps[0] == "S+.u" && count > 0)
//                    {
//                        if(SSum > 0 && SSum > PSum)
//                        {
//                            cout << "decison " << "S+" << goru<< endl;
//                           // percent.push_back(to_string(SSum)); orbtype.push_back(MOLine[5]+MOLine[6]);
//                           // symmorb.push_back(finSym);//"S+"+goru);
//                        }
//                        if(PSum > 0 && PSum > SSum)
//                        {
//                            cout << "decison " << "P+" << goru<< endl;
//                           // percent.push_back(to_string(PSum)); orbtype.push_back(MOLine[5]+MOLine[6]);
//                           // symmorb.push_back(finSym);//"P+"+goru);
//                        }
//                    }
//                    cout << "new orb " << "\n\n";
//                    SSum=0;
//                    PSum=0;
//                    transform(MOLine[3].begin(), MOLine[3].end(), MOLine[3].begin(), ::tolower);
//                    orbitalSymm=MOLine[2]+MOLine[3];
//                    cout << orbitalSymm << " OrbSymm"<< endl;
//                    if(Irreps[0] == "S+.u")
//                    {
//                        if(MOLine[6] == "S")
//                        {
//                            temp = MOLine[4]; temp.pop_back();
//                            cout << "SSum check " << SSum << endl;
//                            cout << "temp " << temp << endl;
//                            SSum+=stod(temp);
//                            //cout << "SSum check " << SSum << endl;
//                            //cout << "temp " << temp << endl;
//                            int dot_find=0;
//                            int col_find=MOLine[3].size();
//                            for(int dot=0;dot<MOLine[3].size();dot++)
//                            {
//                                if(MOLine[3][dot] =='.')
//                                {
//                                    dot_find=dot;
//                                }
//                                if(MOLine[3][dot] ==':')
//                                {
//                                    col_find=dot;
//                                }
//                            }
//                            cout << col_find << "  " << MOLine[3].size() << endl;
//                            goru = MOLine[3];
//                            goru.erase(col_find,MOLine[3].size());
//                            goru.erase(0,dot_find);
//                            finSym="S+"+goru;
//                            cout << MOLine[2]+finSym +"  "+ temp +"  "+ MOLine[5]+MOLine[6] << endl;
//                            //i++;
//                        }
//                        if(MOLine[6] == "P:z" || MOLine[6] == "P:x" || MOLine[6] == "P:y")
//                        {
//                            temp = MOLine[4]; temp.pop_back();
//                            PSum+=stod(temp);
//                            int dot_find=0;
//                            int col_find=MOLine[3].size();
//                            for(int dot=0;dot<MOLine[3].size();dot++)
//                            {
//                                if(MOLine[3][dot] =='.')
//                                {
//                                    dot_find=dot;
//                                }
//                                if(MOLine[3][dot] ==':')
//                                {
//                                    col_find=dot;
//                                }
//                            }
//                            goru = MOLine[3];
//                            goru.erase(col_find,MOLine[3].size());
//                            goru.erase(0,dot_find);
//                            finSym="P+"+goru;
//                            cout << finSym +"  "+ temp +"  "+ MOLine[5]+MOLine[6] << endl;
//                            //i++;
//                        }
//                        else
//                        {
//                            temp = MOLine[4]; temp.pop_back();
//                            symmorb.push_back(orbitalSymm); percent.push_back(temp); orbtype.push_back(MOLine[5]+MOLine[6]);
//                        }
//                    }
//		    else
//		    {
                        temp = MOLine[4]; temp.pop_back();
                        symmorb.push_back(orbitalSymm); percent.push_back(temp); orbtype.push_back(MOLine[5]+MOLine[6]);
//                  }
                    cout << orbitalSymm +"  "+ temp +"  "+ MOLine[5]+MOLine[6] << endl;
                }
                if(MOLine.size() < 9)
                {
 //                   count++;
 //                   //cout << "Irrep " << Irreps[0] << " MOLine[2] " << MOLine[2] << endl;
 //                   if(Irreps[0] == "S+.u")
 //                   {
 //                       if(MOLine[2] == "S")
 //                       {
 //                           temp = MOLine[0]; temp.pop_back();
 //                           SSum += stod(temp);
 //                           cout << "SSum " << SSum << endl;
 //                           //i++;
 //                       }
 //                       if(MOLine[2] == "P:z" || MOLine[2] == "P:x" || MOLine[2] == "P:y")
 //                       {
 //                           temp = MOLine[0]; temp.pop_back();
 //                           PSum += stod(temp);
 //                           cout << "PSum " << PSum << endl;
 //                           //i++;
 //                       }
 //                       else
 //                       {
 //                           cout << "something else inside of line condition " <<endl;
 //                           temp = MOLine[0]; temp.pop_back();
 //                           cout << orbitalSymm <<"  "<< temp <<"  "<< MOLine[1]+MOLine[2] << endl;
 //                           symmorb.push_back(orbitalSymm); percent.push_back(temp); orbtype.push_back(MOLine[1]+MOLine[2]);
 //                       }
 //                 }
 //                 else
 //                 {
                        cout << "something else " <<endl;
                        temp = MOLine[0]; temp.pop_back();
                        cout << orbitalSymm <<"  "<< temp <<"  "<< MOLine[1]+MOLine[2] << endl;
                        symmorb.push_back(orbitalSymm); percent.push_back(temp); orbtype.push_back(MOLine[1]+MOLine[2]);
 //                 }
                }
                i++;

            }
        }
  }
//  if(SSum != 0 && SSum > PSum)
//  {
//      cout << "decison " << "S+" << goru<< endl;
//  //    symmorb.push_back("S+"+goru);
//  }
//  if(PSum !=0 && PSum > SSum)
//  {
//      cout << "decison " << "P+" << goru<< endl;
//   //   symmorb.push_back("P+"+goru);
//  }
  cout <<"symmorb size " << symmorb.size() << endl;
  cout <<"percent size " << percent.size() << endl;
  cout <<"orbType size " << orbtype.size() << endl;
  //for(int sz=1; sz<2;sz++)
  for(int i=0; i<symmorb.size();i++)
  {
      cout << i << endl;
      cout << (symmorb[i]) << " symm " ;
      cout << (percent[i]) << " percent " ;
      cout << (orbtype[i]) << " orbtype " ;
      this->orbitalSymmetry.push_back(symmorb[i]);
      this->percents.push_back(percent[i]);
      this->orbitalType.push_back(orbtype[i]);
  }

}

void ADF_File::show_MOs()
{
  cout << "using show_MOs()" << endl;
  //for(int sz=1; sz<2;sz++)
  cout << "vec size " << orbitalSymmetry.size() << endl;
  for(int i=0;i<this->orbitalSymmetry.size();i++)
  {
    cout << this->orbitalSymmetry[i] << "  " << this->percents[i] << "  " << this->orbitalType[i] << endl;
  }
}

string ADF_File::dotColFind(string lineElement, string SorP)
{
  string finSym;
  string goru;
  int dot_find=0;
  int col_find=lineElement.size();
  int let_find=0;
  string num;
  for(int dot=0;dot<lineElement.size();dot++)
  {
      if(isalpha(lineElement[dot]))
      {
          let_find++;
      }
      if(let_find == 0)
      {
          if(isdigit(lineElement[dot]))
          {
              num +=lineElement[dot];
          }
          
      }
      if(lineElement[dot] =='.')
      {
          dot_find=dot;
      }
      if(lineElement[dot] ==':')
      {
          col_find=dot;
      }
  }
  cout << col_find << "  " << lineElement.size() << endl;
  goru = lineElement;
  goru.erase(col_find,lineElement.size());
  goru.erase(0,dot_find);
  finSym=num+SorP+goru;
return finSym;
}


void ADF_File::check_MOs()
{
  cout << "IN check_MOs() dipoleAllowedIrreps[0] " << dipoleAllowedIrreps[0] << endl;
  if(dipoleAllowedIrreps[0] == "S+.u" || dipoleAllowedIrreps[0] == ("Pi.u"))
  {
    cout << "using check_MOs() Irreps must have Pi.u or S+.u" << endl;
    for(int i=0; i<orbitalSymmetry.size();i++)
    {
        int counter=0;
        double SSum=0; double PSum=0; double DSum=0;
        while(orbitalSymmetry[i] == orbitalSymmetry[i+1])
        {
            cout <<"printing orbital symmetry, type, and percent" << orbitalSymmetry[i] << "  " <<  orbitalType[i] << "  " << percents[i] << endl;
            if(orbitalType[i].find("S")!=string::npos)
            {
                SSum += stod(percents[i]);
            }
            if(orbitalType[i].find("P")!=string::npos)
            {
                PSum += stod(percents[i]);
            }
            if(orbitalType[i].find("D")!=string::npos)
            {
                DSum += stod(percents[i]);
            }

            //counter++;
            i++;
        }
        if(orbitalType[i].find("S")!=string::npos)
        {
            SSum += stod(percents[i]);
        }
        if(orbitalType[i].find("P")!=string::npos)
        {
            PSum += stod(percents[i]);
        }
        if(orbitalType[i].find("D")!=string::npos)
        {
            DSum += stod(percents[i]);
        }
        cout <<"printing orbital symmetry, type, and percent" << orbitalSymmetry[i] << "  " <<  orbitalType[i] << "  " << percents[i] << endl;
        counter++;
        cout <<"printing final orbital symmetry, type, and percent" << orbitalSymmetry[i] << "  " <<  orbitalType[i] << " S " << SSum << " P " << PSum << " D " << DSum  << endl;


        //for(int j=counter;j<=counter;j++)
        //{
        cout << "SSum " << SSum << " PSum " << PSum << endl;
            if(SSum > PSum)// && SSum >= 0.70)
            {
                string temp=dotColFind(orbitalSymmetry[i],"s+");
                cout << "predominantly S " << temp <<endl;
                //this->orbitalSymmetry[i]=temp;
                this->SAOSymmetry.at(counter)=temp;
                this->SAOSymmetry.push_back(temp);
                cout << "Orbital Symmetry " << this->SAOSymmetry.at(counter) << endl;
                cout << orbitalSymmetry[i] << endl;

            }
            if(PSum > SSum)// && PSum >= 0.70)
            {
                string temp=dotColFind(orbitalSymmetry[i],"pi");
                cout << "predominantly P " << temp <<endl;
                //this->orbitalSymmetry[i] = temp;
                this->SAOSymmetry.at(counter)=temp;
                this->SAOSymmetry.push_back(temp);
                cout << "Orbital Symmetry " << this->SAOSymmetry.at(counter) << endl;
                cout << orbitalSymmetry[i] << endl;
            }

        //}
        cout << orbitalSymmetry[i] << " "<< SSum << " SSum " << counter << " counter "<< endl; 
    }
      

  }

//  for(int i=0; i<orbitalSymmetry.size();i++)
//  {
//      cout << (orbitalSymmetry[i]) << " symm " ;
//      cout << (percents[i]) << " percent " ;
//      cout << (orbitalType[i]) << " orbtype " << endl;
//  }

}


void ADF_File::show_SuperAtomic()
{
  cout << "using show_SuperAtomic() " << endl;
  ofstream outFile;
  outFile.open("superatomic.txt");
  vector<string> SAOSymmetryprocessed;
  vector<double> SAOSumprocessed;
  for(int i=0; i<SAOSymmetry.size();i++)
  {
      for(int j=0; j<SAOSymmetry.size();j++)
      {
          if(SAOSymmetry[i] == SAOSymmetry[j] && i != j)
          {
              SAOSymmetry.erase(SAOSymmetry.begin()+j);
              SAOSum.erase(SAOSum.begin()+j);
          }
      }
  }
  for(int i=0; i<this->SAOSymmetry.size();i++)
  {
    outFile << SAOSymmetry[i] << "  " << SAOSum[i] << endl;
    //cout << SAOSymmetry[i] << "  " << SAOSum[i] << endl;
  }
  outFile.close();
}


void ADF_File::calcSuperAtomic(vector<string> supOrb)
{
  cout << "using calcSuperAtomic() " << endl;
  double percentSto;
  double allSto;
  double percentSumm=0;
  double AllPercent=0;
  for(int i=0; i<this->orbitalSymmetry.size();i++)
  {
  cout<<"in calcSAO " << orbitalSymmetry[i] << orbitalSymmetry[i+1]<< endl;
  cout<<"in calcSAO " << percentSumm << endl;
  percentSumm=0;
          while(this->orbitalSymmetry[i] == this->orbitalSymmetry[i+1]) 
          {
		  AllPercent += abs(stod(this->percents[i]));
		  for(int k=0;k<supOrb.size();k++)
              {
                  cout <<"sup orp " << supOrb[k] << " orbtype " << orbitalType[i] << " percent " << percents[i]  << endl;
                  if(this->orbitalType[i] == supOrb[k])
                  {
                      cout << "percent " << percents[i] << endl;
                      percentSumm += stod(this->percents[i]);
                      cout <<"if 1 " <<  this->orbitalSymmetry[i] << "  "<< this->orbitalSymmetry[i+1] << "  "<<this->orbitalType[i]<< "  " << percentSumm << "  total percents " << AllPercent << endl; 
                  }
              }
              i++;
          }
          //if(this->orbitalSymmetry[i] != this->orbitalSymmetry[i+1] )
          //{
              AllPercent += abs(stod(this->percents[i]));
              for(int k=0;k<supOrb.size();k++)
              {
                  if(this->orbitalType[i] == supOrb[k])
                  {
                      cout << "percent no condition" << percents[i] << endl;
                      percentSumm += stod(this->percents[i]);
//                      //cout << this->orbitalSymmetry[i] << "  "<< this->orbitalSymmetry[i+1] << "  "<<this->orbitalType[i]<< "  " << percentSumm << endl;
//                      this->SAOSymmetry.push_back(orbitalSymmetry[i]);
//                      this->SAOSum.push_back(percentSumm);
//                      cout << "if 2 final sum cond1 " << this->orbitalSymmetry[i] << "  " << orbitalType[i] << "  "<< percentSumm << " Total percent " << AllPercent << endl;
                      //percentSto=percentSumm;
                      //percentSumm = 0;
                  }
              }
//              allSto = AllPercent;
//              AllPercent=0;
//              cout << "symm  " << orbitalSymmetry[i] << " part " << percentSto << " whole " << allSto << " ratio " << percentSto/allSto << endl;
//              this->SAOSymmetry.push_back(orbitalSymmetry[i]);
//              this->SAOSum.push_back(percentSto/allSto);
//              percentSto=0; allSto=0;
          //}
//          else if(this->orbitalSymmetry[i] != this->orbitalSymmetry[i+1] && this->orbitalSymmetry[i] != this->orbitalSymmetry[i-1])
//               {
//               percentSumm=0;
//               percentSumm+= stod(this->percents[i]);
//               //cout<<"else if  " << this->orbitalSymmetry[i] << "  " << percentSumm << endl;
//               //cout<<"set to zero" << endl;
//               percentSumm=0;
//               }
//          else if(percentSumm != 0)
          if(percentSumm >=0.0)
              {
//                  this->SAOSymmetry.push_back(this->orbitalSymmetry[i]);
                  //this->SAOSum.push_back(percentSumm);
                  //this->AllOrbSum.push_back(AllPercent);
                  //cout <<"else final sum cond2 " <<  this->orbitalSymmetry[i] << "  " << percentSumm << " All orbs " << AllPercent << endl;
                  percentSto = percentSumm;
                  percentSumm=0;
                  //cout << "under else"<< endl;
                  allSto = AllPercent;
                  AllPercent=0;
                  cout << "Symmetry " << orbitalSymmetry[i]<< "part " << percentSto << " whole " << allSto << " ratio " << percentSto/allSto << endl;
                  if(orbitalSymmetry[i].find("AAA")){orbitalSymmetry[i] = replaceSubstring(orbitalSymmetry[i], "AAA", "a''");}
                  if(orbitalSymmetry[i].find("AA")){orbitalSymmetry[i] = replaceSubstring(orbitalSymmetry[i], "AA", "a'");}
                  if(orbitalSymmetry[i].find("A")){orbitalSymmetry[i] = replaceSubstring(orbitalSymmetry[i], "A", "a");}
                  if(orbitalSymmetry[i].find(":"))
                  {
                    size_t index=0;
                    index=orbitalSymmetry[i].find(":",index);
                    int j=index;
                    while(index<orbitalSymmetry[i].size())
                    {
                      orbitalSymmetry[i].erase(j); j++;
                    }
                  }
                  this->SAOSymmetry.push_back(this->orbitalSymmetry[i]);
                  this->SAOSum.push_back(percentSto/allSto);
                  percentSto=0; allSto=0;
              }
              AllPercent=0;
  }
}

string ADF_File::replaceSubstring(string element, string old , string newstring )
{
  size_t index=0;
  while(true)
  {
      index = element.find(old, index);
      if(index == std::string::npos){break;}
      element.replace(index, newstring.size(), newstring);
      index += 3;
  }
  return element;
}



#include <iostream>
#include <string>
#include <sstream>
#include <fstream> 
#include <vector>
#include <thread>
#include<bits/stdc++.h>//insert function

using namespace std;

class READ_TAPE{
  private:

  public:
  vector<string> fileLines;
  vector<string> tapeLines;
  //vector<string> dipoleMomentVector;
  //vector<double> dipoleX;
  //vector<double> dipoleY;
  //vector<double> dipoleZ;
  vector<double> eigenvectorVector;
  int dipoleLineNum=0;
  int numOfMOs=0;
  int numOfBasisFunctions=0;
  vector<int> numOfEigenvectors;
  vector<double> additivities;


  //CLASS FUNCTIONS DEFINED
  void read_outputFile(string outFileName); 
  void read_TAPEFileNew(string tapeFileName);
  void read_TAPEFileOld(string tapeFileName);
  void read_DipoleMoments();
  void read_Eigenvectors(int state);
  void read_ComputedEigenvectors();
  void calc_dipadd(string tapeFileName, int state);


  //REPRODUCING AGAIN FUNCTIONS AND VARIABLES
  vector<string> symmetries;
  //vector<double> symOcc;
  //vector<double> symVir;
  vector<double> occEnergies;
  vector<double> virEnergies;
  double occSum=0;
  double virSum=0;
  double dipolex, dipoley, dipolez;



  void read_symm(string outFileName);
  void read_orbitals(string outFileName);
  void read_dip_elem(int occ_index, int vir_index, string TAPE);


  //REPRODUCING AGAIN FUNCTIONS AND VARIABLES
  

  vector<string> splitLine(string someLine);
  vector<double> appendVec(string someLine, vector<double> someVector);

  //CLASS CONSTRUCTOR
  READ_TAPE()
  {
    cout << "READ TAPE CONSTRUCTOR " << endl;
  }

  READ_TAPE(string outFileName, string tapeFileName, vector<double> TAPEADDI,int i_state, int f_state)
  {
    read_symm(outFileName);
    read_orbitals(outFileName);
    read_TAPEFileOld(tapeFileName);
    read_ComputedEigenvectors(); 
    for(int state=i_state;state<f_state;state++)
    {
        calc_dipadd(tapeFileName,state);
    }
    ofstream TAPEout;
    TAPEout.open("Tapeadditivities"+to_string(i_state)+"_"+to_string(f_state)+".out");  
    for(int i=0;i<additivities.size(); i++)
    {
        TAPEout << additivities[i] << endl;
    }
    
    cout << "additivities size " << additivities.size() << endl; 
//    exit(3);
//    read_outputFile(outFileName);
//    read_TAPEFileOld(tapeFileName);
//    read_TAPEFileNew(tapeFileName);
//    read_DipoleMoments();
//    read_ComputedEigenvectors(); 
//
//    for(int state=0;state<5;state++)//numOfEigenvectors.size();state++)
//    {
//        read_Eigenvectors(state);//numOfEigenvectors[state]);
//	        calc_dipadd(tapeFileName,state);
//    }
//    ofstream TAPEout;
//    TAPEout.open("Tapeadditivities"+to_string(i_state)+"_"+to_string(f_state)+".out");  
//    for(int i=0;i<additivities.size(); i++)
//    {
//        TAPEout << additivities[i] << endl;
//    }  
//    cout << "additivities size " << additivities.size() << endl; 
    
//    //double addi[additivities.size()];
//    ofstream TAPEout;
//    TAPEout.open("Tapeadditivities.out");  
//    for(int i=0;i<additivities.size(); i++)
//    {
//        TAPEout << additivities[i] << endl;
//    }
//    
//    cout << "additivities size " << additivities.size() << endl; 
    //dipoleRearrange();

  }

};


void READ_TAPE::read_dip_elem(int occ_index, int vir_index, string TAPE)
{
  int nelem = virSum * occSum;
  int index = occ_index + vir_index * occSum;
  //double dipolex, dipoley, dipolez;
  int dipline=0;
  string GrabLine;
 // string runner;
 // ifstream File;
 // File.open(TAPE, ios::in);
 // if(!File){cout << "Unable to open TAPE File to read dipole moments" << TAPE << " !! \n"; exit(1) ;}
  cout << "reading dipole moment" << endl;
  cout << virSum << " " << occSum << " " << occ_index << " " << vir_index << endl;
  cout << tapeLines.size() << endl;
  for(int fi=1; fi<tapeLines.size();fi++)
  {
      //cout << "made it here"<< endl; 
      //GrabLine = tapeLines[fi];
      if(tapeLines[fi].find("dipole elements") != string::npos)
      {
          //cout << tapeLines[fi] << endl;
          fi=fi+2; 
          //cout << "made it here"<< endl; 
          for(int cor=1; cor<4; cor++)
          {
              int elemLine = index/3;
              int elemPos = index%3;
              fi = fi+(elemLine-dipline);
              //for(int lin=0; lin<elemLine - dipline; lin++)
              //{
              //    getline(File,runner);
              //    cout << runner << endl;
              //}
                  //cout << "elemLine " << elemLine << "  elemPos " << elemPos << endl;
                  dipline = elemLine;
                  index += nelem;
                  istringstream dipoleLine(tapeLines[fi]);
                  //cout << "tape Line in dipole elem " << tapeLines[fi] << endl;
                  string element;
                  vector<double> splitDipole;
                  while( dipoleLine >> element)
                  {
                      //cout << "line element " << element << endl;
                      splitDipole.push_back(stod(element));
                  }
                  //dipLine = elemLine;
                  //index += (occOrbs) * (virOrbs);
                  if(cor%3 == 0)
                  {
                      dipolex=(splitDipole[elemPos]);
                      cout << "dip x " << splitDipole[elemPos]<< endl;;
                  }
                  
                  if(cor%3 == 1)
                  {
                      dipoley=(-1*splitDipole[elemPos]);
                      cout << "dip y " << splitDipole[elemPos] << endl;
                  }

                  if(cor%3 == 2)
                  {
                      dipolez=(splitDipole[elemPos]);
                      cout << "dip z " << splitDipole[elemPos]<< endl;
                  }
          }
      }
      
  }


}


void READ_TAPE::read_symm(string outFileName)
{
  vector<int> linesOfInterest;
  string GrabLine;
  string runner;
  ifstream File;
  File.open(outFileName, ios::in);

  if(!File){cout << "Unable to open " << outFileName << " !! \n"; exit(1) ;} /*check if the file is open*/
  cout << outFileName << " is open" << endl;
  int counter = 0;
  while(getline(File,runner))
  {
      GrabLine = runner;
      fileLines.push_back(GrabLine);
      if(GrabLine.find("Irreducible Representations, including subspecies") != string::npos)
      {
          cout << GrabLine << endl;
          linesOfInterest.push_back(counter);
      }
      counter++;
  }
  File.close();
  
  vector<string> lineVec;
  for(int a=linesOfInterest[linesOfInterest.size()-1]+2; a < fileLines.size(); a++)
  {
      if(fileLines[a].size() > 1)
      {
//          cout << fileLines[a] << endl;
          string element;
          istringstream line(fileLines[a]);
          while(line >> element)
          {
              if(element.find(":") != string::npos)
              {
                  size_t found = element.find(":");
                  element.erase(found,element.length());
              }
              cout << element << "  " ;
              lineVec.push_back(element);
          }
          cout << endl;
          //if(numOfMOs<stoi(lineVec[1]))
          //{
          //    numOfMOs = stoi(lineVec[1]);
          //}
      }
      else{break;}
  }
  this->symmetries = lineVec;
//  for(int i=0; i<symmetries.size(); i++)
//  {
//    for(int j=0; j<symmetries.size(); j++)
//    {
//        if(symmetries[i] == symmetries[j] && i!=j)
//        {
//            symmetries[j].erase();
//        }
//    }
//  }
  cout << "number of orbitals " << numOfMOs << endl;

}


void READ_TAPE::read_orbitals(string outFileName)
{
  vector<int> linesOfInterest;
  string GrabLine;
  string runner;
  ifstream File;
  File.open(outFileName, ios::in);

  if(!File){cout << "Unable to open " << outFileName << " !! \n"; exit(1) ;} /*check if the file is open*/
  cout << outFileName << " is open" << endl;
  int counter = 0;
  while(getline(File,runner))
  {
      GrabLine = runner;
      fileLines.push_back(GrabLine);
      if(GrabLine.find("Orbital Energies, all Irreps") != string::npos)
      {
          cout << GrabLine << endl;
          linesOfInterest.push_back(counter);
      }
      counter++;
  }
  File.close();
  
  
  vector<double> symOcc(symmetries.size());
  vector<double> symVir(symmetries.size());
  //vector<double> occEnergies;
  //vector<double> virEnergies;
  //double occSum=0;
  //double virSum=0;
  for(int a=linesOfInterest[linesOfInterest.size()-1]+5; a < fileLines.size(); a++)
  {
      if(fileLines[a].size() > 1)
      {
//          cout << fileLines[a] << endl;
          string element;
          vector<string> lineVec;
          istringstream line(fileLines[a]);
          while(line >> element)
          {
              cout << element << "  " ;
              lineVec.push_back(element);
          }
          cout << endl;
      
          double occupation;
          occupation = stod(lineVec[2]);
          if(occupation > 0.001)
          {
              for(int j=0; j<symmetries.size();j++)
              {
                  if(symmetries[j] == lineVec[0])
                  {
                      symOcc[j] += 1;
                      occEnergies.push_back(stod(lineVec[4]));
                  }
              }
          }
          else
          {
              for(int j=0; j<symmetries.size();j++)
              {   
                  if(symmetries[j] == lineVec[0])
                  {
                      symVir[j] += 1;
                      virEnergies.push_back(stod(lineVec[4]));
                  }
              }
          }
          
      
          //if(numOfMOs<stoi(lineVec[1]))
          //{
          //    numOfMOs = stoi(lineVec[1]);
          //}
      }
      else{break;}
  }
  for(int i=0; i<symOcc.size();i++)
  {
      occSum += symOcc[i]; 
      virSum += symVir[i];
      cout <<"symmetry "<< symmetries[i] <<" occ Sym " << symOcc[i] << " vir Sym " << symVir[i] << endl;
  }
  
  cout << "number of Occ orbitals " << occSum << "  Vir orbitals " << virSum << endl;
  
}


void READ_TAPE::read_outputFile(string outFileName)
{
  vector<int> linesOfInterest;
  string GrabLine;
  string runner;
  ifstream File;
  File.open(outFileName, ios::in);

  if(!File){cout << "Unable to open " << outFileName << " !! \n"; exit(1) ;} /*check if the file is open*/
  cout << outFileName << " is open" << endl;
  int counter = 0;
  while(getline(File,runner))
  { 
      GrabLine = runner;
      fileLines.push_back(GrabLine);
      if(GrabLine.find("Orbital Energies, all Irreps") != string::npos)
      {
          cout << GrabLine << endl;
          linesOfInterest.push_back(counter);
      }
      counter++;
  }
  File.close();

  for(int a=linesOfInterest[linesOfInterest.size()-1]+5; a < fileLines.size(); a++)
  {
      if(fileLines[a].size() > 1)
      {
//          cout << fileLines[a] << endl;
          string element;
          vector<string> lineVec;
          istringstream line(fileLines[a]);
          while(line >> element)
          {
              //cout << element << "  " ;
              lineVec.push_back(element);
          }
          //cout << endl;
          if(numOfMOs<stoi(lineVec[1]))
          {
              numOfMOs = stoi(lineVec[1]);
          }

      }
      else{break;}
  }
  cout << "number of orbitals " << numOfMOs << endl;

}

vector<string> READ_TAPE::splitLine(string someLine)
{
    string element;
    vector<string> line;
    istringstream lineStream(someLine);
    while(lineStream >> element)
    {
        line.push_back(element);
    }
    return line;
}

vector<double> READ_TAPE::appendVec(string someLine, vector<double> someVector)
{
    string element;
    istringstream lineStream(someLine);
    while(lineStream >> element)
    {
        someVector.push_back(stod(element));
    }
    return someVector;
}

void READ_TAPE::read_TAPEFileOld(string tapeFileName)
{
  vector<int> linesOfInterest;
  vector<string> Lines;
  string GrabLine;
  string runner;
  ifstream File;
  File.open(tapeFileName, ios::in);

  if(!File){cout << "Unable to open " << tapeFileName << " !! \n"; exit(1) ;} /*check if the file is open*/
  cout << tapeFileName << " is open" << endl;
  int counter = 0;
  while(getline(File,runner))
  { 
      GrabLine = runner;
      Lines.push_back(GrabLine);
      if(GrabLine.find("dipole elements") != string::npos)
      {
          //cout << GrabLine << endl;
          linesOfInterest.push_back(counter);
          //cout << counter << endl;
      }
      counter++;
  }
  dipoleLineNum = linesOfInterest[0];
  File.close();
  this->tapeLines = Lines;
} 


void READ_TAPE::read_TAPEFileNew(string tapeFileName)
{
  vector<int> linesOfInterest;
  vector<double> contribLines;
  vector<double> contribs;
  vector<double> transDips;
  string GrabLine;
  string runner;
  string numOfContributions;
  string contribNum;
  ifstream File;
  File.open(tapeFileName, ios::in);

  if(!File){cout << "Unable to open " << tapeFileName << " !! \n"; exit(1) ;} /*check if the file is open*/
  cout << tapeFileName << " is open" << endl;
  int counter = 0;
  int contribCounter = 1;
  int coeffCounter = 1;
  while(getline(File,runner))
  { 
      GrabLine = runner;
      //tapeLines.push_back(GrabLine);
      //if(GrabLine.find("dipole elements") != string::npos)
      //{
      //    //cout << GrabLine << endl;
      //    linesOfInterest.push_back(counter);
      //    //cout << counter << endl;
      //}
      if(GrabLine.find("nr of contributions ") != string::npos)
      {
         //cout << GrabLine << "  ";
         contribLines.push_back(counter);
         vector<string> tempLine;
         tempLine = splitLine(GrabLine);
         contribNum = tempLine[tempLine.size()-1];
         //cout << "Excited State " + contribNum << endl;
         int numOfContributions; 
         getline(File,runner);
         getline(File,runner);
         GrabLine=runner;
         numOfContributions = stoi(GrabLine);
         //cout << "Number of Contributions " << numOfContributions << endl;
         while(getline(File,runner))
         {
             GrabLine = runner;
             if(GrabLine.find("contr "+ contribNum +" ") != string::npos)
             {
                 //cout << GrabLine << endl;
                 getline(File,runner);
                 //getline(File,runner);
                 //GrabLine = runner;
                 //cout << GrabLine << endl;
                 while(getline(File,runner)) 
                 {
                     GrabLine = runner;
                     if(GrabLine.find("Excitations SS") != string::npos)
                     {
                         break;
                     }    
                     //getline(File,runner);
                     //cout << GrabLine << "  " << "vec line" << endl;
                     istringstream lineStream(GrabLine);
                     string element;
                     while(lineStream >> element)
                     {
                         contribs.push_back(stod(element));
                     } 
                 }
             }
             if(GrabLine.find("contr transdip "+ contribNum + " ") != string::npos)
             {
                 //cout << "Transition dipole moment State " + contribNum << endl;
                 getline(File,runner);
                 //getline(File,runner);
                 //GrabLine = runner;
                 //cout << GrabLine << endl;
                 while(getline(File,runner)) 
                 {
                     GrabLine = runner;
                     if(GrabLine.find("Excitations SS") != string::npos)
                     {
                         break;
                     }    
                     //getline(File,runner);
                     //cout << GrabLine << "  vec line" << endl;
                     istringstream lineStream(GrabLine);
                     string element;
                     while(lineStream >> element)
                     {
                         transDips.push_back(stod(element));
                     } 
                 }
               break;
             }
         }
      if(contribNum.size() != 0)
      {
          cout << "\nState "+contribNum << "  "<< contribs.size() << endl;
          for(int i=0;i<contribs.size();i++)
          {
              if(i%3 == 0 && i>0){cout << endl;}
              cout << contribs[i] << "  " ;
          }
          cout << "\nTransition Dipole moments"<< " " << transDips.size() << endl;
          for(int i=0;i<transDips.size();i++)
          {
             if(i%3 == 0 & i>0){cout << endl;}
             cout << transDips[i] << "  " ; 
          }
          cout << endl;
         double xaddabs=0;
         double yaddabs=0;
         double zaddabs=0;
         double xadd=0;
         double yadd=0;
         double zadd=0;
         double totadd;
         double totaddabs;
         for(int i=0; i<contribs.size();i++)
         {
             xaddabs += abs(contribs[i] *transDips[((i*3))]);
             yaddabs += abs(contribs[i] *transDips[(i*3)+1]);
             zaddabs += abs(contribs[i] *transDips[(i*3)+2]);
             xadd +=  contribs[i] *transDips[((i*3))];
             yadd +=  contribs[i] *transDips[(i*3)+1];
             zadd +=  contribs[i] *transDips[(i*3)+2];
             //COUTS HERE WILL HELP IN TROUBLESHOOTING CALCULATIONS
             cout << "contrib:  x, y, z " << contribs[i] <<": "<< transDips[((i*3))] << ", " << transDips[((i*3))+1] << ", "<< transDips[(i*3)+2]<< endl;
             cout << "x add & x add abs " << xadd << " " << xaddabs << " " << endl;
             cout << "y add & y add abs " << yadd << " " << yaddabs << " " << endl;
             cout << "z add & z add abs " << zadd << " " << zaddabs << " " << endl;
         }
         totadd = (  pow(pow(xadd,2) + pow(yadd,2) + pow(zadd,2),0.5) );
         totaddabs = ( pow(pow(xaddabs,2) + pow(yaddabs,2) + pow(zaddabs,2),0.5) );
         cout << "additivity " << totadd/totaddabs << endl;
         additivities.push_back(totadd/totaddabs);  

      }
      contribs.clear();
      transDips.clear(); 
//         getline(File,runner);
//         GrabLine = runner;
//         cout << GrabLine << endl;
//         contribCounter++;
      }
      counter++;
  }
//  dipoleLineNum = linesOfInterest[0];
  File.close();
} 


void READ_TAPE::read_DipoleMoments()
{
//  string lengthString = tapeLines[dipoleLineNum+1];
//  istringstream pik(lengthString);
//  vector<string> pikElements;
//  string element;
//  while(pik >> element)
//  {
//      pikElements.push_back(element);
//  }
//  for(int i=0; i<pikElements.size(); i++)
//  {
//      //cout << pikElements[i] << "  ";
//  }
//  //cout << endl;
//  
//  numOfBasisFunctions = stoi(pikElements[1])/3;
//  int virOrbs = 438;//numOfBasisFunctions/numOfMOs;
//  int occOrbs = 303;//numOfMOs;
//  int extraline = numOfBasisFunctions % 3; 
//  int linesDown = numOfBasisFunctions + extraline + dipoleLineNum;
//  cout << "extra lines " << extraline << endl;
//  cout << "number of basis functions " << numOfBasisFunctions << endl;
//  cout << "number of elements " << occOrbs * virOrbs << " occ orbs " << occOrbs << " vir orbs " << virOrbs << endl;
//  //EXTRACTING DIPOLE MOMENT VECTORS INTO A SINGLE VECTOR
//  for(int a=dipoleLineNum+2; a < linesDown; a++)
//  {
//          //cout << tapeLines[a] << endl; 
//          //string element;
//          //vector<string> lineVec;
//          //istringstream line(tapeLines[a]);
//          //while(line >> element)
//          //{
//              //cout << element << "  " ;
//          dipoleMomentVector.push_back(tapeLines[a]);
//          //}
//          //cout << endl;
//  }
//  for(int i=0; i<occOrbs; i++)
//  {
//      for(int j=0; j<virOrbs; j++)
//      {
//          int dipLine=0;
//          int index= i + j * (occOrbs);
//          //cout << "dipline " << dipLine << endl;
//          //cout << "index  " << index << endl;
//          //cout << "dipole lines size " << dipoleMomentVector.size() << endl;
//          for(int cor=1; cor<4; cor++)
//          {
//              int elemLine = index/3;
//              int elemPos = index%3;
//              cout << "elemLine " << elemLine << "  elemPos " << elemPos << endl;
//              string element;
//              vector<double> splitDipole;
//              istringstream dipoleLine(dipoleMomentVector[elemLine]);
//              while( dipoleLine >> element)
//              {
//                  splitDipole.push_back(stod(element));
//              }
//
//              dipLine = elemLine;
//              index += (occOrbs) * (virOrbs);
//              if(cor%3 == 0)
//              {
//                  dipoleX.push_back(splitDipole[elemPos]);
//                  cout << "dip x " << splitDipole[elemPos]<< endl;;
//              }
//              
//              if(cor%3 == 1)
//              {
//                  dipoleY.push_back(-1*splitDipole[elemPos]);
//                  cout << "dip y " << splitDipole[elemPos] << endl;
//              }
//
//              if(cor%3 == 2)
//              {
//                  dipoleZ.push_back(splitDipole[elemPos]);
//                  cout << "dip z " << splitDipole[elemPos]<< endl;
//              }
//          }
//     }
//  }
//      for(int j=0; j<orbindex;j++)//dipoleMomentVector.size();i++)
//      {
    //      if(i <= numOfBasisFunctions)
    //      {
              //dipoleY.push_back(stod(dipoleMomentVector[i])*-1); 
              //cout << " y dipole " << stod(dipoleMomentVector[i])*-1 << "  " ;
    //      }
    //      else if(i <= numOfBasisFunctions * 2)
    //      { 
              //dipoleZ.push_back(stod(dipoleMomentVector[i+(159)])); 
              //cout << " z dipole " << stod(dipoleMomentVector[i]) << "  ";
    //      }
    //      else if(i <= numOfBasisFunctions *3)
    //      {
              //dipoleX.push_back(stod(dipoleMomentVector[i+(319)]));
              //cout << " x dipole " << stod(dipoleMomentVector[i]) << endl;;
    //      }
      //}
  
//  cout << "size x, y, z " << dipoleX.size() << "  " << dipoleY.size() << "  " << dipoleZ.size() << endl;
}



//void READ_TAPE::dipoleRearrange()
//{
//  for(int i=0; i<dipoleMomentVector.size();i++)
//  {
//      dipoleY.push_back(stod(dipoleMomentVector[i])*-1); i++;
//      dipoleZ.push_back(stod(dipoleMomentVector[i])); i++;
//      dipoleX.push_back(stod(dipoleMomentVector[i]));
//  }
//  cout << "size x, y, z " << dipoleX.size() << "  " << dipoleY.size() << "  " << dipoleZ.size() << endl;
//}

void READ_TAPE::read_ComputedEigenvectors()
{
  int counter=0;
  int lineOfInterest=0;
  for(int i=0;i<tapeLines.size();i++)
  {
      if(tapeLines[i].find("eigenvector ") != string::npos)
      {
          counter++;
//          cout << tapeLines[i]<< " counter number " << counter << " line number " << lineOfInterest<< endl;
          numOfEigenvectors.push_back(counter);
      }
  }
}


void READ_TAPE::read_Eigenvectors(int state)
{
  string lengthString = tapeLines[2];
  istringstream pik(lengthString);
  vector<string> pikElements;
  string element;
  while(pik >> element)
  {
      pikElements.push_back(element);
  }
  for(int i=0; i<pikElements.size(); i++)
  {
      //cout << pikElements[i] << "  ";
  }
  //cout << endl;
  
  numOfBasisFunctions = stoi(pikElements[1])/3;
 // numOfBasisFunctions = stoi(pikElements[1])/3;
  int counter=0;
  vector<double> vecSave;
  int lineOfInterest=0;
  for(int i=0;i<tapeLines.size();i++)
  {
      if(tapeLines[i].find("eigenvector "+to_string(state)+" ") != string::npos)
      {
          cout << "eigenvector " +to_string(state)+" found" << endl;
          istringstream basisNum(tapeLines[i+1]);//WORKING HERE
          
          lineOfInterest = counter;
          //cout << tapeLines[i]<< " counter number " << counter << " line number " << lineOfInterest<< endl;
      }
      counter++;
  }

  int extraline = numOfBasisFunctions % 3; 
  int linesDown = virSum*occSum;//(numOfBasisFunctions/3) + extraline -1;//WORKING HERE 
  //cout << "extra lines " << extraline << endl;
  //cout << "number of basis functions " << numOfBasisFunctions << endl;
  //cout << "lines down " << linesDown << endl;
  //EXTRACTING EIGENVECTORS INTO A SINGLE VECTOR
  int a=0;
  while(a < linesDown)
  {
      //cout << tapeLines[lineOfInterest+2+a] << endl; 
      string element;
      //vector<string> lineVec;
      istringstream line(tapeLines[lineOfInterest+a+2]);
      while(line >> element)
      {
          //cout << element << "  " ;
  //        if( isdigit(element[0]) )
          {
              vecSave.push_back(stod(element));
          }
  //        else
  //        {
  //            cout << "Warning pushing back a zero in Eigenvector read from TAPE File " << endl;
  //            vecSave.push_back(0.00);
  //        }
      }
      //cout << endl;
      a++;
  }
  eigenvectorVector = vecSave;
  cout << "size of eigenvectorVector " << eigenvectorVector.size() << endl;

}


void READ_TAPE::calc_dipadd(string tapeFileName, int state)
{
  int occOrbs = occSum;//numOfMOs;
  int virOrbs = virSum;//dipoleX.size()/numOfMOs;
  double absSumX=0;
  double absSumY=0;
  double absSumZ=0;
  double sumX=0;
  double sumY=0;
  double sumZ=0;
  double totalSum=0;
  double totalAbsSum=0;
  int excstate=1;
  int count=0;
  cout <<"calculating state " << state << endl;
  read_Eigenvectors(numOfEigenvectors[state]);
  cout <<" size of eigenvector " << eigenvectorVector.size() << endl;
  int occ_index=0;
  int vir_index=0;
  int index = 0;
  
      //while(index < (occSum * virSum))
      //{
          //index = occ_index + vir_index * occSum;
          //read_dip_elem(occ_index, vir_index, index, tapeFileName);
          cout << "made it here "<< eigenvectorVector.size()<< endl;
//      for(int j=0; j<virOrbs;j++)//dipoleMomentVector.size();i++)
//      for(int i=0; i< dipoleX.size();i++)
//      {
//    if(i == virOrbs*count)
//    {
//        count++;
  //COUT STATEMENTS IN THIS BLOCK CAN BE HELPFUL FOR TROUBLESHOOTING CALCULATIONS
//        cout << "State num " << count << endl;
//    }
               //std::thread Th1([]{
               for(int i=0;i<eigenvectorVector.size();i++)
               {
               if(abs(eigenvectorVector[i]) > 0.000001)
               {
                     //index = occ_index + vir_index * occSum;
                     read_dip_elem(occ_index, vir_index, tapeFileName);
                   //if(dipolex != 0 || dipoley != 0 || dipolez !=0)
                     //{ 
                     cout <<" occ_index " << occ_index <<" vir_index " << vir_index << " coeff " << eigenvectorVector[i] << endl;
                     absSumX += abs(eigenvectorVector[i] * dipolex);
                     sumX += eigenvectorVector[i] * dipolex;
                     absSumY += abs(eigenvectorVector[i] * dipoley);
                     sumY += eigenvectorVector[i] * dipoley;
                     absSumZ += abs(eigenvectorVector[i] * dipolez);
                     sumZ += eigenvectorVector[i] * dipolez;
                     cout << " dip x " << dipolex << " dip y " << dipoley << " dip z " << dipolez << endl;
                     cout << "Sum " << sqrt(pow(sumX,2) + pow(sumY,2) + pow(sumZ,2)) << " abs Sum " << sqrt(pow(absSumX,2) + pow(absSumY,2) + pow(absSumZ,2)) << endl;
                     //}
               }
               count +=1;
               vir_index +=1;
               if(vir_index >= virSum)
               {
                   vir_index = 0;
                   occ_index +=1;
               }
               }

               //});//END OF THREAD

          index += count;
          //cout << index << endl;
           
      //}     
       
  cout << "State num " << state << endl;
  totalAbsSum = sqrt(pow(absSumX,2) + pow(absSumY,2) + pow(absSumZ,2));
  totalSum = sqrt(pow(sumX,2) + pow(sumY,2) + pow(sumZ,2));
  cout <<"x, y, z sum " << sumX << ",  " << sumY << ",  " << sumZ << endl
       <<"absx, absy, absz sum " << absSumX << ",  " << absSumY << ",  " << absSumZ << endl
       << "additivity " <<  totalSum/ totalAbsSum << endl;
  additivities.push_back(totalSum/totalAbsSum);
  eigenvectorVector.clear();
    
  
}




//int main(int argc, char ** argv)
//{
//  string name = argv[1];
//  string name2 = argv[2];
//  READ_TAPE test(name, name2); 
//  cout << "additivity size " << test.additivities.size() << endl; 
//
//}


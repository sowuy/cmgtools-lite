#include <iostream>


float nJetsnBJets( int nJets, int nBJets){ 
  if (nJets > 4)  nJets = 4;
  if (nBJets > 4) nBJets = 4;
  return nJets + nJets*(nJets-1)/2 + nBJets; 
}


float nJetsnBJetsTTZ( int nJets, int nBJets){ 

  if (nJets < 2) return -1;

  if (nBJets == 0){
    if      (nJets == 2) return 0;
    else if (nJets == 3) return 1;
    else                 return 2; // > 3
  }
  else if (nBJets == 1){
    if      (nJets == 2) return 3;
    else if (nJets == 3) return 4;
    else                 return 5; // > 3
  }
  else if (nBJets >= 2){
    if      (nJets == 2) return 6;
    else if (nJets == 3) return 7;
    else                 return 8; // > 3
  }
  else return -1;


  //  cout << nJets << " " << nBJets << " " << nJets + nJets*(nJets-1)/2 + nBJets<< endl;
  //return nJets + nJets*(nJets-1)/2 + nBJets; 


}

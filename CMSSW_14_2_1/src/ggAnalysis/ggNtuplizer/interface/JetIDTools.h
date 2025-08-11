#ifndef GGNTUPLIZER_JETIDTOOLS_H
#define GGNTUPLIZER_JETIDTOOLS_H

#include <cmath>
#include <string>
#include <iostream>

struct JetIDVars {
  float NHF;  // Neutral Hadron Fraction
  float NEMF; // Neutral EM Fraction
  float MUF;  // Muon Fraction
  float CHF;  // Charged Hadron Fraction
  float CEMF; // Charged EM Fraction
  int NumConst; // Number of Constituents
  int NumNeutralParticle; // Number of Neutral Particles
  int CHM;  // Charged Multiplicity
};

// Declare the Jet ID checking function
bool getJetID(JetIDVars vars,
              std::string jettype = "CHS",
              std::string year = "2018",
              double eta = 0,
              bool tightLepVeto = true,
              bool UltraLegacy = false,
              bool isRun3 = false);

#endif // GGNTUPLIZER_JETIDTOOLS_H


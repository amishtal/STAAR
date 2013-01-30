/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Residue.hpp
//  Date: 16 Jan 2011
//  Version: 1.0
//  Description: Holds the class declaration for Residue
//
/***************************************************************************************************/
//
/***************************************************************************************************/
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//  Redistributions of source code must retain the above copyright notice,
//  this list of conditions and the following disclaimer.
//  Redistributions in binary form must reproduce the above copyright notice,
//  this list of conditions and the following disclaimer in the documentation
//  and/or other materials provided with the distribution.
//  Neither the name of the University of Tennessee nor the names of its contributors
//  may be used to endorse or promote products derived from this software
//  without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
//  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
//  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS
//  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
//  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
//  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
//  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*************************************************************************************************/


#ifndef __RESIDUE_HPP__
#define __RESIDUE_HPP__

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "Atom.hpp"


// Defines to set indexes for the vector that holds plane information
// This is only to make reading the code easier so that when you look
// through the code you don't see just numbers without any explanation
// But, these are for PHE
#define CG_PLANE_COORD_PTT  0 
#define CD1_PLANE_COORD_PTT 1
#define CD2_PLANE_COORD_PTT 2
#define CE2_PLANE_COORD_PTT 2
#define C_2_PLANE_COORD_PTT 2

// And these are for GLU and ASP.  Don't ask why I used different 
// indexes than I did above, because I'm not quite sure
#define C__PLANE_COORD_AG  2
#define O_1_PLANE_COORD_AG 0
#define O_2_PLANE_COORD_AG 1

#define HYDROGEN_BOND_DISTANCE 0.632469

class Residue{
private:
  // calculates the centers of the AA and sets the plane
  // information order to calculate the angle later
  void centerPHEorTYR();
  void centerTRP();

  // Calculates the CM of PHE using alternate locations
  void centerPHEorTYR_altloc();

  // Calculates the center of mass for the PO4, 2HP, and PI ligands
  //void centerPO4();
  void centerPO4or2HPorPI();
  // and the center of charge of PO4, 2HP, and PI ligands
  //void centerPO4();
  void centerPO4or2HPorPI_charge();

  // Calculates the center of mass for the 2PO and PO3 ligands
  void center2POorPO3();
  // and the center of charge of 2PO and PO3 ligands
  void center2POorPO3_charge();

  // DO NOT USE THESE!!!! If you use this, puppies will be harmed!
  // This is only here for legacy reasons because this is what the old
  // STAAR code used.  It is wrong.
  void centerASP();
  void centerGLU();

  // These find the oxygens that will later be used for distances for ASP and GLU
  void centerASP_oxygen();
  void centerASP_oxygen_altloc();
  void centerGLU_oxygen();
  void centerGLU_oxygen_altloc();

  void centerPHEorTYR_simplified();
  void centerASP_charge();
  void centerGLU_charge();

  void printPHEorTYR(FILE* output);
  void printASP(FILE* output);
  void printGLU(FILE* output);

  // These functions mark the alternate locations of the atoms that were
  // not used in the calculation of the center of mass
  void markAltLocAtomsPHEorTYR(int index);
  void markAltLocAtomsASP(int index);
  void markAltLocAtomsGLU(int index);
  void markAltLocAtomsPO4or2HPorPI(int index);
  void markAltLocAtoms2POorPO3(int index);

  void unmarkAltLocAtomsPHEorTYR();
  void unmarkAltLocAtomsASP();
  void unmarkAltLocAtomsGLU();
  void unmarkAltLocAtomsPO4or2HPorPI();
  void unmarkAltLocAtoms2POorPO3();

  string makeConectPHEorTYR();
  string makeConectPHEorTYR_altloc(int c);
  string makeConectASP();
  string makeConectASP_altloc(int c);
  string makeConectGLU();
  string makeConectGLU_altloc(int c);
  string makeConectPO4or2HPorPI();
  string makeConect2POorPO3();
  string makeConectCarbonRing(int c);

public:
  // constructor
  Residue();

  // destructor
  ~Residue();

  // Finds the alternate locations
  void determineAltLoc(vector<char>&altloc_ids);

  bool findCarbonRings();

  // calculate the center of the AA. Calls the individual
  // functions above depending on AA
  void calculateCenter(bool center);

  void calculateAnglesPreHydrogens(Residue aa2,
                                   int index1,
                                   int index2,
                                   float* angle,
                                   float* angle1,
                                   float* angleP);
  bool calculateDistancesAndAnglesPostHydrogens(Residue aa2,
                                                Coordinates closestOxygen,
                                                float threshold,
                                                float* dist,
                                                float* distOxy,
                                                float* distOxy2,
                                                float* angle,
                                                float* angle1,
                                                float* angleP);

  // Prints out only the atoms that we need to create benzene or formate
  void printNeededAtoms(FILE* output);

  // Marks alt loc atoms that aren't used for center calculations
  void markAltLocAtoms(int index);

  // Unmarks alt loc atoms that aren't used for center calculations
  // Needed so that coordinates aren't invalid when they are 
  // looked at as a different pair
  void unmarkAltLocAtoms();

  // Make CONECT lines.  This avoids mono/di-atomic molecules
  string makeConect(int c);

  // Fix to remove excess hydrogens that Babel sometimes adds
  // to the GLU and ASP residues
  bool removeExcessHydrogens(vector<string> conect);

  vector< vector<Atom*> > altlocs;

  // this holds pointers to the ATOM strings
  vector<Atom*> atom;

  // holds the calculated centers that will be examined
  vector<Coordinates> center;

  // Holds the carbon atoms involved in carbon rings.
  //vector< vector<Atom*> > carbonRings;
  // Holds the centers (and plane info) of the carbon rings
  // stored in the above vector.
  //vector<Coordinates> carbonRingCenters;

  // Holds information for carbon rings, grouping alternate
  // locations in a nested vector.
  // So carbonRingCenters[0] is a set of alternate locations
  // for a single carbon ring and may contain only a single
  // element.
  vector< vector< vector<Atom*> > > carbonRings;
  vector< vector<Coordinates> > carbonRingCenters;

  // holds the residue name
  string residue;

  // true of this AA isn't important, false otherwise
  bool skip;

  // true if there is an alt loc in ATOM line
  bool altLoc;

  bool corrected;

  string line;

  friend ostream& operator<<(ostream& output, const Residue& p);

};

#endif

/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Residue.cpp
//  Date: 16 Jan 2011
//  Date Modified: 7 Feb 2011
//  Version: 1.0
//  Description: Class implementations for Residue
//
//  Updates: Fixed the non-center of charge calculations for ASP and GLU (7 Feb 2011)
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

#include "Residue.hpp"
#include "Geometry.hpp"
#include "CoutColors.hpp"

Residue::Residue()
{
  atom.clear();
  center.clear();
  skip = false;
  altLoc = false;
  corrected = false;
}

Residue::~Residue()
{
  atom.clear();
  center.clear();
  skip = false;
  altLoc = false;
  corrected = false;
}

// All combinations of centers are calculated for every possible
// combination of locations.  For instance, if we have something like
// the following:
// ATOM    816  N   GLU    99      30.629  13.024  26.769 ...
// ATOM    817  CA  GLU    99      29.580  12.059  27.155 ...
// ATOM    818  C   GLU    99      30.037  10.639  26.931 ...
// ATOM    819  O   GLU    99      29.810   9.732  27.750 ...
// ATOM    820  CB  GLU    99      28.301  12.393  26.389 ...
// ATOM    821  CG  GLU    99      27.763  13.769  26.785 ...
// ATOM    822  CD AGLU    99      26.791  14.298  25.745 ...
// ATOM    823  CD BGLU    99      27.482  13.867  28.279 ...
// ATOM    824  OE1AGLU    99      26.375  13.533  24.848 ...
// ATOM    825  OE1BGLU    99      26.910  12.902  28.837 ...
// ATOM    826  OE2BGLU    99      27.834  14.887  28.915 ...
// ATOM    827  OE2AGLU    99      26.437  15.492  25.828 ...
// There would be a total of 8 combinations of possible centers:
// CG, CD A, OE1 A, OE2 A
// CG, CD A, OE1 A, OE2 B
// CG, CD A, OE1 B, OE2 A
// CG, CD A, OE1 B, OE2 B
// CG, CD B, OE1 A, OE2 A
// CG, CD B, OE1 A, OE2 B
// CG, CD B, OE1 B, OE2 A
// CG, CD B, OE1 B, OE2 B
// Thus, we store all of these possible centers in a vector to 
// examined later when we calculate the distances.
//////////////////////////////////////////////////////////////////////

// Calculate the centers for PHE or TYR
// these AA need the following atoms:
//   CG, CZ, CD1, CD2, CE1, CE2
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void Residue::centerPHEorTYR()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 6 atoms, the size of the outer vector
  // is 6 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; CZ = 1 ; CD1 = 2 ; CD2 = 3 ; CE1 = 4 ; CE2 = 5
  vector< vector <Atom*> > temp;
  temp.resize(6);

  // Push all of the important atoms in their respective vectors
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ' )
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " CZ ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD2")
        {
          temp[3].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE1")
        {
          temp[4].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE2")
        {
          temp[5].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }
  
  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 ||
     temp[4].size() == 0 || temp[5].size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the PHE ring at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize(temp[0].size() * temp[1].size() * temp[2].size() *
                temp[3].size() * temp[4].size() * temp[5].size());

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CZ
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the CD1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++) // go through all the CD2
                {
                  for(unsigned int m = 0; m < temp[4].size(); m++) // go through all the CE1
                    {
                      for(unsigned int n = 0; n < temp[5].size(); n++, index++) // go through all the CE2
                        {
                          center[index].plane_info.resize(6);
                          center[index].plane_info[ CG_PLANE_COORD_PTT] = &(temp[0][i]->coord);
                          center[index].plane_info[CD1_PLANE_COORD_PTT] = &(temp[2][k]->coord);
                          center[index].plane_info[CD2_PLANE_COORD_PTT] = &(temp[3][l]->coord);
                          center[index].plane_info[3] = &(temp[1][i]->coord);
                          center[index].plane_info[4] = &(temp[4][m]->coord);
                          center[index].plane_info[5] = &(temp[5][n]->coord);

                          // This is just an average of the all the coordinates
                          center[index] = ( temp[0][i]->coord +
                                            temp[1][j]->coord +
                                            temp[2][k]->coord +
                                            temp[3][l]->coord +
                                            temp[4][m]->coord +
                                            temp[5][n]->coord )/(6);

                          // And this is just so we know what combination of
                          // alternate locations (may be able to take out later)
                          string t;
                          t.insert(t.end(),1,temp[0][i]->altLoc);
                          t.insert(t.end(),1,temp[1][j]->altLoc);
                          t.insert(t.end(),1,temp[2][k]->altLoc);
                          t.insert(t.end(),1,temp[3][l]->altLoc);
                          t.insert(t.end(),1,temp[4][m]->altLoc);
                          t.insert(t.end(),1,temp[5][n]->altLoc);
                          center[index].altLoc = t;
#ifdef DEBUG
                          cout << "CHECK: " << residue << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                        }
                    }
                }
            }
        }
    }
}

void Residue::centerPHEorTYR_altloc()
{
  int numaltlocs = altlocs.size();
  center.resize(numaltlocs);
  if(numaltlocs > 1) this->altLoc = true;
  for(unsigned int al = 0; al < numaltlocs; al++)
    {
      int C_count=0;
      center[al].skip = false;
      center[al].plane_info.resize(6);
      center[al].set(0,0,0);
      // Push all of the important atoms in their respective vectors
      for(unsigned int i =0; i< altlocs[al].size(); i++)
        {
          if(altlocs[al][i]->name == " CG ")
            {
              center[al].plane_info[ CG_PLANE_COORD_PTT] = &(altlocs[al][i]->coord);
              center[al] += altlocs[al][i]->coord;
              C_count++;
            }
          else if(altlocs[al][i]->name == " CZ ")
            {
              center[al].plane_info[3] = &(altlocs[al][i]->coord);
              center[al] += altlocs[al][i]->coord;
              C_count++;
            }
          else if(altlocs[al][i]->name == " CD1")
            {
              center[al].plane_info[CD1_PLANE_COORD_PTT] = &(altlocs[al][i]->coord);
              center[al] += altlocs[al][i]->coord;
              C_count++;
            }
          else if(altlocs[al][i]->name == " CD2")
            {
              center[al].plane_info[CD2_PLANE_COORD_PTT] = &(altlocs[al][i]->coord);
              center[al] += altlocs[al][i]->coord;
              C_count++;
            }
          else if(altlocs[al][i]->name == " CE1")
            {
              center[al].plane_info[4] = &(altlocs[al][i]->coord);
              center[al] += altlocs[al][i]->coord;
              C_count++;
            }
          else if(altlocs[al][i]->name == " CE2")
            {
              center[al].plane_info[5] = &(altlocs[al][i]->coord);
              center[al] += altlocs[al][i]->coord;
              C_count++;
            }
          else
            {
              // This means that this atom is not useful so we flag it 
              altlocs[al][i]->skip = true;
            }
        }
  
      // Error check.  If we don't have at least one of each
      // of the atoms, we throw a warning, set an ignore flag
      // for future reference, and go to next altLoc
      if(C_count != 6)
        {
          skip = true;
#ifndef DISABLE_WARNING
          cout << cyan << "WARNING" << reset << ": Could not find all atoms in the PHE ring at "
               << altlocs[al][0]->resSeq << " altLoc: " << altlocs[al][0]->altLoc << endl;
#endif
        }
      center[al] /= 6;
    }
}

// Calculate the centers for TRP
// these AA need the following atoms:
//   CG, CH2, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void Residue::centerTRP()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 9 atoms, the size of the outer vector
  // is 9 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0  ; CH1 = 1 ; CD1 = 2 ; CD2 = 3 ; NE1 = 4 ; CE2 = 5 ;
  // CE3 = 6 ; CZ2 = 7 ; CZ3 = 8
  vector< vector <Atom*> > temp;
  temp.resize(9);

  // Push all of the important atoms in their respective vectors
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ' )
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " CH2")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD2")
        {
          temp[3].push_back(atom[i]);
        }
      else if(atom[i]->name == " NE1")
        {
          temp[4].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE2")
        {
          temp[5].push_back(atom[i]);
        }
      else if(atom[i]->name == " CE3")
        {
          temp[6].push_back(atom[i]);
        }
      else if(atom[i]->name == " CZ2")
        {
          temp[7].push_back(atom[i]);
        }
      else if(atom[i]->name == " CZ3")
        {
          temp[8].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 ||
     temp[4].size() == 0 || temp[5].size() == 0 ||
     temp[6].size() == 0 || temp[7].size() == 0 ||
     temp[8].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the TRP ring at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }
  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size()*
                 temp[4].size() * temp[5].size() * temp[6].size() * temp[7].size() * temp[8].size());

  // And go through all combinations of the alternate locations
  // (Sorry that it is really freaking ugly)
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CH2
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the CD1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++) // go through all the CD2
                {
                  for(unsigned int m = 0; m < temp[4].size(); m++) // go through all the NE1
                    {
                      for(unsigned int n = 0; n < temp[5].size(); n++) // go through all the CE2
                        {
                          for(unsigned int o = 0; o < temp[6].size(); o++) // go through all the CE3
                            {
                              for(unsigned int p = 0; p < temp[7].size(); p++) // go through all the CZ2
                                {
                                  for(unsigned int q = 0; q < temp[8].size(); q++, index++) // go through all the CZ2
                                    {
                                      center[index].plane_info.resize(3);
                                      center[index].plane_info[CG_PLANE_COORD_PTT]  = &temp[0][i]->coord;
                                      center[index].plane_info[CD1_PLANE_COORD_PTT] = &temp[2][k]->coord;
                                      center[index].plane_info[CD2_PLANE_COORD_PTT] = &temp[3][l]->coord;

                                      // This is a weighted average of the atoms weighted
                                      // by their mass
                                      center[index] = ( temp[0][i]->coord * MASS_C +
                                                        temp[1][j]->coord * MASS_C +
                                                        temp[2][k]->coord * MASS_C +
                                                        temp[3][l]->coord * MASS_C +
                                                        temp[4][m]->coord * MASS_N +
                                                        temp[5][n]->coord * MASS_C +
                                                        temp[6][o]->coord * MASS_C +
                                                        temp[7][p]->coord * MASS_C +
                                                        temp[8][q]->coord * MASS_C )/(MASS_C*8+MASS_N);

                                      // And this is just so we know what combination of
                                      // alternate locations (may be able to take out later)
                                      string t;
                                      t.insert(t.end(),1,temp[0][i]->altLoc);
                                      t.insert(t.end(),1,temp[1][j]->altLoc);
                                      t.insert(t.end(),1,temp[2][k]->altLoc);
                                      t.insert(t.end(),1,temp[3][l]->altLoc);
                                      t.insert(t.end(),1,temp[4][m]->altLoc);
                                      t.insert(t.end(),1,temp[5][n]->altLoc);
                                      t.insert(t.end(),1,temp[6][o]->altLoc);
                                      t.insert(t.end(),1,temp[7][p]->altLoc);
                                      t.insert(t.end(),1,temp[8][q]->altLoc);
                                      center[index].altLoc = t;
#ifdef DEBUG
                                      cout << "CHECK: TRP" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// DO NOT USE THIS!!!! If you use this, a fairy will lose its wings!
// This is only here for legacy reasons because this is what the old
// STAAR code used.  It is wrong.
// Calculate the centers of charge for ASP
// these AA need the following atoms:
//   CG, CB, OD1, OD2
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void Residue::centerASP()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 4 atoms, the size of the outer vector
  // is 4 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; CB = 1 ; OD1 = 2 ; OD2 = 3 ;
  vector< vector <Atom*> > temp;
  temp.resize(4);

  // Push all of the important atoms in their respective vectors
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " H  ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD2")
        {
          temp[3].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the ASP side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize(temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size());

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CB
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the OD1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++, index++) // go through all the OD2
                {
                  center[index].plane_info.resize(3);
                  center[index].plane_info[C__PLANE_COORD_AG]  = &temp[0][i]->coord;
                  center[index].plane_info[O_1_PLANE_COORD_AG] = &temp[2][k]->coord;
                  center[index].plane_info[O_2_PLANE_COORD_AG] = &temp[3][l]->coord;

                  // This is just an average of the all the coordinates
                  center[index] = ( temp[0][i]->coord * CHARGE_C +
                                    temp[1][j]->coord * CHARGE_H + // treating CB as H in formic acid
                                    temp[2][k]->coord * CHARGE_O +
                                    temp[3][l]->coord * CHARGE_O ) * -1;
                  
                  // And this is just so we know what combination of
                  // alternate locations (may be able to take out later)
                  string t;
                  t.insert(t.end(),1,temp[0][i]->altLoc);
                  t.insert(t.end(),1,temp[1][j]->altLoc);
                  t.insert(t.end(),1,temp[2][k]->altLoc);
                  t.insert(t.end(),1,temp[3][l]->altLoc);
                  center[index].altLoc = t;
#ifdef DEBUG
                  cout << "CHECK: ASP" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                }
            }
        }
    }
}

// DO NOT USE THIS!!!! If you use this, puppies will be harmed!
// This is only here for legacy reasons because this is what the old
// STAAR code used.  It is wrong.  
// Calculate the centers of charge for GLU
// these AA need the following atoms:
//   CG, CD, OE1, OE2
// all of the other ones don't affect the center
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void Residue::centerGLU()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 4 atoms, the size of the outer vector
  // is 4 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; CD = 1 ; OE1 = 2 ; OE2 = 3
  vector< vector <Atom*> > temp;
  temp.resize(4);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " H  ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " CD ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE1")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE2")
        {
          temp[3].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the GLU side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize(temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size());

    // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int i = 0; i < temp[0].size(); i++) // go through all the CG
    {
      for(unsigned int j = 0; j < temp[1].size(); j++) // go through all the CD
        {
          for(unsigned int k = 0; k < temp[2].size(); k++) // go through all the OE1
            {
              for(unsigned int l = 0; l < temp[3].size(); l++, index++) // go through all the OE2
                {

                  center[index].plane_info.resize(3);
                  center[index].plane_info[C__PLANE_COORD_AG]  = &(temp[1][j]->coord);
                  center[index].plane_info[O_1_PLANE_COORD_AG] = &(temp[2][k]->coord);
                  center[index].plane_info[O_2_PLANE_COORD_AG] = &(temp[3][l]->coord);
                  
                  // This is a weighted average weighted by the charges.
                  // Since the charges sum to -1, we just need to negate the signs
                  center[index] = ( temp[0][i]->coord * CHARGE_H + // CG treated as H in formic acid
                                    temp[1][j]->coord * CHARGE_C +
                                    temp[2][k]->coord * CHARGE_O +
                                    temp[3][l]->coord * CHARGE_O ) * -1;
                  
                  // And this is just so we know what combination of
                  // alternate locations (may be able to take out later)
                  string t;
                  t.insert(t.end(),1,temp[0][i]->altLoc);
                  t.insert(t.end(),1,temp[1][j]->altLoc);
                  t.insert(t.end(),1,temp[2][k]->altLoc);
                  t.insert(t.end(),1,temp[3][l]->altLoc);
                  center[index].altLoc = t;
#ifdef DEBUG
                  cout << "CHECK: GLU" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                }
            }
        }
    }
}

// Calculate the centers for ASP
// these AA need the following atoms:
//   OD1 and OD2
// all of the other ones don't affect the center
// This is is really simple: all it does is add the
// OD1 and OD2 coords to the center vector.
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void Residue::centerASP_oxygen()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 3 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CG = 0 ; OD1 = 1 ; OD2 = 2
  vector< vector <Atom*> > temp;
  temp.resize(3);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " CG ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD1")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OD2")
        {
          temp[2].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if( temp[1].size() == 0 ||
      temp[2].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the ASP side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }
  unsigned int offset = temp[0].size() * temp[1].size() * temp[2].size();
  // Allocate size for all of the possible centers
  center.resize(offset * 2);

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the OD1
    {      
      for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the CD
        {
          for(unsigned int l = 0; l < temp[2].size(); l++, index++) // go through all the OD2
            {
              
              center[index].plane_info.resize(3);
              center[index].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
              center[index + offset].plane_info.resize(3);
              center[index + offset].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index + offset].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index + offset].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
                  
              center[index] = temp[1][k]->coord;
              center[index + offset] = temp[2][l]->coord;
                  
              // And this is just so we know what combination of
              // alternate locations (may be able to take out later)
              string t;
              t.insert(t.end(),1,temp[0][j]->altLoc);
              t.insert(t.end(),1,temp[1][k]->altLoc);
              t.insert(t.end(),1,temp[2][l]->altLoc);
              center[index].altLoc = t;
              center[index + offset].altLoc = t;
#ifdef DEBUG
              cout << "CHECK: ASP" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
            }
        }
    }
}

void Residue::centerASP_oxygen_altloc()
{

  int numaltlocs = altlocs.size();
  center.resize(numaltlocs*2);
  if(numaltlocs > 1) this->altLoc = true;
  // Push all of the important atoms in their respective vectors  
  for(unsigned int al = 0; al < numaltlocs; al++)
    {
      int atom_count = 0;
      center[al].skip = false;
      center[al].plane_info.resize(3);
      center[al].set(0,0,0);
      center[al+numaltlocs].skip = false;
      center[al+numaltlocs].plane_info.resize(3);
      center[al+numaltlocs].set(0,0,0);
      for(unsigned int i =0; i< altlocs[al].size(); i++)
        {
          if(altlocs[al][i]->name == " CG ")
            {
              center[al].plane_info[C__PLANE_COORD_AG]  = &(altlocs[al][i]->coord);
              center[al+numaltlocs].plane_info[C__PLANE_COORD_AG]  = &(altlocs[al][i]->coord);
              atom_count++;
            }
          else if(altlocs[al][i]->name == " OD1")
            {
              center[al].plane_info[O_1_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al+numaltlocs].plane_info[O_1_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al] = altlocs[al][i]->coord;
              atom_count++;
            }
          else if(altlocs[al][i]->name == " OD2")
            {
              center[al].plane_info[O_2_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al+numaltlocs].plane_info[O_2_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al+numaltlocs] = altlocs[al][i]->coord;
              atom_count++;
            }
          else
            {
              // This means that this atom is not useful so we flag it 
              altlocs[al][i]->skip = true;
            }
        }

      // Error check.  If we don't have at least one of each
      // of the atoms, we throw a warning, set an ignore flag
      // for future reference, and leave the function
      if(atom_count != 3)
        {
          skip = true;
#ifndef DISABLE_WARNING
          cout << cyan << "WARNING" << reset << ": Could not find all atoms in the ASP side chain at "
               << atom[0]->resSeq << " altloc: " << altlocs[al][0]->altLoc << endl;
#endif
        }
    }
}

// Calculate the centers for GLU
// these AA need the following atoms:
//   OE1 and OE2
// all of the other ones don't affect the center
// This is is really simple: all it does is add the
// OE1 and OE2 coords to the center vector.
// This also sets the plane coordinates used to calculate 
// the angle later in the program
void Residue::centerGLU_oxygen()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 3 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // CD = 0 ; OE1 = 1 ; OE2 = 2
  vector< vector <Atom*> > temp;
  temp.resize(3);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " CD ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE1")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " OE2")
        {
          temp[2].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0)
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the GLU side chain at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }
  unsigned int offset = temp[0].size() * temp[1].size() * temp[2].size();
  // Allocate size for all of the possible centers
  center.resize(offset * 2);

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the OE1
    {
      for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the CD
        {
          for(unsigned int l = 0; l < temp[2].size(); l++, index++) // go through all the OE2
            {
              
              center[index].plane_info.resize(3);
              center[index].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
              center[index + offset].plane_info.resize(3);
              center[index + offset].plane_info[C__PLANE_COORD_AG]  = &(temp[0][j]->coord);
              center[index + offset].plane_info[O_1_PLANE_COORD_AG] = &(temp[1][k]->coord);
              center[index + offset].plane_info[O_2_PLANE_COORD_AG] = &(temp[2][l]->coord);
                  
              center[index] = temp[1][k]->coord;
              center[index + offset] = temp[2][l]->coord;
                  
              // And this is just so we know what combination of
              // alternate locations (may be able to take out later)
              string t;
              t.insert(t.end(),1,temp[0][j]->altLoc);
              t.insert(t.end(),1,temp[1][k]->altLoc);
              t.insert(t.end(),1,temp[2][l]->altLoc);
              center[index].altLoc = t;
              center[index + offset].altLoc = t;
#ifdef DEBUG
              cout << "CHECK: GLU" << atom[0]->resSeq << " : " << center[index] << endl;
#endif
            }
        }
    }
}

void Residue::centerGLU_oxygen_altloc()
{

  int numaltlocs = altlocs.size();
  center.resize(numaltlocs*2);
  if(numaltlocs > 1) this->altLoc = true;
  // Push all of the important atoms in their respective vectors  
  for(unsigned int al = 0; al < numaltlocs; al++)
    {
      int atom_count = 0;
      center[al].skip = false;
      center[al].plane_info.resize(3);
      center[al].set(0,0,0);
      center[al+numaltlocs].skip = false;
      center[al+numaltlocs].plane_info.resize(3);
      center[al+numaltlocs].set(0,0,0);
      for(unsigned int i =0; i< altlocs[al].size(); i++)
        {
          if(altlocs[al][i]->name == " CD ")
            {
              center[al].plane_info[C__PLANE_COORD_AG]  = &(altlocs[al][i]->coord);
              center[al+numaltlocs].plane_info[C__PLANE_COORD_AG]  = &(altlocs[al][i]->coord);
              atom_count++;
            }
          else if(altlocs[al][i]->name == " OE1")
            {
              center[al].plane_info[O_1_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al+numaltlocs].plane_info[O_1_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al] = altlocs[al][i]->coord;
              atom_count++;
            }
          else if(altlocs[al][i]->name == " OE2")
            {
              center[al].plane_info[O_2_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al+numaltlocs].plane_info[O_2_PLANE_COORD_AG] = &(altlocs[al][i]->coord);
              center[al+numaltlocs] = altlocs[al][i]->coord;
              atom_count++;
            }
          else
            {
              // This means that this atom is not useful so we flag it 
              altlocs[al][i]->skip = true;
            }
        }

      // Error check.  If we don't have at least one of each
      // of the atoms, we throw a warning, set an ignore flag
      // for future reference, and leave the function
      if(atom_count != 3)
        {
          skip = true;
#ifndef DISABLE_WARNING
          cout << cyan << "WARNING" << reset << ": Could not find all atoms in the GLU side chain at "
               << atom[0]->resSeq << " altloc: " << altlocs[al][0]->altLoc << endl;
#endif
        }
    }
}

void Residue::centerPHEorTYR_simplified()
{
  center.resize(1);

  // Temp variable that will have the center coordinates
  Coordinates tempCenter(0.0, 0.0, 0.0);

  // Here, we are just getting 2 C atoms. Doesn't matter
  // which ones we choose as long as they are
  // diametrically opposed. Based off of 
  // http://comp.chem.nottingham.ac.uk/dichrocalc/files/atomlabels-sidechains.png
  // the following are diametically opposed:
  //   CG  - CZ
  //   CD2 - CE1
  //   CD1 - CE2
  center[0].set(0.0, 0.0, 0.0);
  center[0].plane_info.resize(3);
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->name == " CE2")
        {
          center[0] += atom[i]->coord;
          center[0].plane_info[CE2_PLANE_COORD_PTT] = &atom[i]->coord;
        }
      else if(atom[i]->name == " CD1")
        {
          center[0] += atom[i]->coord;
          center[0].plane_info[CD1_PLANE_COORD_PTT] = &atom[i]->coord;
        }
      else if(atom[i]->name == " CG ")
        {
          center[0].plane_info[CG_PLANE_COORD_PTT] = &atom[i]->coord;
        }
    }
  center[0] /= 2;
}

void Residue::centerASP_charge()
{
  center.resize(1);
  Coordinates tempCenter(0.0, 0.0, 0.0);
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->name == " CG ")
        {
          center[0] = atom[i]->coord;
          tempCenter += atom[i]->coord;
        }
      else if(atom[i]->name == " H  ")
        {
          tempCenter -= atom[i]->coord;
        }
    }
  center[0] += (tempCenter * HYDROGEN_BOND_DISTANCE) / tempCenter.norm();
}

void Residue::centerGLU_charge()
{
  center.resize(1);
  Coordinates tempCenter(0.0, 0.0, 0.0);
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->name == " CD ")
        {
          center[0] = atom[i]->coord;
          tempCenter += atom[i]->coord;
        }
      else if(atom[i]->name == " H  ")
        {
          tempCenter -= atom[i]->coord;
        }
    }
  center[0] += (tempCenter * HYDROGEN_BOND_DISTANCE) / tempCenter.norm();
}


// Calculate the centers for PO4, 2HP, PI
// these AA need the following atoms:
//   P, O1, O2, O3, and O4
void Residue::centerPO4or2HPorPI()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 5 atoms, the size of the outer vector
  // is 5 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1 = 1 ; O2 = 2 ; O3 = 3 ; O4 = 4
  vector< vector <Atom*> > temp;
  temp.resize(5);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " O1 ")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 ")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 ")
        {
          temp[3].push_back(atom[i]);
        }
      else if(atom[i]->name == " O4 ")
        {
          temp[4].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 ||
     temp[4].size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << " chain " << atom[0]->chainID << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size() * temp[4].size() );

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  float totalmass = MASS_P + MASS_O * 4;
  for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the P
    {
      for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the O1
        {
          for(unsigned int l = 0; l < temp[2].size(); l++) // go through all the O2
            {
              for(unsigned int m = 0; m < temp[3].size(); m++) // go through all the O3
                {
                  for(unsigned int n = 0; n < temp[4].size(); n++, index++) // go through all the O4
                    {
                      center[index] = temp[0][j]->coord * MASS_P +
                        ( temp[1][k]->coord + temp[2][l]->coord +
                          temp[3][m]->coord + temp[4][n]->coord ) * MASS_O;

                      center[index] /= totalmass;

                      center[index].plane_info.resize(5);
                      center[index].plane_info[0] = &(temp[0][j]->coord);
                      center[index].plane_info[1] = &(temp[1][k]->coord);
                      center[index].plane_info[2] = &(temp[2][l]->coord);
                      center[index].plane_info[3] = &(temp[3][m]->coord);
                      center[index].plane_info[4] = &(temp[4][n]->coord);

#ifdef DEBUG
                      cout << "CHECK: " << residue  << " " << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                    }
                }
            }
        }
    }
}

void Residue::centerPO4or2HPorPI_charge()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 6 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1 = 1 ; O2 = 2 ; O3 = 3 ; O4 = 4
  Atom* P = NULL;
  vector< Atom* > O;
  vector< Atom* > H;

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          P = atom[i];
        }
      else if(atom[i]->name == " O1 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O4 ")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " H  ")
        {
          H.push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check to see if we have 1 P, 4 O, and at least 1 H
  if(P == NULL || O.size() != 4 || H.size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << " chain " << atom[0]->chainID << endl;
#endif
      cout << cyan << "Skipping " << reset << " because there weren't enough hydrogens in " 
           << residue << " " << atom[0]->resSeq << " chain " << atom[0]->chainID << endl;
      return;
    }

  center.resize( 1 );
  // Here is an estimate of the center by doing a charge weighted average
  // This could be wrong, but it is just something to put in here for the
  // time being until I talk with Dr. Hinde.
  center[0] = P->coord * CHARGE_P;
  Coordinates temp(0.0, 0.0, 0.0);
  temp  = O[0]->coord;
  temp += O[1]->coord;
  temp += O[2]->coord;
  temp += O[3]->coord;
  center[0] += temp*CHARGE_O;
  temp.set(0.0,0.0,0.0);
  for(unsigned int i=0; i<H.size(); i++)
    {
      temp += H[i]->coord;
    }
  center[0] += temp * CHARGE_H;

  // And since the formal charge of PO4 is -3, divide by it
  if(residue == "PO4")
    center[0] /= -3;
  else if(residue == "2HP")
    center[0] *= -1;
  else if(residue == " PI")
    {
      center[0] /= -2;
    }
}

// Calculate the centers for 2PO and PO3
// these AA need the following atoms:
//   P, O1(P), O2(P), O3(P)
void Residue::center2POorPO3()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 5 atoms, the size of the outer vector
  // is 5 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1(P) = 1 ; O2(P) = 2 ; O3(P) = 3 
  vector< vector <Atom*> > temp;
  temp.resize(4);

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          temp[0].push_back(atom[i]);
        }
      else if(atom[i]->name == " O1 " || atom[i]->name == " O1P")
        {
          temp[1].push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 " || atom[i]->name == " O2P")
        {
          temp[2].push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 " || atom[i]->name == " O3P")
        {
          temp[3].push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check.  If we don't have at least one of each
  // of the atoms, we throw a warning, set an ignore flag
  // for future reference, and leave the function
  if(temp[0].size() == 0 || temp[1].size() == 0 ||
     temp[2].size() == 0 || temp[3].size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  // Allocate size for all of the possible centers
  center.resize( temp[0].size() * temp[1].size() * temp[2].size() * temp[3].size() );

  // And go through all combinations of the alternate locations
  unsigned int index = 0;
  float totalmass = MASS_P + MASS_O * 4;
  for(unsigned int j = 0; j < temp[0].size(); j++) // go through all the P
    {
      for(unsigned int k = 0; k < temp[1].size(); k++) // go through all the O1
        {
          for(unsigned int l = 0; l < temp[2].size(); l++) // go through all the O2
            {
              for(unsigned int m = 0; m < temp[3].size(); m++) // go through all the O3
                {
                  center[index] = temp[0][j]->coord * MASS_P +
                    ( temp[1][k]->coord + temp[2][l]->coord +
                      temp[3][m]->coord  ) * MASS_O;

                  center[index] /= totalmass;

                  center[index].plane_info.resize(4);
                  center[index].plane_info[0] = &(temp[0][j]->coord);
                  center[index].plane_info[1] = &(temp[1][k]->coord);
                  center[index].plane_info[2] = &(temp[2][l]->coord);
                  center[index].plane_info[3] = &(temp[3][m]->coord);

#ifdef DEBUG
                  cout << "CHECK: " << residue  << " " << atom[0]->resSeq << " : " << center[index] << endl;
#endif
                }
            }
        }
    }
}

void Residue::center2POorPO3_charge()
{
  // This makes a vector of a vector of ATOM* 
  // since we are looking for 6 atoms, the size of the outer vector
  // is 3 while the size of the inner vectors are dependent on the 
  // number of alternate locations for each one. Indexes are:
  // P = 0 ; O1 = 1 ; O2 = 2 ; O3 = 3 ; O4 = 4
  Atom* P = NULL;
  vector< Atom* > O;
  vector< Atom* > H;

  // Push all of the important atoms in their respective vectors  
  for(unsigned int i =0; i< atom.size(); i++)
    {
      if(atom[i]->altLoc != ' ')
        {
          altLoc = true;
        }
      if(atom[i]->name == " P  ")
        {
          P = atom[i];
        }
      else if(atom[i]->name == " O1 " || atom[i]->name == " O1P")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O2 " || atom[i]->name == " O2P")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " O3 " || atom[i]->name == " O3P")
        {
          O.push_back(atom[i]);
        }
      else if(atom[i]->name == " H  ")
        {
          H.push_back(atom[i]);
        }
      else
        {
          // This means that this atom is not useful so we flag it 
          atom[i]->skip = true;
        }
    }

  // Error check to see if we have 1 P, 4 O, and at least 1 H
  if(P == NULL || O.size() != 3 || H.size() == 0 )
    {
      skip = true;
#ifndef DISABLE_WARNING
      cout << cyan << "WARNING" << reset << ": Could not find all atoms in the " << residue <<  " at "
           << atom[0]->resSeq << endl;
#endif
      return;
    }

  center.resize( 1 );
  // Here is an estimate of the center by doing a charge weighted average
  // This could be wrong, but it is just something to put in here for the
  // time being until I talk with Dr. Hinde.
  center[0] = P->coord * CHARGE_P;
  Coordinates temp(0.0, 0.0, 0.0);
  temp  = O[0]->coord;
  temp += O[1]->coord;
  temp += O[2]->coord;
  center[0] += temp*CHARGE_O;
  temp.set(0.0,0.0,0.0);
  for(unsigned int i=0; i<H.size(); i++)
    {
      temp += H[i]->coord;
    }
  center[0] += temp * CHARGE_H;

  // And since the formal charge of PO4 is -3, divide by it
  if(residue == "2PO")
    center[0] /= -2;
  else if(residue == "PO3")
    center[0] /= -3;
}

template<class T>
int find_index(vector<T>& vec, T to_find)
{
  for(int i=0; i<vec.size(); i++)
    {
      if(vec[i] == to_find) return i;
    }
  return vec.size();
}

void Residue::determineAltLoc(vector<char>&altloc_ids)
{
  int count = altloc_ids.size();
  if(count == 0)
    {
      count = 1;
      altlocs.resize(1);
      for(int i=0; i<this->atom.size(); i++)
        {
          this->atom[i]->skip = false;
          altlocs[0].push_back(this->atom[i]);
        }
    }
  else
    {
      altlocs.resize(count);

      for(int i=0; i<this->atom.size(); i++)
        {
          if(this->atom[i]->altLoc == ' ')
            {
              for(int j=0; j<count; j++)
                {
                  this->atom[i]->skip = false;
                  altlocs[j].push_back(this->atom[i]);
                }
            }
          else
            {
              this->atom[i]->skip = false;
              int idx = find_index(altloc_ids, this->atom[i]->altLoc);
              altlocs[idx].push_back(this->atom[i]);
            }
        }
    }

  // for(int j=0; j<count; j++)
  //   {
  //     for(int i=0; i<altlocs[j].size(); i++)
  //       cout << *altlocs[j][i] << endl;
  //   }

}

// Calculates the center of the amino acid
void Residue::calculateCenter(bool centerOfCharge)
{
//cout << "calculateCenter: residue=" << residue << ", centerOfCharge=" << centerOfCharge << ", altlocs.size()=" << altlocs.size() << endl;

  if( !centerOfCharge )
    {
      // TRP
      if ( residue == "TRP" )
        {
          centerTRP();
        }
      // PHE or TYR
      else if ( residue == "PHE" || residue == "TYR" )
        {
          centerPHEorTYR_altloc();
        }
      // ASP
      else if (residue == "ASP")
        {
          centerASP_oxygen_altloc();
        }
      // GLU
      else if (residue == "GLU")
        {
          centerGLU_oxygen_altloc();
        }
      else if (residue == "PO4" || residue == "2HP" || residue == " PI")
        {
          centerPO4or2HPorPI();
        }
      else if (residue == "2PO" || residue == "PO3")
        {
          center2POorPO3();
        }
      // This is for all of other amino acids out there that we 
      // don't support
      else
        {
          //cerr << red << "ERROR" << reset << ": " << residue << " is not yet supported." << endl;
          skip = true;
        }
    }
  // Now we are going to calculate the center of charges for the formate
  else
    {
      // PHE or TYR
      if ( residue == "PHE" || residue == "TYR" )
        {
          // So, this does a simplified center of mass calculation that
          // Dr. Hinde used.  See function notes above for details.
          centerPHEorTYR_simplified();
        }
      // ASP
      else if (residue == "ASP")
        {
          centerASP_charge();
        }
       // GLU
      else if (residue == "GLU")
        {
          centerGLU_charge();
        }
      else if (residue == "PO4" || residue == "2HP" || residue == " PI")
        {
          centerPO4or2HPorPI_charge();
        }
      else if (residue == "2PO" || residue == "PO3")
        {
          center2POorPO3();
        }
    }
}


/*
  Given a vector containing a carbon ring (6 carbon atoms in a hexagon
  configuration), find all atoms that are at most 2 bonds away from
  a ring atom. Return these atoms as a vector.
*/
vector<Atom*> Residue::findAdditionalAtoms(vector<Atom*> ring)
{
  vector<Atom*> temp;

  char ring_alt_loc = ring[0]->altLoc;

  for(int i=0; i<ring.size(); i++)
    {
      Atom *ring_carbon = ring[i];
      for(int j=0; j<atom.size(); j++)
        {
          Atom *other_atom= atom[j];

          if (find(ring.begin(), ring.end(), other_atom) != ring.end() ||
              other_atom->altLoc != ring_alt_loc)
            {
              // This atom is part of the ring
              continue;
            }

          if (ring_carbon->isBonded(*other_atom))
            {
              temp.push_back(other_atom);
            }
        }
    }

  int n_found_atoms = temp.size();
  for(int i=0; i<n_found_atoms; i++)
    {
      Atom *one_bond_atom = temp[i];
      for (int j=0; j<atom.size(); j++)
        {
          Atom *other_atom = atom[j];

          if (find(ring.begin(), ring.end(), other_atom) != ring.end() ||
              find(temp.begin(), temp.end(), other_atom) != temp.end() ||
              other_atom->altLoc != ring_alt_loc)
            {
              // This atom is part of the ring or one bond away from the ring
              continue;
            }
          if (one_bond_atom->isBonded(*other_atom))
            {
              temp.push_back(other_atom);
            }
        }
    }

  return temp;
}

bool Residue::findCarbonRings()
{
  bool foundRings = false;
  const float TARGET_DIST = 2.8;

  vector< vector<Atom*> > pairs;
  vector<Coordinates> midpoints;
  //cout << residue << endl;

  // Loop through all pairs of atoms, only considering carbon atoms.
  // If two carbon atoms are the correct distance apart to be directly
  // across from each other in a benzene ring, then store that pair
  // and their midpoint.
  for(int i=0; i<atom.size(); i++)
    {
      Atom * firstCarbon = atom[i];
      if (firstCarbon->element != " C")
        {
          continue;
        }

      for(int j=i+1; j<atom.size(); j++)
        {
          Atom * secondCarbon = atom[j];
          if (secondCarbon->element != " C")
            {
              // We're only looking for carbon-only rings.
              continue;
            }
          if (firstCarbon->altLoc != secondCarbon->altLoc)
            {
              // Don't compare atoms that aren't part of the same
              // alternate location set.
              continue;
            }
          float dist = firstCarbon->coord.distance(secondCarbon->coord);
/*
cout << "Checking distance between" << endl;
cout << "  " << firstCarbon->line << endl;
cout << "  " << secondCarbon->line << endl;
cout << "-- " << dist << endl << endl;
*/
          //cout << "Distance: " << dist << endl;

          
          if (dist >= TARGET_DIST - 0.11 && dist <= TARGET_DIST + 0.11)
            {
              vector<Atom*> pair;
              pair.push_back(firstCarbon);
              pair.push_back(secondCarbon);
              pairs.push_back(pair);

              midpoints.push_back(firstCarbon->coord.midpoint(secondCarbon->coord));
            }
        }
/*
cout << endl << "Found the following pairs for " << residue << ":" << endl;
for (int I=0; I < pairs.size(); I++) {
  cout << "  " << pairs[I][0]->line << endl;
  cout << "  " << pairs[I][1]->line << endl;
  cout << "--" << endl;
}
*/

    }

  // Collect pairs whose midpoints are approximately equal
  // and that are part of the same location of the ligand.
  vector< vector<int> > potential_rings;
  for(int i=0; i<midpoints.size(); i++)
    {
      vector<int> index_set;
      index_set.push_back(i);
      for(int j=i+1; j<midpoints.size(); j++)
        {
          float dist = midpoints[i].distance(midpoints[j]);
          if (dist <= 0.3 && pairs[i][0]->altLoc == pairs[j][0]->altLoc)
            {
              index_set.push_back(j);
            }
          //cout << dist << endl;
        }
      potential_rings.push_back(index_set);
    }

  // A set of three pairs means a ring. Store the constituent
  // atoms, find its center, maybe compute it's plane info, and
  // set the appropriate flag.
  for(int i=0; i<potential_rings.size(); i++)
    {
      vector<int> index_set = potential_rings[i];
      if (index_set.size() == 3)
        {
          foundRings = true;
          vector<Atom*> newRing;
          Coordinates ringCenter;
          for(int j=0; j<3; j++)
            {
              ringCenter += midpoints[index_set[j]];
              vector<Atom*> pair = pairs[index_set[j]];
              newRing.push_back(pair[0]);
              newRing.push_back(pair[1]);
            }
          ringCenter /= 3.0;
          ringCenter.skip = false;

          // Get the plane info for this ring. (Not quite sure if
          // this is correct...)
          // In the following structure,
          //                C1 -- C2
          //               /        \
          //              C6        C3
          //               \        /
          //                C5 -- C4
          // we want either (C1, C3, C5) or (C2, C5, C6), I think.
          ringCenter.plane_info.resize(3);
          ringCenter.plane_info[CE2_PLANE_COORD_PTT] = &pairs[index_set[0]][0]->coord;
          ringCenter.plane_info[CD1_PLANE_COORD_PTT] = &pairs[index_set[0]][1]->coord;
          Coordinates* coord1 = &pairs[index_set[1]][0]->coord;
          Coordinates* coord2 = &pairs[index_set[1]][1]->coord;
          float dist1 = pairs[index_set[0]][1]->coord.distance(*coord1);
          if (dist1 >= 1.35 && dist1 <= 1.45)
            {
              ringCenter.plane_info[CG_PLANE_COORD_PTT] = coord1;
            }
          else
            {
              ringCenter.plane_info[CG_PLANE_COORD_PTT] = coord2;
            }
//cout << "RING:" << endl;
//for (int I=0; I<6; I++)
//  cout << newRing[I]->line << endl;

          //carbonRings.push_back(newRing);
          //carbonRingCenters.push_back(ringCenter);
          //center.push_back(ringCenter);

          // Handle alternate locations for carbon rings.
          if (newRing[0]->altLoc == ' ')
            {
              // If the atoms of the ring don't have a specified location,
              // then we can go ahead and put them in their own vectors.
              vector< vector<Atom*> > newRingVector;
              newRingVector.push_back(newRing);
              carbonRings.push_back(newRingVector);

              vector<Coordinates> ringCenterVector;
              ringCenterVector.push_back(ringCenter);
              carbonRingCenters.push_back(ringCenterVector);

              vector<Atom*> newAdditionalAtoms = findAdditionalAtoms(newRing);
              vector< vector<Atom*> > newAdditionalAtomsVector;
              newAdditionalAtomsVector.push_back(newAdditionalAtoms);
              additionalAtoms.push_back(newAdditionalAtomsVector);

/*
if (newRing.size() == 0)
  cout << "What the fucking fuck?" << endl;
cout << "Found a ring:" << endl;
for (int I; I<newRing.size(); I++)
  cout << "\t" << newRing[I]->line << endl;
cout << "  with additional atoms:" << endl;
for (int I; I<newAdditionalAtoms.size(); I++)
  cout << "\t" << newAdditionalAtoms[I]->line << endl;
cout << endl;
*/
            }
          else
            {
              // Otherwise, we need to look through all the other carbon
              // rings to see if there is a matching alternate.
              bool foundAltLoc = false;
              for (int j=0; j<carbonRings.size(); j++)
                {
                  vector< vector<Atom*> > potential_alt_set = carbonRings[j];

                  if (newRing[0]->altLoc == potential_alt_set[0][0]->altLoc)
                    {
                      // If the atoms have the same alternate location
                      // then they must be part of completely separate
                      // rings.
                      continue;
                    }

                  for (int k=0; k<6; k++)
                    {
                      string carbonName = newRing[k]->name;
                      bool foundCarbon = false;
                      for (int ii=0; ii<6; ii++)
                        {
                          if (carbonName == potential_alt_set[0][ii]->name)
                            {
                              foundCarbon = true;	
                              break;
                            }
                        }

                      if (!foundCarbon)
                        {
                          // Abort searching this current ring if we don't find one
                          // of the needed carbons.
                          break;
                        }
                      if (foundCarbon && k == 5)
                        {
                          // We've found all 6 carbons, so this is an alternate!
                          foundAltLoc = true;
                          break;
                        }
                    }

                  if (foundAltLoc)
                   {
                     // We've found a matching set of alternates for this carbon ring,
                     // so store it and move on.
                     carbonRings[j].push_back(newRing);
                     carbonRingCenters[j].push_back(ringCenter);
                     additionalAtoms[j].push_back(findAdditionalAtoms(newRing));
                     break;
                   }
                }
              if (!foundAltLoc)
                {
                  // We haven't found an alternate location for this ring, so make
                  // a new alternate set for it and store it. This shouldn't happen,
                  // but some authors have stored alternate locations as separate
                  // ligands, and we have decided to ignore this.
                  vector< vector<Atom*> > newRingVector;
                  newRingVector.push_back(newRing);
                  carbonRings.push_back(newRingVector);

                  vector<Coordinates> ringCenterVector;
                  ringCenterVector.push_back(ringCenter);
                  carbonRingCenters.push_back(ringCenterVector);

                  vector<Atom*> newAdditionalAtoms = findAdditionalAtoms(newRing);
                  vector< vector<Atom*> > newAdditionalAtomsVector;
                  newAdditionalAtomsVector.push_back(newAdditionalAtoms);
                  additionalAtoms.push_back(newAdditionalAtomsVector);

/*
cout << "Found a ring:" << endl;
for (int I=0; I<newRing.size(); I++)
  cout << "\t" << newRing[I]->line << endl;
cout << "  with additional atoms:" << endl;
for (int I=0; I<newAdditionalAtoms.size(); I++)
  cout << "\t" << newAdditionalAtoms[I]->line << endl;
cout << endl;
*/
                }
           }

        }
    }
  //cout << "Found " << pairs.size() << " pairs" << endl;

  if (foundRings)
    {
      skip = false;
    }
 

  return foundRings;
}

void Residue::calculateAnglesPreHydrogens(Residue aa2,
                                            int index1,
                                            int index2,
                                            float* angle,
                                            float* angle1,
                                            float* angleP)
{
  Residue aa1 = *this;
  Coordinates planeP;
  Coordinates planeProject;

  // Calculate the equation of the plane
  float det = getPlaneEquation( *(aa1.center[index1].plane_info[ CG_PLANE_COORD_PTT]),
                                *(aa1.center[index1].plane_info[CD1_PLANE_COORD_PTT]),
                                *(aa1.center[index1].plane_info[C_2_PLANE_COORD_PTT]),
                                &planeP);

  // Calculate the angle between a plane and a line
  *angle = angleBetweenPlaneAndLine( planeP,
                                     aa1.center[index1],
                                     aa2.center[index2]);
  
  // Calculate the "plane project" coordinates
  planeProjectCoordinate( planeP, 
                          aa2.center[index2],
                          -1 * det, 
                          &planeProject);
  
  // find the second angle betwen the CG ATOM and the "plane project"
  *angle1 = findAngle(*aa1.center[index1].plane_info[ CG_PLANE_COORD_PTT],
                      aa1.center[index1],
                      planeProject);
  
  // finally, calculate the angle between the planes of AA1 and AA2
  *angleP = calculateAngleBetweenPlanes( planeP,
                                         aa2,
                                         index2 );

}

bool Residue::calculateDistancesAndAnglesPostHydrogens(Residue aa2,
                                                         Coordinates closestOxygen,
                                                         float threshold,
                                                         float* dist,
                                                         float* distOxy,
                                                         float* distOxy2,
                                                         float* angle,
                                                         float* angleOxy,
                                                         float* angleOxy2)
{
  Residue aa1 = *this;

//cout << "calculateDistancesAndAnglesPostHydrogens for " << residue << " and " << aa2.residue << endl; 
//cout << "length of aa1.center[0].plane_info: " << aa1.center[0].plane_info.size() << endl;
//cout << aa1.center[0].plane_info[0] << " " << endl;


  // These are the 3 points in the benzene ring determined in centerPHEorTYR_simplified()
  Coordinates dBenzene1 = *aa1.center[0].plane_info[1] - *aa1.center[0].plane_info[0];
  Coordinates dBenzene2 = *aa1.center[0].plane_info[2] - *aa1.center[0].plane_info[0];

  // These values are just to match the Perl script
  float a = dBenzene1.x;
  float b = dBenzene1.y;
  float c = dBenzene1.z;
  float d = dBenzene2.x;
  float e = dBenzene2.y;
  float f = dBenzene2.z;

  // Get the perpendicular vector
  float xp = b * f - c * e;
  float yp = c * d - a * f;
  float zp = a * e - b * d;
  Coordinates perp(xp, yp, zp);  

  // Calculate the distance between the centers
  // This is the vector pointing from the benzene center to formate center of charge
  Coordinates distance = aa2.center[0] - aa1.center[0];

  
  float num = dotProduct(perp, distance);
  float perpnorm = perp.norm();
  float distFromMassToChg = distance.norm();
  float denom = perpnorm * distFromMassToChg;

  if(distFromMassToChg > threshold)
    {
      cout << gray << "Note" << reset << ": post hydrogen distance, " << distFromMassToChg <<" > "<< threshold << "A.  Skipping." << endl;
      return false;
    }

  if(denom == 0)
    {
      cerr << red << "Error" << reset << ": denom is zero.  Skipping residue" << endl;
      return false;
    }  

  // We already have one of the oxygens as an input param
  // Now let's find the other one
  Coordinates otherOxygen;
  for(int i = 0; i < aa2.atom.size(); i++)
    {
      // Only look at the oxygen atoms
      if(aa2.atom[i]->element == "O")
        {
          // Make sure this oxygen is different than the closest one
          if(aa2.atom[i]->coord.x != closestOxygen.x ||
             aa2.atom[i]->coord.y != closestOxygen.y ||
             aa2.atom[i]->coord.z != closestOxygen.z)
            {
              otherOxygen = aa2.atom[i]->coord;
            }
        }
    }

  // Vector between benzene center of mass and closest oxygen
  Coordinates oD  = closestOxygen - aa1.center[0];
  // Vector between benzene center of mass and the other oxygen
  Coordinates oD2 = otherOxygen - aa1.center[0];
  
  // Oxygen dot product
  float oxy_numerator  = dotProduct(perp, oD); 
  float oxy_numerator2 = dotProduct(perp, oD2); 

  // distance from beneze center and the oxygens
  float distFromCenterToOxy  = oD.norm();
  float distFromCenterToOxy2 = oD2.norm();

  // Denominators
  float oxy_denom  = perpnorm * distFromCenterToOxy;
  float oxy_denom2 = perpnorm * distFromCenterToOxy2;

  if(oxy_denom == 0)
    {
      cerr << red << "Error" << reset << ": oxy_denom are zero.  Skipping residue" << endl;
      return false;
    }
  
  float u = num / denom;
  float uOxy  = oxy_numerator  / oxy_denom;
  float uOxy2 = oxy_numerator2 / oxy_denom2;

  // Force u to be [-1,1]
  if(u > 1)       u =  1;
  else if(u < -1) u = -1;

  if(uOxy > 1)       uOxy =  1;
  else if(uOxy < -1) uOxy = -1;

  if(uOxy2 > 1)       uOxy2 =  1;
  else if(uOxy2 < -1) uOxy2 = -1;

  // Get the angle and change it to degrees
  // take absolute vale to put it in [0,90]
  *dist      = distFromMassToChg;
  *distOxy   = distFromCenterToOxy;
  *distOxy2  = distFromCenterToOxy2;
  *angle     = fabs( 90 - acos(u)     * 180/3.14159 );
  *angleOxy  = fabs( 90 - acos(uOxy)  * 180/3.14159 );
  *angleOxy2 = fabs( 90 - acos(uOxy2) * 180/3.14159 );
  return true;
}

void Residue::markAltLocAtomsPHEorTYR(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " &&
          &this->atom[i]->coord != (this->center[index].plane_info[CG_PLANE_COORD_PTT]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CD1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[CD1_PLANE_COORD_PTT]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CD2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[CD2_PLANE_COORD_PTT]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CZ " &&
               &this->atom[i]->coord != (this->center[index].plane_info[3]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CE1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[4]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " CE2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[5]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void Residue::markAltLocAtomsASP(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " &&
          &this->atom[i]->coord != (this->center[index].plane_info[C__PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OD1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_1_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OD2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_2_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void Residue::markAltLocAtomsGLU(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CD " &&
          &this->atom[i]->coord != (this->center[index].plane_info[C__PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OE1" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_1_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " OE2" &&
               &this->atom[i]->coord != (this->center[index].plane_info[O_2_PLANE_COORD_AG]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void Residue::markAltLocAtomsPO4or2HPorPI(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " &&
          &this->atom[i]->coord != (this->center[index].plane_info[0]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O1 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[1]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O2 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[2]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O3 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[3]))
        {
          this->atom[i]->skip = true;
        }
      else if( this->atom[i]->name == " O4 " &&
               &this->atom[i]->coord != (this->center[index].plane_info[4]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void Residue::markAltLocAtoms2POorPO3(int index)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " &&
          &this->atom[i]->coord != (this->center[index].plane_info[0]))
        {
          this->atom[i]->skip = true;
        }
      else if( (this->atom[i]->name == " O1 " || this->atom[i]->name == " O1P") &&
               &this->atom[i]->coord != (this->center[index].plane_info[1]))
        {
          this->atom[i]->skip = true;
        }
      else if( (this->atom[i]->name == " O2 "  || this->atom[i]->name == " O2P") &&
               &this->atom[i]->coord != (this->center[index].plane_info[2]))
        {
          this->atom[i]->skip = true;
        }
      else if( (this->atom[i]->name == " O3 "  || this->atom[i]->name == " O2P") &&
               &this->atom[i]->coord != (this->center[index].plane_info[3]))
        {
          this->atom[i]->skip = true;
        }
    }  
}

void Residue::markAltLocAtoms(int index)
{
  if(residue == "PHE" || residue == "TYR")
    {
      markAltLocAtomsPHEorTYR(index);
    }
  else if(residue == "ASP")
    {
      markAltLocAtomsASP(index);
    }
  else if(residue == "GLU")
    {
      markAltLocAtomsGLU(index);
    }
  else if(residue == "PO4" || residue == "2HP" || residue == " PI")
    {
      markAltLocAtomsPO4or2HPorPI(index);
    }
  else if(residue == "2PO" || residue == "PO3")
    {
      markAltLocAtoms2POorPO3(index);
    }
}

void Residue::unmarkAltLocAtomsPHEorTYR()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CD1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CD2" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CZ " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CE1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " CE2" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void Residue::unmarkAltLocAtomsASP()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OD1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OD2" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void Residue::unmarkAltLocAtomsGLU()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CD " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OE1" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " OE2" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void Residue::unmarkAltLocAtomsPO4or2HPorPI()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O1 " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O2 " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O3 " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O4 " )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void Residue::unmarkAltLocAtoms2POorPO3()
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O1 "  || this->atom[i]->name == " O1P" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O2 "  || this->atom[i]->name == " O2P" )
        {
          this->atom[i]->skip = false;
        }
      else if( this->atom[i]->name == " O3 "  || this->atom[i]->name == " O3P" )
        {
          this->atom[i]->skip = false;
        }
    }  
}

void Residue::unmarkAltLocAtoms()
{
  if(residue == "PHE" || residue == "TYR")
    {
      unmarkAltLocAtomsPHEorTYR();
    }
  else if(residue == "ASP")
    {
      unmarkAltLocAtomsASP();
    }
  else if(residue == "GLU")
    {
      unmarkAltLocAtomsGLU();
    }
  else if(residue == "PO4" || residue == "2HP" || residue == " PI")
    {
      unmarkAltLocAtomsPO4or2HPorPI();
    }
  else if(residue == "2PO" || residue == "PO3")
    {
      unmarkAltLocAtoms2POorPO3();
    }
}

string Residue::makeConectPHEorTYR()
{
  string serials[6];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CD1"  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CD2"  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CZ " && !this->atom[i]->skip )
        {
          serials[3] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CE1" && !this->atom[i]->skip )
        {
          serials[4] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " CE2" && !this->atom[i]->skip )
        {
          serials[5] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[2] + serials[1] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[4] + serials[0] + "                                                 \n";
  conect += "CONECT" + serials[2] + serials[5] + serials[0] + "                                                 \n";
  conect += "CONECT" + serials[3] + serials[5] + serials[4] + "                                                 \n";
  conect += "CONECT" + serials[4] + serials[1] + serials[3] + "                                                 \n";
  conect += "CONECT" + serials[5] + serials[2] + serials[3] + "                                                 \n";
  return conect;
}

string Residue::makeConectPHEorTYR_altloc(int c)
{
  string serials[6];
  for(int i=0; i<this->altlocs[c].size(); i++)
    {
      if( this->altlocs[c][i]->name == " CG " && !this->altlocs[c][i]->skip )
        {
          serials[0] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " CD1"  && !this->altlocs[c][i]->skip )
        {
          serials[1] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " CD2"  && !this->altlocs[c][i]->skip )
        {
          serials[2] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " CZ " && !this->altlocs[c][i]->skip )
        {
          serials[3] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " CE1" && !this->altlocs[c][i]->skip )
        {
          serials[4] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " CE2" && !this->altlocs[c][i]->skip )
        {
          serials[5] = this->altlocs[c][i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[2] + serials[1] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[4] + serials[0] + "                                                 \n";
  conect += "CONECT" + serials[2] + serials[5] + serials[0] + "                                                 \n";
  conect += "CONECT" + serials[3] + serials[5] + serials[4] + "                                                 \n";
  conect += "CONECT" + serials[4] + serials[1] + serials[3] + "                                                 \n";
  conect += "CONECT" + serials[5] + serials[2] + serials[3] + "                                                 \n";
  return conect;
}

string Residue::makeConectGLU()
{
  string serials[3];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CD " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OE1"  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OE2"  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  return conect;
}

string Residue::makeConectGLU_altloc(int c)
{
  string serials[3];
  for(int i=0; i<this->altlocs[c].size(); i++)
    {
      if( this->altlocs[c][i]->name == " CD " && !this->altlocs[c][i]->skip )
        {
          serials[0] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " OE1"  && !this->altlocs[c][i]->skip )
        {
          serials[1] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " OE2"  && !this->altlocs[c][i]->skip )
        {
          serials[2] = this->altlocs[c][i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  return conect;
}


string Residue::makeConectASP()
{
  string serials[3];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " CG " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OD1"  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " OD2"  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  return conect;
}

string Residue::makeConectASP_altloc(int c)
{
  string serials[3];
  for(int i=0; i<this->altlocs[c].size(); i++)
    {
      if( this->altlocs[c][i]->name == " CG " && !this->altlocs[c][i]->skip )
        {
          serials[0] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " OD1"  && !this->altlocs[c][i]->skip )
        {
          serials[1] = this->altlocs[c][i]->line.substr(6,5);
        }
      else if( this->altlocs[c][i]->name == " OD2"  && !this->altlocs[c][i]->skip )
        {
          serials[2] = this->altlocs[c][i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + "                                                 \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  return conect;
}


string Residue::makeConectPO4or2HPorPI()
{
  string serials[5];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( this->atom[i]->name == " P  " && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O1 "  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O2 "  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O3 "  && !this->atom[i]->skip )
        {
          serials[3] = this->atom[i]->line.substr(6,5);
        }
      else if( this->atom[i]->name == " O4 "  && !this->atom[i]->skip )
        {
          serials[4] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + serials[3] + serials[4] + "\n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[3] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[4] + serials[0] + "                                                       \n";
  return conect;
}

string Residue::makeConect2POorPO3()
{
  string serials[5];
  for(int i=0; i<this->atom.size(); i++)
    {
      if( (this->atom[i]->name == " P  ") && !this->atom[i]->skip )
        {
          serials[0] = this->atom[i]->line.substr(6,5);
        }
      else if( (this->atom[i]->name == " O1 " || this->atom[i]->name == " O1P" )  && !this->atom[i]->skip )
        {
          serials[1] = this->atom[i]->line.substr(6,5);
        }
      else if( (this->atom[i]->name == " O2 " || this->atom[i]->name == " O1P" )  && !this->atom[i]->skip )
        {
          serials[2] = this->atom[i]->line.substr(6,5);
        }
      else if( (this->atom[i]->name == " O3 " || this->atom[i]->name == " O1P" )  && !this->atom[i]->skip )
        {
          serials[3] = this->atom[i]->line.substr(6,5);
        }
    }
  string conect = "CONECT" + serials[0] + serials[1] + serials[2] + serials[3] + "      \n";
  conect += "CONECT" + serials[1] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[2] + serials[0] + "                                                       \n";
  conect += "CONECT" + serials[3] + serials[0] + "                                                       \n";
  return conect;
}

string Residue::makeConectCarbonRing(int c)
{
  string serials[6];
  string conect = "";

  for(int i=0; i<this->altlocs[c].size(); i++)
    {
      conect += "CONECT" + altlocs[c][i]->line.substr(6,5);

      int nConnected = 0;
      for(int j=0; j<this->altlocs[c].size(); j++)
        {
          float dist = altlocs[c][i]->coord.distance(altlocs[c][j]->coord);
//cout << dist << " ";
          if (dist >= 1.2 && dist <= 1.6)
            {
              nConnected++;
              conect += altlocs[c][j]->line.substr(6,5);
            }
        }
//cout << endl;

      if (nConnected != 2)
        {
          cout << "Problem with carbon ring: One carbon is connected to " << nConnected;
          cout << " others. " << endl;
        }

      conect += "      \n";
      
    }

/*
  for(int i=0; i<this->carbonRings[c].size(); i++)
    {
      conect += "CONECT" + carbonRings[c][i]->line.substr(6,5);

      int nConnected = 0;
      for(int j=0; j<this->carbonRings[c].size(); j++)
        {
          float dist = carbonRings[c][i]->coord.distance(carbonRings[c][j]->coord);
//cout << dist << " ";
          if (dist >= 1.2 && dist <= 1.6)
            {
              nConnected++;
              conect += carbonRings[c][j]->line.substr(6,5);
            }
        }
//cout << endl;

      if (nConnected != 2)
        {
          cout << "Problem with carbon ring: One carbon is connected to " << nConnected;
          cout << " others. " << endl;
        }

      conect += "      \n";
      
    }
*/
  
  return conect;
}

string Residue::makeConect(int c)
{
  if(residue == "PHE" || residue == "TYR")
    {
      return makeConectPHEorTYR_altloc(c);
    }
  else if(residue == "ASP")
    {
      return makeConectASP();
    }
  else if(residue == "GLU")
    {
      return makeConectGLU_altloc(c);
    }
  else if(residue == "PO4" || residue == "2HP" || residue == " PI")
    {
      return makeConectPO4or2HPorPI();
    }
  else if(residue == "2PO" || residue == "PO3")
    {
      return makeConect2POorPO3();
    }
  else if(carbonRings.size() != 0)
    {
      return makeConectCarbonRing(c);
    }
}

// This checks the validity of the hydrogens for 
// GLU and ASP residues.  Sometimes babel will add 4 hydrogens:
// 2 to C and 1 to each O.  We can figure out why, but we came
// up with this hackish fix.  We take the 2 that are connected
// to the carbon, average them, and use that as the coordinates 
// for the hydrogen that we are looking for.  We then throw all
// of the other hydrogens away.
bool Residue::removeExcessHydrogens(vector<string> conect)
{
  if( !(residue == "GLU" || residue == "ASP") )
    return false;
  
  vector<Atom*>::iterator it;
  int carbonSerialNumber;
  int hydrogenCount = 0;
  for(it = atom.begin(); it != atom.end(); ++it)
    {
      if( (*it)->name == " CG " || (*it)->name == " CD " )
        {
          carbonSerialNumber = (*it)->serialNumber;
        }
      else if( (*it)->name == " H  " )
        {
          hydrogenCount++;
        }
    }
  corrected = false;
  if( hydrogenCount > 1 )
    {
      corrected = true;
      Coordinates avg(0,0,0);
      vector<string>::iterator conectit;
      Atom* lastHydrogen;
      for(it = atom.begin(); it < atom.end(); ++it)
        {
          for(conectit = conect.begin(); conectit < conect.end(); ++conectit)
            {
              int number;
              int cnumber;
              if(!from_string<int>(number,(*conectit).substr(7,5),dec))
                {
                  cerr << red << "Error" << reset << ": failed to convert serial number into an int"  << endl;
                }
              string tempstr = (*conectit).substr(11,5);
              if(tempstr != "     ")
                {
                  if(!from_string<int>(cnumber,(*conectit).substr(11,5),dec))
                    {
                      cerr << red << "Error" << reset << ": failed to convert serial number into an int "  << endl;
                    }

                  if( (*it)->name == " H  " && (*it)->serialNumber == number 
                      && cnumber == carbonSerialNumber )
                    {
                      avg += (*it)->coord;
                    }
                }
            }
          if((*it)->name == " H  " )
            {
              lastHydrogen = (*it);
              vector<Atom*>::iterator tempit = it-1;
              atom.erase(it);
              it = tempit;
            }
        }
      avg /= 2;
      lastHydrogen->coord = avg;
      char cstr[25];
      sprintf(cstr,"%8.3lf%8.3lf%8.3lf",avg.x, avg.y, avg.z);
      string temp(cstr);
      lastHydrogen->line = lastHydrogen->line.substr(0,30) + 
        temp.substr(0,24) + lastHydrogen->line.substr(54,26);
      atom.push_back(lastHydrogen);
    }
  return corrected;
}

void Residue::printPHEorTYR(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " CD1" || 
          atom[i]->name == " CD2" || 
          atom[i]->name == " CE1" || 
          atom[i]->name == " CE2" || 
          atom[i]->name == " CZ " ||
          atom[i]->name == " CG ") && !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void Residue::printASP(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " OD1" || 
         atom[i]->name == " OD2" || 
         atom[i]->name == " CG " ) && !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void Residue::printGLU(FILE* output)
{
  for(int i=0; i < this->atom.size(); i++)
    {
      if((atom[i]->name == " OE1" || 
         atom[i]->name == " OE2" || 
         atom[i]->name == " CD " )&& !atom[i]->skip)
        {
          this->atom[i]->print(output);
        }
    }
}

void Residue::printNeededAtoms(FILE* output)
{

  if(residue == "PHE" || residue == "TYR")
    {
      printPHEorTYR(output);
    }
  else if(residue == "ASP")
    {
      printASP(output);
    }
  else if(residue == "GLU")
    {
      printGLU(output);
    }
  else
    {
      cerr << red << "Error" << reset << ": Unsupported residue" << endl;
    }
}

ostream& operator<<(ostream& output, const Residue& p) 
{
  for(int i=0; i < p.atom.size(); i++)
    {
      output << *(p.atom[i]) << endl;
    }
  return output;  // for multiple << operators.
}

/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Atom.cpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: This file contains all the class function implementations for
//               the Atom class. Includes constructors, destructor, print
//               function, and PDB ATOM line parser. This not only works for
//               ATOM lines, but it also works for HETATM lines since they have
//               the same format.
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

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "Atom.hpp"
#include "Coordinates.hpp"
#include "CoutColors.hpp"

// Constructor setting everything to initial values
Atom::Atom()
{
  serialNumber = 0;
  name = "";
  altLoc = '\0';
  residueName = "";
  chainID = '-';
  resSeq = 0;
  iCode = ' ';
  occupancy = 0.0;
  tempFactor = 0.0;
  element = "";
  charge = "";
  failure = false;
  skip = false;
}

// Constructor that parses an ATOM line from a PDB file
Atom::Atom(string line, int num)
{
  parseAtom(line, num);
}

// Constructor that parses an ATOM line from a PDB file
Atom::Atom(char* linecs, int num)
{
  string line(linecs);
  parseAtom(line, num);

}

// Destructor setting everything to initial values
Atom::~Atom()
{
  serialNumber = 0;
  name = "";
  altLoc = '\0';
  residueName = "";
  chainID = '-';
  resSeq = 0;
  iCode = ' ';
  occupancy = 0.0;
  tempFactor = 0.0;
  element = "";
  charge = "";
  failure = false;
  skip = false;
}

// Returns true if the parsing failed, false otherwise
bool Atom::fail()
{
  return failure;
}

// Print out all the values space separated
void Atom::print()
{
  cout << fixed << setprecision(3)
       << serialNumber << " " << name   << " " << altLoc << " " << residueName << " "
       << chainID      << " " << resSeq << " " << iCode  << " " << coord.x << " "
       << coord.y << " " << coord.z << " " << occupancy << " " << tempFactor << " "
       << element << " " << charge << " " << endl;
}

// Parse the ATOM line of a PDB file
void Atom::parseAtom(string line, int num)
{
  failure = false;
  this->line = line;

  if(this->line[line.length()-1] == '\r')
    {
      this->line.erase(line.length()-1,1);
    }

  // Error check to ensure the file is formatted correctly
  if(this->line.length() != 80)
    {
      cerr << red << "Error" << reset << ": Possible malformed PDB file on ATOM/HETATM line " << num << "." << endl;
cout << line << endl;
      failure = true;
    }

  // Get the serial number and error check
  if(!from_string<int>(serialNumber,this->line.substr(6,5),dec))
    {
      cerr << red << "Error" << reset << ": failed to convert serial number into an unsigned int on line " << num << endl;
      failure = true;
    }

  // Grab the name, alternate location, residue name, and chain ID
  name = this->line.substr(12,4);
  altLoc = this->line[16];
  residueName = this->line.substr(17,3);
  chainID = this->line[21];
  if(chainID == ' ')
    {
      chainID = 'A';
    }
  // Grab the residue sequence number
  if(!from_string<int>(resSeq,this->line.substr(22,4),dec))
    {
      cerr << red << "Error" << reset << ": failed to convert residue sequence into an unsigned int " << num << endl;
      failure = true;
    }

  // Grab the insertion code
  iCode = line[26];

  // Grab the coordinates, occupancy and temperature factor
  if(!from_string<float>(coord.x,this->line.substr(30,8),dec))
    {
      cerr << red << "Error" << reset << ": failed to convert x coordinate into a double" << endl;
      failure = true;
    }
  if(!from_string<float>(coord.y,this->line.substr(38,8),dec))
    {
      cerr << red << "Error" << reset << ": failed to convert y coordinate into a double" << endl;
      failure = true;
    }
  if(!from_string<float>(coord.z,this->line.substr(46,8),dec))
    {
      cerr << red << "Error" << reset << ": failed to convert z coordinate into a double" << endl;
      failure = true;
    }
  if(!from_string<double>(occupancy,this->line.substr(54,6),dec))
    {
      cerr << red << "Error" << reset << ": failed to convert occupancy into a double" << endl;
      failure = true;
    }
  if(!from_string<double>(tempFactor,this->line.substr(60,8),dec))
    {
      cerr << red << "Error" << reset << ": failed to convert temperature factor into a double" << endl;
      failure = true;
    }

  // Store the element and charge
  element = this->line.substr(76,2);
  charge = this->line.substr(78,2);

  // Set the atom number
  if (element == " C")
    {
      element_num = C;
    }
  else if (element == " H")
    {
      element_num = H;
    }
  else if (element == " N")
    {
      element_num = N;
    }
  else if (element == " O")
    {
      element_num = O;
    }
  else if (element == " F")
    {
      element_num = F;
    }
  else if (element == " P")
    {
      element_num = P;
    }
  else if (element == " S")
    {
      element_num = S;
    }
  else if (element == "CL" || element == "Cl")
    {
      element_num = Cl;
    }
  else if (element == "BR" || element == "Br")
    {
      element_num = Br;
    }
  else
    {
      cerr << "Warning: Element unaccounted for: " << element << endl;
    }

}
// Outputs the ATOM line into a file
void Atom::print(FILE* output)
{
  // fprintf(output,
  //      "ATOM  %5u %-4s%c%3s %c%4u%c   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s%-s\n",
  //      this->serialNumber,
  //      this->name.c_str(),
  //      this->altLoc,
  //      this->residueName.c_str(),
  //      this->chainID,
  //      this->resSeq,
  //      this->iCode,
  //      this->coord.x,
  //      this->coord.y,
  //      this->coord.z,
  //      this->occupancy,
  //      this->tempFactor,
  //      this->element.c_str(),
  //      this->charge.c_str());
  fprintf(output, "%s\n",this->line.c_str());
}

// Outputs the ATOM line to a stream
ostream& operator<<(ostream& output, const Atom& p) 
{
  // output << "ATOM  "
  //     << setw(5) << right << p.serialNumber << " "
  //     << setw(4) << left  << p.name
  //     << p.altLoc
  //     << setw(3) << p.residueName << " "
  //     << p.chainID
  //     << setw(4) << right << p.resSeq
  //     << p.iCode << "   "
  //     << fixed
  //     << setw(8) << right << setprecision(3) << p.coord.x
  //     << setw(8) << right << setprecision(3) << p.coord.y
  //     << setw(8) << right << setprecision(3) << p.coord.z
  //     << setw(6) << right << setprecision(2) << p.occupancy
  //     << setw(6) << right << setprecision(2) << p.tempFactor
  //     << "          ";
  // output << setw(2) << right << p.element
  //     << setw(2) << left  << p.charge;
  output << p.line;

  return output;  // for multiple << operators.
}

bool Atom::isBonded(Atom &other)
{
  float dist = coord.distance(other.coord);
  float upper_bound = COVALENT_RADIUS[element_num] +
                      COVALENT_RADIUS[other.element_num] +
                      BOND_TOLERANCE_FACTOR;

  if (dist <= upper_bound)
    {
      return true;
    }
  else
    {
      return false;
    }
}

/*
  If the two atoms have the same resSeq, then the 'smaller'
  atom has the smaller chainID. Otherwise, if the atoms are
  both HETATMs or ATOMs, then the 'smaller' atom has the
  smaller resSeq.  Finally, the ATOMs are 'smaller' than the
  HETATMs.
*/
bool Atom::operator<(const Atom &rhs) const
{
  bool ret_val;
  if(this->resSeq == rhs.resSeq)
    {
      return (this->chainID < rhs.chainID);
    }
  else
    {
      if (line.substr(0, 6) == rhs.line.substr(0, 6))
        {
          return this->resSeq < rhs.resSeq;
        }
      else if (line.substr(0, 4) == "ATOM")
        {
          return true;
        }
      else
        {
          return false;
        }
    }
}

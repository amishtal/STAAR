/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: Options.cpp
//  Date: 10 Jan 2011
//  Version: 1.0
//  Description: Implementation file for the command line parser and Options class
//
//  Updates: Added the ability to input residues of interest through the cmd line (11 Feb 2011)
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

#include "Options.hpp"
#include "CoutColors.hpp"

// Initialize the Options to empty stuff
Options::Options()
{
  center = false;
  pdbfile = NULL;
  outputfile = NULL;
  failure = false;
  sameChain = false;
  threshold = 15.0;
}

// Intialize options then parse the cmd line arguments
Options::Options( int argc, char **argv )
{
  center = false;
  pdbfile = NULL;
  outputfile = NULL;
  failure = false;
  sameChain = false;
  threshold = 15.0;
  parseCmdline( argc, argv );
}

// Empty destructor
Options::~Options(){}

// Print a help page
void printHelp()
{
  cerr << "-h or --help          " << "Displays this message" << endl;
  cerr << "-p or --pdbdir        " << "Specifies the folder for PDB files" << endl;
  cerr << "-o or --out           " << "Specifies the output file" << endl;
  cerr << "-r or --residues      " << "Set the residues that we are going to analyze" << endl;
  cerr << "                      " << " residues are set as follows (include quotations):" << endl;
  cerr << "                      " << " \"PHE,TYR,TRP;GLU,ASP\"" << endl;
  cerr << "                      " << " at least 1, but as many as you want" << endl;
  cerr << "-c or --center        " << "Specifies whether to calculate center of charge" << endl;
  cerr << "-t or --threshold     " << "Set the distance threshold between amino acids" << endl;
  cerr << "-s or --samechain     " << "Look for interactions in same chain only" << endl;
  cerr << "-l or --ligand        " << "Look for interactions with ligands" << endl;
  cerr << "                      " << " add in the residue name looking for in HETATM:" << endl;
  cerr << "                      " << " \"2HP,WO2\"" << endl;
}

// Return true of cmd line parsing failed, false otherwise
bool Options::fail()
{
  return failure;
}

// Parses the cmd line options
void Options::parseCmdline( int argc, char **argv )
{
  
  int c;
  // Set the options that are possible with this program
  static struct option long_options[] = 
    {
      {"help",		required_argument, 0, 'h'},
      {"pdbdir",	required_argument, 0, 'p'},
      {"out",		required_argument, 0, 'o'}, 
      {"center",	no_argument,       0, 'c'},
      {"threshold",	required_argument, 0, 't'},
      {"samechain",     no_argument,       0, 's'},
      {"residues",      required_argument, 0, 'r'},
      {0, 0, 0, 0}
    };
  int option_index;

  // Go through the options and set them to variables
  while( !( ( c = getopt_long(argc, argv, "hp:o:ct:sr:", long_options, &option_index) ) < 0 ) )
    {
    switch(c)
      {
      case 'h':
	// the user is asking for help here
        printHelp();
	exit(1);
        break;

      case 'p':
	// the user is specifying the PDB file/directory
        this->pdbfile = optarg;
        break;

      case 'o':
	// the user is wanting to set the output file, but we need to make
	// sure that is ia a file and not a directory
	if( !isDirectory( optarg ) )
	  {
	    this->outputfile = optarg;
	  }
	else
	  {
	    cerr << red << "Error" << reset << ": Output file must be a file, not a directory" << endl;
	    printHelp();
	    exit(1);
	  }
        break;

      case 'c':
	// user wants to do the center of charge calculations, but really
	// there is no option at this point because we need to fix the 
	// setting of the plane variables
        this->center = true;
        break;

      case 't':
	// user wants to set the threshold value but we need to make
	// sure that an actual number was inputted and not just a string
	if ( sscanf(optarg, "%f", &(this->threshold)) != 1 )
	  {
	    cerr << red << "Error" << reset << ": Must insert a valid number for the threshold" << endl;
	    printHelp();
	    exit(1);
	  }
        break;

      case 's':
	// user wants to look for interactions within the same chains only
	this->sameChain = true;
        break;

      case 'r':
	{
	  vector<string> temp = split( (string)optarg, ';' );
	  if(temp.size() != 2)
	    {
	      cerr << red << "Error" << reset << ": Residue string must have 1 semicolon. No more, no less." << endl;
	      printHelp();
	      exit(1);
	    }
	  this->residue1 = split( temp[0], ',' );
	  this->residue2 = split( temp[1], ',' );
	}
	break;

      default:
	printHelp();
	exit(1);
	break;
      }
    }

  // Error check to make sure there wasn't extra crap on the command line
  // and all necessary arguments were set
  if( optind < argc ){
    cerr << red << "Error" << reset << ": There are some extra commands inputted. Please cleanup the command line!" << endl;
    // print some help info
    printHelp();
    failure=true;
  }else if( !pdbfile ){
    cerr << red << "Error" << reset << ": Must specify the PDB list file with -p or --pdblist" <<  endl;
    printHelp();
    failure=true;
  }else if( !outputfile ){
    cerr << red << "Error" << reset << ": Must specify the op file with -o or --op" <<  endl;
    printHelp();
    failure=true;
  }

}



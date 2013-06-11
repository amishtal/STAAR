/****************************************************************************************************/
//  COPYRIGHT 2011, University of Tennessee
//  Author: David Jenkins (david.d.jenkins@gmail.com)
//  File: staar.cpp
//  Date: 12 Jan 2011
//  Version: 1.0
//  Description: This is the STAAR program. It will take a PDB file/directory parsing through it
//               looking for possible anion-quadrople interactions.  It will output into a
//               specified file all possible interactions that are within a specified distance
//               threshold and the angles at which the two residues are interacting at.
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

#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <float.h>

#include "Utils.hpp"
#include "Options.hpp"
#include "PDB.hpp"
#include "Seqres.hpp"
#include "Geometry.hpp"
#include "Residue.hpp"
#include "Coordinates.hpp"
#include "CoutColors.hpp"
#include "StatisticsRecorder.hpp"

#define MAX_STR_LENGTH 1024

// Parses through a single PDB file checking for interactions
bool processSinglePDBFile(const char* filename,
                          Options& opts,
                          ofstream& output_file,
                          const char* chains,
                          StatisticsRecorder* stats);

// Read PDB names from a list and parses them from the specified directory
bool processPDBList(Options& opts, StatisticsRecorder* stats);

// Reads the PDB names and chains from a list and parses them from the
// specified directory
bool processPDBChainList(Options& opts, StatisticsRecorder* stats);

// Traverses through a directory of PDB files processing each one
bool processPDBDirectory(Options& opts, StatisticsRecorder* stats);

// Searches through all of the chains looking for interactions
void searchChainInformation(PDB & PDBfile,
                            unsigned int chain1,
                            unsigned int chain2,
                            string residue1,
                            string residue2,
                            Options & opts,
                            ofstream& output_file);

void searchLigandsInformation(PDB & PDBfile,
                              Residue & ligand,
                              unsigned int chain1,
                              string residue1,
                              Options & opts,
                              ofstream& output_file);

void searchCarbonRingLigandsInformation(PDB & PDBfile,
                                        Residue & ligand,
                                        unsigned int chain1,
                                        string residue1,
                                        Options & opts,
                                        ofstream& output_file,
                                        StatisticsRecorder *stats);

// Finds the closest distance among all of the centers
// associated with each amino acid
double findClosestDistance(Residue& aa1,
                           Residue& aa2,
                           float threshold,
                           unsigned int* closest_index1,
                           unsigned int* closest_index2);

// Finds the closest interaction among all of the possible
// amino acid centers
bool findBestInteraction( Residue& aa1,
                          Residue& aa2,
                          float threshold,
                          PDB& PDBfile,
                          char* gamessfolder,
                          bool ligand,
                          ofstream& output_file);

// Writes the INP files
void outputINPfile(string input_filename,
                   char* filename,
                   Residue& aa1h,
                   Residue& aa2h);

void write_output_head(ofstream& out);

int main(int argc, char* argv[]){
  printHeader();
  int return_value;
  double start = getTime();

  StatisticsRecorder* stats = new StatisticsRecorder();

  // Parse and error check the command line arguments
  Options opts(argc, argv);

  // Looks like something is wrong with the cmd line args
  // or the user input the -h flag
  if(opts.fail())
    {
      return 1;
    }

  // Lets check to see what kind of input we received for the
  // PDB file.  If we have a list file of PDBs, just go through
  // list with the specified directory. If a directory, parse 
  // all files in the directory. Otherwise just parse the single 
  // file
  if( opts.pdblist )
    {
      return_value = processPDBList(opts, stats);
    }
  else if( opts.chain_list )
    {
      return_value = processPDBChainList(opts, stats);
    }
  else if( isDirectory(opts.pdbfile) )
    {
      // go through each file in the directory
      return_value = processPDBDirectory(opts, stats);
    }
  else
    {
      ofstream output_file(opts.outputfile);
      if(!output_file)
        {
          cerr << red << "Error" << reset << ": Failed to open output file," << opts.outputfile << endl;
          return 1;
        }
      write_output_head(output_file);
      return_value = processSinglePDBFile(opts.pdbfile, opts, output_file, NULL, stats);
      output_file.close();
    }

  cout << "Time taken: " << getTime() - start << "s" << endl;
  cout << "Total interactions: " << stats->getTotalNumberOfInteractions() << endl;
  cout << "Average interactions per protein: " << stats->averageInteractionsPerProtein() << endl;

  vector<int> ring_interaction_counts = stats->carbonRingLigandInteractionCounts();
  vector<int> ring_counts = stats->carbonRingLigandCounts();
  if (ring_counts.size() != ring_interaction_counts.size())
    {
      ring_interaction_counts.resize(ring_counts.size());
    }
  cout << "Interactions involving ligands with carbon rings:" << endl;
  cout << setw(7) << "# Rings" << "  " << setw(14) << "# Interactions" << "  " << setw(15) << "# Total Ligands" << endl;
  for (int i=0; i < ring_counts.size(); i++)
    {
      cout << setw(7) << i+1 << "  " << setw(14) << ring_interaction_counts[i] << "  " << setw(15) << ring_counts[i] << endl;
    }

  return !return_value;
}

bool processSinglePDBFile(const char* filename,
                          Options& opts,
                          ofstream& output_file,
                          const char* chains,
                          StatisticsRecorder* stats)
{
  int numRes1 = opts.residue1.size();
  int numRes2 = opts.residue2.size();

  // Read in the PDB file
  // This function actually reads and stores more information than
  // we really need, but the files are relatively small so it isn't
  // taking up much RAM (from what I saw, < 5MB each PDB file)
  // with a few exceptions
  PDB PDBfile_whole(filename, opts.resolution);

  if( PDBfile_whole.fail() )
    {
      //cerr << red << "Error" << reset << ": Parsing PDB file failed!" << endl;
      PDBfile_whole.printFailure();
      return false;
    }


  stats->setActiveProtein(PDBfile_whole.idCode);
//cout << stats->getActiveProteinName() << endl;

  for(unsigned int model=0; model < PDBfile_whole.models.size(); model++)
    {
      PDB PDBfile = PDBfile_whole.models[model];
      PDBfile.filename = PDBfile_whole.filename;
      PDBfile.resolution = PDBfile_whole.resolution;

      if( opts.findCarbonRings || opts.numLigands )
        {
          PDBfile.keepLigands = true;
        }
      else
        {
          PDBfile.keepLigands = false;
        }

      PDBfile.setResiduesToFind(&opts.residue1, &opts.residue2);
      if(opts.numLigands)
        {
          PDBfile.setLigandsToFind(&opts.ligands);
        }

      PDBfile.populateChains(false);

      if( opts.numLigands )
        {
          PDBfile.findLigands( opts.ligands );
        }
      if(opts.findCarbonRings)
        {
          PDBfile.findCarbonRingsInLigands();
        }

      // Searching for interations within each chain
      for(unsigned int i = 0; i < PDBfile.chains.size(); i++)
        {
          // Check if we are supposed to look at certain chains
          if( chains )
            {
              if( !strchr(chains,PDBfile.chains[i].id) )
                {
                  continue;
                }
            }
          // if we want to look for interactions between the ith chain
          // and each of the other chains, we set the indices to loop
          // though all the chains
          unsigned int start = 0;
          unsigned int end   = PDBfile.chains.size();

          // otherwise we just set the indices to go through the ith chain
          if(opts.sameChain)
            {
              start = i;
              end   = i+1;
            }

          for(unsigned int j = start; j < end; j++)
            {
              for(int ii = 0; ii < numRes1; ii++)
                {
                  for(int jj = 0; jj < numRes2; jj++)
                    {
                      // Search for different residue combinations
                      searchChainInformation(PDBfile,
                                             i,
                                             j,
                                             opts.residue1[ii],
                                             opts.residue2[jj],
                                             opts,
                                             output_file);
                    }
                }
            }
      
          // Go through the ligands, if there are any
          for(unsigned int j = 0; j<PDBfile.ligands.size(); j++)
            {
              for(int ii=0; ii<numRes1; ii++)
                {
                  searchLigandsInformation(PDBfile,
                                           *PDBfile.ligands[j],
                                           i,
                                           opts.residue1[ii],
                                           opts,
                                           output_file);
                }
            }

          // Go through the ligands containing carbon rings
          for(unsigned int j = 0; j<PDBfile.carbonRingLigands.size(); j++)
            {
              for(int ii=0; ii<numRes2; ii++)
                {
                  searchCarbonRingLigandsInformation(PDBfile,
                                                     *PDBfile.carbonRingLigands[j],
                                                     i,
                                                     opts.residue2[ii],
                                                     opts,
                                                     output_file,
                                                     stats);
                }
            }
        }
    }
  return true;
}

bool processPDBList(Options& opts, StatisticsRecorder* stats)
{
  // Open the list file
  ifstream listfp(opts.pdblist);
  if( !listfp )
    {
      cerr << red << "Error" << reset << ": Failed to open list file, " << opts.pdblist << endl;
      perror("\t");
      return false;
    }

  // Open the output file
  ofstream output_file(opts.outputfile);
  if( !output_file )
    {
      cerr << red << "Error" << reset << ": Failed to open output file" << endl;
      perror("\t");
    }

  // Write the header to output file
  write_output_head(output_file);

  string line;

  // Go through each line of the PDB list file
  while(getline(listfp, line))
    {
      // Get and print the file name
      string filename(opts.pdbfile);
      filename += "/" + line + opts.extension;
      cout << purple << line << opts.extension << endl;

      // Process the PDB file
      processSinglePDBFile(filename.c_str(), opts, output_file, NULL, stats);
    }
  cout << endl;
  output_file.close();
  listfp.close();
  return true;
}

bool processPDBChainList(Options& opts, StatisticsRecorder* stats)
{
  // Open up the chain list file
  ifstream listfp(opts.chain_list);
  if( !listfp )
    {
      cerr << red << "Error" << reset << ": Failed to open list file, " << opts.chain_list << endl;
      perror("\t");
      return false;
    }

  // Open up the output file
  ofstream output_file(opts.outputfile);
  if( !output_file )
    {
      cerr << red << "Error" << reset << ": Failed to open output file" << endl;
      perror("\t");
    }

  // Write the header to the ouput file
  write_output_head(output_file);

  string line;

  // Go through each line of the list file
  while(getline(listfp, line))
    {
      // Store the PDB directory path first
      string filename(opts.pdbfile);

      // Separate the line into files
      vector<string> fields = split(line,'\t');
      if( fields.size() < 2)
        {
          cerr << "Chain file is malformed!" << endl;
          return false;
        }

      // Create and output the file name
      filename += "/" + fields[0] + opts.extension;
      cout << purple << fields[0] << opts.extension << endl;

      // Process the file
      processSinglePDBFile(filename.c_str(), opts, output_file, fields[1].c_str(), stats);
    }
  cout << endl;
  output_file.close();
  listfp.close();
  return true;
}

bool processPDBDirectory(Options& opts, StatisticsRecorder* stats)
{
  DIR* directory;
  struct dirent* filename;
  int numberOfFiles;
  int count=0;
  ofstream output_file(opts.outputfile);
  if( !output_file )
    {
      cerr << red << "Error" << reset << ": Failed to open output file" << endl;
      perror("\t");
    }
  write_output_head(output_file);

  // Open the directory for traversal
  if( (directory = opendir( opts.pdbfile )) )
    {
      // Go through each file of the directory
      while( (filename = readdir( directory )) )
        {
          if( strcmp(filename->d_name, ".") != 0 && strcmp(filename->d_name, "..") != 0 )
            {
              count++;
              // Checkpoint
              cout << purple << filename->d_name << endl;

              //create the full path to the file
              char fullFilePath [MAX_STR_LENGTH];
              sprintf(fullFilePath, "%s/%s", opts.pdbfile, filename->d_name);

              // perform some work on the current file
              processSinglePDBFile(fullFilePath, opts, output_file, NULL, stats);
            }
        }
    }
  else
    {
      perror(opts.pdbfile);
      return false;
    }
  closedir(directory);
  output_file.close();
  return true;
}

void searchChainInformation(PDB & PDBfile,
                            unsigned int chain1,
                            unsigned int chain2,
                            string residue1,
                            string residue2,
                            Options & opts,
                            ofstream& output_file)
{

  Chain* c1 = &(PDBfile.chains[chain1]);
  Chain* c2 = &(PDBfile.chains[chain2]);
#ifdef DEBUG
  unsigned int length_chain1 = c1->seqres[0]->numberOfResidues;
  unsigned int length_chain2 = c2->seqres[0]->numberOfResidues;

  cout << purple << "Checking for res1= "<< residue1 << ", res2= "<< residue2 << endl;
       << " in chains " << c1->id << " and " << c2->id << "..." << endl;

  cout << purple << "length(Chain1)= " << length_chain1
       << " length(Chain2)= " << length_chain2 << endl;
#endif
  // Go through each AA in the first chain
  for(unsigned int i = 0; i < c1->aa.size(); i++)
    {
      // If this chain has a residue that we are looking for,
      // let's do some analysis!
      if( c1->aa[i].residue.compare(residue1) == 0 && !(c1->aa[i].skip))
        {
          for(unsigned int j = 0; j < c2->aa.size(); j++)
            {              
              if(c2->aa[j].residue == residue2 && !(c2->aa[j].skip))
                {
                  // Find the best interaction out of all the centers
                  // for this AA pair
                  findBestInteraction( c1->aa[i],
                                       c2->aa[j],
                                       opts.threshold,
                                       PDBfile,
                                       opts.gamessfolder,
                                       false,
                                       output_file);
                }
            }
        }
    }
#ifdef DEBUG
  cout << purple << "\tDone" << endl;
#endif
}

void searchLigandsInformation(PDB & PDBfile,
                              Residue & ligand,
                              unsigned int chain1,
                              string residue1,
                              Options & opts,
                              ofstream& output_file)
{
  Chain* c1 = &(PDBfile.chains[chain1]);
  for(int i=0; i<c1->aa.size(); i++)
    {
      if(c1->aa[i].residue == residue1 && !(c1->aa[i].skip))
        {
          findBestInteraction(c1->aa[i],
                              ligand,
                              opts.threshold,
                              PDBfile,
                              opts.gamessfolder,
                              true,
                              output_file);
        }
    }
}

void searchCarbonRingLigandsInformation(PDB & PDBfile,
                                        Residue & ligand,
                                        unsigned int chain1,
                                        string residue2,
                                        Options & opts,
                                        ofstream& output_file,
                                        StatisticsRecorder *stats)
{
  // Due to the way centers are handled, I have to save the
  // original centers and replace them with the center of
  // carbon rings one at a time.

///*
  vector<Coordinates> original_center = ligand.center;
  vector< vector<Atom*> > original_altlocs = ligand.altlocs;
  vector< Atom* > original_atom;

  string original_residue = ligand.residue;

  ligand.center.clear();
  ligand.altlocs.clear();

  Chain* c1 = &(PDBfile.chains[chain1]);

  stats->setActiveLigand(ligand, ligand.getResidueNumber(), QPOLE);

  // Set skip flag on all ligand atoms
  for(int i=0; i < ligand.atom.size(); i++)
    {
      ligand.atom[i]->skip = true;
    }

  for(int i=0; i<c1->aa.size(); i++)
    {
      if(c1->aa[i].residue == residue2 && !(c1->aa[i].skip))
        {
          stats->setActiveResidue(c1->aa[i], c1->aa[i].getResidueNumber(), ANION);

          for (int j=0; j<ligand.carbonRingCenters.size(); j++)
            {
              // Ignore all the other carbon ring centers except for the jth.
              for (int k=0; k<ligand.carbonRingCenters[j].size(); k++)
                {
                  ligand.center.push_back(ligand.carbonRingCenters[j][k]);
                  ligand.center[k].skip = false;

                  vector<Atom*> temp;
                  for (int i=0; i < ligand.carbonRings[j][k].size(); i++)
                    {
                      temp.push_back(ligand.carbonRings[j][k][i]);
                    }
                  for (int i=0; i < ligand.additionalAtoms[j][k].size(); i++)
                    {
                      temp.push_back(ligand.additionalAtoms[j][k][i]);
                    }
                  for (int i=0; i < temp.size(); i++)
                    {
                      temp[i]->skip = false;	
                    }

                  ligand.altlocs.push_back(temp);
                }

              // Append ring number to residue name 
              ligand.residue = original_residue;
              stringstream ring_num;
              ring_num << j;
              ligand.residue = original_residue + "_" + ring_num.str();

              bool found;
              found = findBestInteraction(ligand,
                                          c1->aa[i],
                                          opts.threshold,
                                          PDBfile,
                                          opts.gamessfolder,
                                          false,
                                          output_file);
              if (found)
                {
                  stats->recordNewInteraction(0,0);
                }
              ligand.center.clear();
              ligand.altlocs.clear();

//              cout << "Done" << endl;
            }
        }
    }

  ligand.center  = original_center;
  ligand.altlocs = original_altlocs;
  ligand.residue = original_residue;

  return;
//*/
}

// Finds the closest distance among all of the centers
// associated with each amino acid
double findClosestDistance(Residue& aa1,
                           Residue& aa2,
                           float threshold,
                           unsigned int* closest_index1,
                           unsigned int* closest_index2)
{
  double dist;
  double closest = FLT_MAX;

  // Go through all combination of distances looking
  // for the closet pair
  //cout << "aa1.center.size() = " << aa1.center.size() << endl;
//  cout << "Residue: " << aa1.residue << ", aa1.altlocs.size() = " << aa1.altlocs.size() << endl;
  //cout << "aa2.center.size() = " << aa2.center.size() << endl;
//  cout << "Residue: " << aa2.residue << ", aa2.altlocs.size() = " << aa2.altlocs.size() << endl;
  for( unsigned int i = 0; i < aa1.center.size(); i++ )
    {
      if(aa1.center[i].skip) continue;
      for( unsigned int j = 0; j < aa2.center.size(); j++ )
        {
          if(aa2.center[j].skip) continue;
          dist = aa1.center[i].distance(aa2.center[j]);
          // flag it if is the closest and within the threshold
          if( dist < closest && dist < threshold )
            {
              *closest_index1 = i;
              *closest_index2 = j;
              closest = dist;
            }
        }
    }
  return closest;
}

bool findBestInteraction( Residue& aa1,
                          Residue& aa2,
                          float threshold,
                          PDB & PDBfile,
                          char* gamessfolder,
                          bool ligand,
                          ofstream& output_file)
{
  float closestDist = FLT_MAX;
  unsigned int closestDist_index1 = 0;
  unsigned int closestDist_index2 = 0;
  float dist;
  float angle;
  float angle1;
  float angleP;
  float angleh;
  float angleOxy;
  float angleOxy2;
  Residue aa1h;
  Residue aa2h;
  char output_filename[1024] = "N/A";
  static int numOutputted = 0;

  // Go through all combination of distances looking
  // for the closet pair
  closestDist = findClosestDistance(aa1,
                                    aa2,
                                    threshold,
                                    &closestDist_index1,
                                    &closestDist_index2);

  // AND WE HAVE A WINNER! 
  // if we found something output it to the file,
  // add hydrogens, and try the process again
  if( closestDist != FLT_MAX )
    {
      // Just some codes that were in the original STAAR
      char code1 = 'I';
      if( aa1.atom[0]->chainID != aa2.atom[0]->chainID )
        code1 = 'X';
      char code2 = 'S';
      //if( aa1.center.size() > 1 || aa2.center.size() > 1 )
      if( aa1.altLoc || aa2.altLoc )
        {
          code2 = 'M';
        }

      PDB pairWithHydrogen;

      // Set the residues and ligands to find
      pairWithHydrogen.setResiduesToFind(PDBfile.residue1, PDBfile.residue2);
      pairWithHydrogen.setLigandsToFind(PDBfile.ligandsToFind);
      

      // Add the hydrogens
      pairWithHydrogen.addHydrogensToPair(aa1,aa2,closestDist_index1,closestDist_index2);

      // Set the filename
      pairWithHydrogen.filename = PDBfile.filename;

      // Separate the pair into 2 variables
//      aa1h = aa1;
//      aa2h = aa2;

      pairWithHydrogen.getPair(aa1.atom[0]->resSeq,
                               aa2.atom[0]->resSeq,
                               &aa1h,
                               &aa2h,
                               ligand);


      if( aa2h.skip == true || aa1h.skip == true )
        {
          return false;
        }

      float dist;
      float distOxy;
      float distOxy2;
      if(!aa1h.calculateDistancesAndAnglesPostHydrogens(aa2h,
                                                        aa2.center[closestDist_index2],
                                                        threshold,
                                                        &dist,
                                                        &distOxy,
                                                        &distOxy2,
                                                        &angleh,
                                                        &angleOxy,
                                                        &angleOxy2))
        {
          // No interaction was found
          return false;
        }

      // calculate the angles of this interaction
      aa1.calculateAnglesPreHydrogens(aa2,
                                      closestDist_index1,
                                      closestDist_index2,
                                      &angle,
                                      &angle1,
                                      &angleP);
      // Construct the atom number strings
      stringstream res1_atoms;
      vector<Atom *> res1_atoms_vector = aa1.altlocs[0];//[closestDist_index1];
      sort(res1_atoms_vector.begin(), res1_atoms_vector.end());
      for (int i=0; i < res1_atoms_vector.size(); i++)
        {
          if (i > 0) res1_atoms << " ";
          res1_atoms << res1_atoms_vector[i]->serialNumber;
        }

      stringstream res2_atoms;
      vector<Atom *> res2_atoms_vector = aa2.altlocs[0];//[closestDist_index2];
      sort(res2_atoms_vector.begin(), res2_atoms_vector.end());
      for (int i=0; i < res2_atoms_vector.size(); i++)
        {
          if (i > 0) res2_atoms << " ";
          res2_atoms << res2_atoms_vector[i]->serialNumber;
        }

      // Get the number of carbon rings each residue has
      int res1_nrings = aa1.carbonRings.size();
      if (res1_nrings == 0 && aa1.residue == "PHE") res1_nrings = 1;

      int res2_nrings = aa2.carbonRings.size();
      if (res2_nrings == 0 && aa2.residue == "PHE") res2_nrings = 1;

      // The following is pretty hackish.  For the time being since we don't
      // have an agreement on how to deal with the PO4 ligands completely, 
      // we will take 1 H out at a time and output the GAMESS input file.
      // Thus, for each PO4, we will have 3 GAMESS input files.
      // Good thing there aren't too many of these
      if(aa2h.residue == "PO4")
        {
          vector<Atom*>::iterator atom_iterator;
          int count = 0;
          for( atom_iterator  = aa2h.atom.begin();
               atom_iterator != aa2h.atom.end();
               ++atom_iterator )
            {
              if((*atom_iterator)->name == " H  ")
                {
                  (*atom_iterator)->skip = true;
                  if( gamessfolder )
                    {
                      numOutputted++;
                      sprintf(output_filename, "%s/gamessinp-%d.inp", gamessfolder, numOutputted);
                      outputINPfile(PDBfile.filename, output_filename, aa1h, aa2h);
                    }
                  else
                    {
                      sprintf(output_filename, "N/A");
                    }

                  // and we finally output some results!
                  output_file << aa1.residue                        << ","
                              << aa2.residue                        << ","
                              << closestDist                        << ","
                              << angle                              << ","
                              << angleP                             << ","
                              << angle1                             << ","
                              << aa1.atom[0]->resSeq                << ","
                              << aa2.atom[0]->resSeq                << ","
                              << code1 << code2                     << ","
                              << PDBfile.filename                   << ","
                              << PDBfile.resolution                 << ","
                              << output_filename                    << ","
                              << aa1.atom[0]->chainID               << ","
                              << aa2.atom[0]->chainID               << ","
                              << aa1.center[closestDist_index1]     << ","
                              << aa2.center[closestDist_index2]     << ","
                              << aa1h.center[0]                     << ","
                              << aa2h.center[0]                     << ","
                              << dist                               << ","
                              << distOxy                            << ","
                              << distOxy2                           << ","
                              << angleh                             << ","
                              << angleOxy                           << ","
                              << angleOxy2                          << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << " "                                << ","
                              << res1_atoms.str()                   << ","
                              << res2_atoms.str()                   << ","
                              << res1_nrings                        << ","
                              << res2_nrings                        << endl;
                  (*atom_iterator)->skip = false;
                  if(count == 2) break;
                  ++count;
                }
            }
        }
      else 
        {
          // Set skip flags appropriately for carbon ring ligands
          if (aa1h.residue != "PHE")
            {
              for (int i=0; i<aa1h.atom.size(); i++)
                {
                  aa1h.atom[i]->skip = true;
                }
              for (int i=0; i<aa1h.altlocs[0].size(); i++)
                {
                  aa1h.altlocs[0][i]->skip = false;
                }
            }


          // Here, we are outputting the
          if( gamessfolder )
            {
              numOutputted++;
              sprintf(output_filename, "%s/gamessinp-%d.inp", gamessfolder, numOutputted);
              outputINPfile(PDBfile.filename, output_filename, aa1h, aa2h);
            }
          else
            {
              sprintf(output_filename, "N/A");
            }

          // and we finally output some results!
          output_file << aa1.residue                        << ","
                      << aa2.residue                        << ","
                      << closestDist                        << ","
                      << angle                              << ","
                      << angleP                             << ","
                      << angle1                             << ","
                      << aa1.atom[0]->resSeq                << ","
                      << aa2.atom[0]->resSeq                << ","
                      << code1 << code2                     << ","
                      << PDBfile.filename                   << ","
                      << PDBfile.resolution                 << ","
                      << PDBfile.model_number               << ","
                      << output_filename                    << ","
                      << aa1.atom[0]->chainID               << ","
                      << aa2.atom[0]->chainID               << ","
                      << aa1.center[closestDist_index1]     << ","
                      << aa2.center[closestDist_index2]     << ","
                      << aa1h.center[0]                     << ","
                      << aa2h.center[0]                     << ","
                      << dist                               << ","
                      << distOxy                            << ","
                      << distOxy2                           << ","
                      << angleh                             << ","
                      << angleOxy                           << ","
                      << angleOxy2                          << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << " "                                << ","
                      << res1_atoms.str()                   << ","
                      << res2_atoms.str()                   << ","
                      << res1_nrings                        << ","
                      << res2_nrings                        << endl;
        }

      return true;
    }
  else
    {
      return false;
    }
}

void outputINPfile(string input_filename, char* filename, Residue& aa1h, Residue& aa2h)
{
  ofstream inpout(filename);

  inpout << INPheader << endl;

  int qpole_mol_size = 0;
  for (int i=0; i < aa1h.atom.size(); i++)
    {
      if ( !aa1h.atom[i]->skip )
        {
          qpole_mol_size += 1;
        }
    }
//  inpout << " $MOROKM IATM(1)=" << aa1h.atom.size() << ",";
  inpout << " $MOROKM IATM(1)=" << qpole_mol_size << ",";
  if(aa2h.residue == "PO4")
    {
      inpout<< aa2h.atom.size() - 1 << " ICHM(1)=0,-1" << " $END" << endl;
    }
  else
    {
      int anion_mol_size = 0;
      for (int i=0; i < aa2h.atom.size(); i++)
        {
          if ( !aa2h.atom[i]->skip )
            {
              anion_mol_size += 1;
            }
        }
//      inpout<< aa2h.atom.size() << " ICHM(1)=0,-1" << " $END" << endl;
      inpout<< anion_mol_size << " ICHM(1)=0,-1" << " $END" << endl;
    }
  inpout << " $DATA" << endl;
  inpout << input_filename << endl;
  inpout << "C1" << endl;
  for(int i=0; i<aa1h.atom.size(); i++)
    {
      if ( !aa1h.atom[i]->skip ) 
        {
          if( aa1h.atom[i]->element == " H" )
            {
              inpout << "H      1.0     ";
            }
          else if( aa1h.atom[i]->element == " C")
            {
              inpout << "C      6.0     ";
            }
          else if( aa1h.atom[i]->element == " N")
            {
              inpout << "N      7.0     ";
            }
          else if( aa1h.atom[i]->element == " O")
            {
              inpout << "O      8.0     ";
            }
          else if( aa1h.atom[i]->element_num == F)
            {
              inpout << "F      9.0     ";
            }
          else if( aa1h.atom[i]->element_num == P)
            {
              inpout << "P     15.0     ";
            }
          else if( aa1h.atom[i]->element_num == S)
            {
              inpout << "S     16.0     ";
            }
          else if( aa1h.atom[i]->element_num == Cl)
            {
              inpout << "Cl    17.0     ";
            }
          else if( aa1h.atom[i]->element_num == Br)
            {
              inpout << "Br    35.0     ";
            }
          else
            {
              cout << "Warning!! Unhandled atom when creating GAMESS input file: "
                   << aa1h.atom[i]->element << endl;
            }
          inpout << aa1h.atom[i]->coord << endl;
        }
    }
  for(int i=0; i<aa2h.atom.size(); i++)
    {
      if( !aa2h.atom[i]->skip )
        {
          if( aa2h.atom[i]->element == " H")
            {
              inpout << "H      1.0     ";
            }
          else if( aa2h.atom[i]->element == " C")
            {
              inpout << "C      6.0     ";
            }
          else if( aa2h.atom[i]->element == " O")
            {
              inpout << "O      8.0     ";
            }
          else if( aa2h.atom[i]->element == " P")
            {
              inpout << "P     15.0     ";
            }
          inpout << aa2h.atom[i]->coord << endl;
        }
    }

  inpout << " $END" << endl;
  inpout.close();
}

void write_output_head(ofstream& out)
{
  //out <<"#res1,res2,dist,angle,angleP,angle1,loc1,loc2,code,pdbID,resolution,model,gamessinput,chain1,chain2,center1,,,center2,,,center1h,,,center2h,,,dist,distOxy,distOxy2,angleh,angleOxy,angleOxy2,gamessoutput,electrostatic(Hartree),electrostatic(kcal/mol),exchangerep(Hartree),exchangerep(kcal/mole),polarization(Hartree),polarization(kcal/mole),chargexfer(Hartree),chargexfer(kcal/mol),highordercoup(Hartree),highordercoup(kcal/mole),totalinter(Hartree),totalinter(kcal/mole)" << endl;
  out <<"#res1,res2,dist,angle,angleP,angle1,loc1,loc2,code,pdbID,resolution,model,gamessinput,chain1,chain2,center1,,,center2,,,center1h,,,center2h,,,dist,distOxy,distOxy2,angleh,angleOxy,angleOxy2,gamessoutput,electrostatic(Hartree),electrostatic(kcal/mol),exchangerep(Hartree),exchangerep(kcal/mole),polarization(Hartree),polarization(kcal/mole),chargexfer(Hartree),chargexfer(kcal/mol),highordercoup(Hartree),highordercoup(kcal/mole),totalinter(Hartree),totalinter(kcal/mole),res1Atoms,res2Atoms,res1_#rings,res2_#rings" << endl;
}

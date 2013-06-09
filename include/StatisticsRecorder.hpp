
#ifndef __STATISTICSRECORDER_HPP__
#define __STATISTICSRECORDER_HPP__

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Residue.hpp"

using namespace std;

enum MoleculeRole {ANION, QPOLE};
enum MoleculeType {RESIDUE, LIGAND};

class StatisticsRecorder
{
  private:
    class Record;
    class LigandRecord;
    class ResidueRecord;
    class ProteinRecord;
    class InteractionRecord;

    uint n_total_interactions;
    
    uint n_proteins;
    uint n_ligands;

    map<string, ProteinRecord *> proteins;
    map<string, ResidueRecord *> residues;
    map<string, LigandRecord  *> ligands;
    vector<InteractionRecord  *> interactions;


    // The protein that is currently being processed
    // in the main STAAR code.
    string active_protein;

    // The quadrupole that is currently being processed
    // in the main STAAR code.
    string active_qpole;
    // The type of molecule the active quadrupole is,
    // either RESIDUE or LIGAND.
    MoleculeType active_qpole_type;

    // The anion that is currently being processed
    // in the main STAAR code.
    string active_anion;
    // The type of molecule the active anion is,
    // either RESIDUE or LIGAND.
    MoleculeType active_anion_type;

    bool isNewProtein(string name);
    bool isNewResidue(string name);
    bool isNewLigand(string name);

  public:
    void setActiveProtein(string name);
    void setActiveModel(string name);
    void setActiveResidue(Residue &residue, uint residue_num, MoleculeRole role);
    void setActiveLigand(Residue &ligand, uint residue_num, MoleculeRole role);

    string getActiveProteinName();

    void recordNewInteraction(int qpole_center_idx, int anion_center_idx);

    int getTotalNumberOfInteractions();
    double averageInteractionsPerProtein();
    vector<int> carbonRingLigandCounts();
    vector<int> carbonRingLigandInteractionCounts();

};

#endif

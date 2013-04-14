
#ifndef __STATISTICSRECORDER_HPP__
#define __STATISTICSRECORDER_HPP__

#include <iostream>
#include <map>
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
    /*
    class Record {
      private:
      public:
        string name;

        bool operator==(Record &);

        
    }
    class LigandRecord : public Record {
        unsigned int n_qpoles;
        unsigned int n_anions;
    }

    class ResidueRecord : public Record {
        unsigned int n_qpoles;
        unsigned int n_anions;

    }

    class ProteinRecord : public Record {

    }
    */

    unsigned int n_total_interactions;
    
    unsigned int n_proteins;
    unsigned int n_ligands;

    map<string, ProteinRecord *> proteins;
    map<string, ResidueRecord *> residues;
    map<string, LigandRecord  *> ligands;
    vector<InteractionRecord  *> interactions;

    //ProteinRecord* active_protein;
    string active_protein;
//    LigandRecord * active_ligand;
//    ResidueRecord* active_residue;

    string active_qpole;
    MoleculeType active_qpole_type;

    string active_anion;
    MoleculeType active_anion_type;

    bool isNewProtein(string name);
    bool isNewResidue(string name);
    bool isNewLigand(string name);

  public:
    void setActiveProtein(string name);
    void setActiveModel(string name);
    void setActiveResidue(Residue &residue, MoleculeRole role);
    void setActiveLigand(Residue &ligand, MoleculeRole role);

    string getActiveProteinName();

    void recordNewInteraction(int qpole_center_idx, int anion_center_idx);

    int getTotalNumberOfInteractions();
    double averageInteractionsPerProtein();
    vector<int> carbonRingLigandCounts();

};

#endif

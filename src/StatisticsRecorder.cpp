
#include "StatisticsRecorder.hpp"

/*
  record:
    number of proteins
    number of ligands 
    number of unique ligands 
    number of ligands of interest (with a carbon ring or anion)
    number of unique ligands of interest 
    number of regions of interest per ligands 
    number of potential interactions (total, per residue)
    
*/

class StatisticsRecorder::Record
{
  private:

  public:
    string name;

    bool operator==(Record &rhs)
    {
      return (name == rhs.name);
    }

};

class StatisticsRecorder::LigandRecord : public StatisticsRecorder::Record
{
  private:
  public:
    uint n_qpoles;
    uint n_anions;

    // Set of protein names and residue numbers
    // that uniquely identify an instance of this
    // ligand.
    set< pair<string, uint> > instances;

    // Record a new instance. Returns whether or not this instance
    // has already been recorded.
    bool addInstance(string protein_name, uint residue_num)
    {
      return (instances.insert(make_pair(protein_name, residue_num))).second;
    }


    LigandRecord(Residue &ligand, string protein, uint residue_num)
    {
      // Strip a possible added ring number
      size_t pos = ligand.residue.find_last_of('_');
      name = ligand.residue.substr(0, pos-1);

      n_qpoles = ligand.carbonRings.size();
      n_anions = 0;
      if (n_qpoles == 0)
        n_anions = 1;

      addInstance(protein, residue_num);
    }

};

class StatisticsRecorder::ResidueRecord : public StatisticsRecorder::Record
{
  private:
  public:
    uint n_qpoles;
    uint n_anions;

    // Set of protein names and residue numbers
    // that uniquely identify an instance of this
    // residue.
    set< pair<string, uint> > instances;

    // Record a new instance. Returns whether or not this instance
    // has already been recorded.
    bool addInstance(string protein_name, uint residue_num)
    {
      return (instances.insert(make_pair(protein_name, residue_num))).second;
    }

    ResidueRecord(Residue &residue)
    {
      name = residue.residue;
      if (name == "PHE" || name == "TYR" || name == "TRP")
        {
          n_qpoles = 1;
          n_anions = 0;
        }
      else if (name == "GLU" || name == "ASP")
        {
          n_qpoles = 0;
          n_anions = 1;
        }
    }
};

class StatisticsRecorder::ProteinRecord : public StatisticsRecorder::Record
{
  private:
  public:
    ProteinRecord(string init_name) 
    {
      name = init_name;
    }
};

class StatisticsRecorder::InteractionRecord
{
  private:
  public:
    string protein_name;
    string qpole_name;
    string anion_name;

    MoleculeType qpole_type;
    MoleculeType anion_type;

    InteractionRecord(string p_name, string q_name, MoleculeType q_type,
                                     string a_name, MoleculeType a_type)
    {
      protein_name = p_name;
      qpole_name = q_name;
      qpole_type = q_type;
      anion_name = a_name;
      anion_type = a_type;
    }
};




bool StatisticsRecorder::isNewProtein(string name)
{
  if (proteins.empty())
    {
      return true;
    }

  map<string, ProteinRecord *>::iterator it = proteins.find(name);

  if (it == proteins.end())
    {
      return true;
    }
  else
    {
      return false;
    }
}

bool StatisticsRecorder::isNewResidue(string name)
{
  if (residues.empty())
    {
      return true;
    }

  map<string, ResidueRecord *>::iterator it = residues.find(name);

  if (it == residues.end())
    {
      return true;
    }
  else
    {
      return false;
    }
}

bool StatisticsRecorder::isNewLigand(string name)
{
  if (ligands.empty())
    {
      return true;
    }

  map<string, LigandRecord *>::iterator it = ligands.find(name);

  if (it == ligands.end())
    {
      return true;
    }
  else
    {
      return false;
    }

}


void StatisticsRecorder::setActiveProtein(string name)
{
  if (isNewProtein(name))
    {
      ProteinRecord *newProtein = new ProteinRecord(name);
      proteins.insert(pair<string, ProteinRecord *>(name, newProtein));
      active_protein = name;
    }
  else
    {
      active_protein = name;
    }
}

string StatisticsRecorder::getActiveProteinName()
{
  return active_protein;
}

void StatisticsRecorder::setActiveModel(string name)
{

}

void StatisticsRecorder::setActiveResidue(Residue &residue, uint residue_num, MoleculeRole role)
{
  string name = residue.residue;

  if (isNewResidue(name))
    {
      ResidueRecord *newResidue = new ResidueRecord(residue);
      residues.insert(pair<string, ResidueRecord *>(name, newResidue));
    }

  if (role == ANION)
    {
      active_anion = name;
      active_anion_type = RESIDUE;
    }
  else if (role == QPOLE)
    {
      active_qpole = name;
      active_qpole_type = RESIDUE;
    }
}

void StatisticsRecorder::setActiveLigand(Residue &ligand, uint residue_num, MoleculeRole role)
{
  // String a possible ring number at the end
  size_t pos = ligand.residue.find_last_of('_');
  string name = ligand.residue.substr(0, pos-1);

  if (isNewLigand(name))
    {
      cout << "NEW LIGAND: " << name << endl;
      LigandRecord *newLigand = new LigandRecord(ligand, active_protein, residue_num);
      ligands.insert(pair<string, LigandRecord *>(name, newLigand));
    }
  else
    {
      map<string, LigandRecord *>::iterator it1 = ligands.find(name);
      LigandRecord *ligand = it1->second;

      pair<string, uint> temp_pair = make_pair(active_protein, residue_num);
      set< pair<string, uint> >::iterator it2 = ligand->instances.find(temp_pair);
      if (it2 == ligand->instances.end()) // New instance
        {
          ligand->instances.insert(temp_pair);
        }
    }

  if (role == ANION)
    {
      active_anion = name;
      active_anion_type = LIGAND;
    }
  else if (role == QPOLE)
    {
      active_qpole = name;
      active_qpole_type = LIGAND;
    }
}

void StatisticsRecorder::recordNewInteraction(int qpole_center_idx, int anion_center_idx)
{
  n_total_interactions++;

  cerr << "New interaction between QPOLE " << active_qpole << " and ANION " << active_anion << endl;
  interactions.push_back(new InteractionRecord(active_protein, active_qpole, active_qpole_type,
                                                               active_anion, active_anion_type));
}

int StatisticsRecorder::getTotalNumberOfInteractions()
{
  return interactions.size();
}

double StatisticsRecorder::averageInteractionsPerProtein()
{
  return (double)interactions.size() / (double)proteins.size();
}


/**
  Counts how many instances of each ligand that contains carbon rings
   have been recorded. Counts are binned by the number of carbon rings.
**/
vector<int> StatisticsRecorder::carbonRingLigandCounts()
{
  vector<int> counts;


  map<string, LigandRecord *>::iterator it;

  for (it = ligands.begin(); it != ligands.end(); it++)
    {
      int n_rings = 0;
      LigandRecord *ligand = it->second;

      n_rings = ligand->n_qpoles;

      if (counts.size() < n_rings)
        {
          counts.resize(n_rings);
        }
      counts[n_rings-1] += ligand->instances.size();
    }

  return counts;
}

/**
  Counts how many interactions have been recorded that involve
   ligands with carbon rings. Counts are binned by the number of
   carbon rings.
  
**/
vector<int> StatisticsRecorder::carbonRingLigandInteractionCounts()
{
  vector<int> counts; 

  vector<InteractionRecord *>::iterator it;

  for (it = interactions.begin(); it != interactions.end(); it++)
    {
      int n_rings = 0;
      InteractionRecord *ir = *it;

      if (ir->qpole_type == LIGAND)
        {
          LigandRecord *ligand = ligands[ir->qpole_name];
          n_rings = ligand->n_qpoles;

          if (counts.size() < n_rings)
            {
              counts.resize(n_rings);
            }
          counts[n_rings-1]++;
        }
    }

  return counts;
}

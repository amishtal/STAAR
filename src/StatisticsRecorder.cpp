
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


/*
    StatisticsRecorder::Record::Record(string init_name)
    {
      name = name;
    }
*/

    bool operator==(Record &rhs)
    {
      return (name == rhs.name);
    }

};

class StatisticsRecorder::LigandRecord : public StatisticsRecorder::Record
{
  private:
  public:
    unsigned int n_qpoles;
    unsigned int n_anions;

    LigandRecord(Residue &ligand)
    {
      // Strip a possible added ring number
      size_t pos = ligand.residue.find_last_of('_');
      name = ligand.residue.substr(0, pos-1);

      n_qpoles = ligand.carbonRings.size();
      n_anions = 0;
      if (n_qpoles == 0)
        n_anions = 1;
    }

};

class StatisticsRecorder::ResidueRecord : public StatisticsRecorder::Record
{
  private:
  public:
    unsigned int n_qpoles;
    unsigned int n_anions;

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

void StatisticsRecorder::setActiveResidue(Residue &residue, MoleculeRole role)
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

void StatisticsRecorder::setActiveLigand(Residue &ligand, MoleculeRole role)
{
  // String a possible ring number at the end
  size_t pos = ligand.residue.find_last_of('_');
  string name = ligand.residue.substr(0, pos-1);

  if (isNewLigand(name))
    {
      cout << "NEW LIGAND: " << name << endl;
      LigandRecord *newLigand = new LigandRecord(ligand);
      ligands.insert(pair<string, LigandRecord *>(name, newLigand));
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

  // Increment a global counter?
  // Increment a counter for each ligand/residue/protein record?
  // Record each interaction? What information to keep?
}

int StatisticsRecorder::getTotalNumberOfInteractions()
{
  return interactions.size();
}

double StatisticsRecorder::averageInteractionsPerProtein()
{
  return (double)interactions.size() / (double)proteins.size();
}

vector<int> StatisticsRecorder::carbonRingLigandCounts()
{
  vector<int> counts; // <n_rings, count>

  vector<InteractionRecord *>::iterator it;

  for (it = interactions.begin(); it != interactions.end(); it++)
    {
      string ligand_name;
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

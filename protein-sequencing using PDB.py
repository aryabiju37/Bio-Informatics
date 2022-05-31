# -*- coding: utf-8 -*-

from Bio.PDB.MMCIFParser import MMCIFParser  # Tip: This module might be useful for parsing...
import numpy as np


############# Exercise 2: Protein Data Bank #############
# General remark: In our exercise every structure will have EXACTLY ONE model.
# This is true for nearly all X-Ray structures. NMR structures have several models.
class PDB_Parser:
    CIF_PARSER = MMCIFParser()  # parser object for reading in structure in CIF format

    def __init__(self, path):
        """
            Initialize every PDB_Parser with a path to a structure-file in CIF format.
            An example file is included in the repository (7ahl.cif).
            Tip: Store the parsed structure in an object variable instead of parsing it
            again & again ...

            """

        self.structure = self.CIF_PARSER.get_structure("PDB_id", path)

        # Parse the structure once and re-use it in the functions below


    # 2.8 Chains    
    def get_number_of_chains(self):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
            Return:
                Number of chains in this structure as integer.
        """

        n_chains = 0
        chains = []
        for chain in self.structure.get_chains():
            chains.append(chain)
        n_chains = len(chains)
        return n_chains

    # 2.9 Sequence  
    def get_sequence(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the amino acid sequence (single-letter alphabet!) of a given chain (chain_id)
                in a Biopython.PDB structure as a string.
        """

        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        sequences = ""
        for model in self.structure.get_list():
            chains = model[chain_id]

        sequences = ""
        aa_sequence = ""
        lst_residue = self.remove_hetero_water(chain_id)
        for residue in lst_residue:
            sequences+=residue.get_resname()
        for i in range(0, len(sequences), 3):
            if (i < len(sequences) - 2):
                if(sequences[i] + sequences[i + 1] + sequences[i + 2] in d):
                    aa_sequence += d[sequences[i] + sequences[i + 1] + sequences[i + 2]]

        return aa_sequence

    # 2.10 Water molecules
    def get_number_of_water_molecules(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the number of water molecules of a given chain (chain_id)
                in a Biopython.PDB structure as an integer.
        """

        n_waters = 0
        sequences = ""
        for model in self.structure:
            chains = model[chain_id]

        for residue in chains:
            sequences += residue.get_resname()

        for i in range(0, len(sequences), 3):
            if (i < len(sequences) - 2):
                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("HOH", "HHO", "OHH")):
                    n_waters += 1

        return n_waters

    # 2.11 C-alpha distance
    def get_ca_distance(self, chain_id_1, index_1, chain_id_2, index_2):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id_1 : String (usually in ['A','B', 'C' ...]. The number of chains
                                depends on the specific protein and the resulting structure)
                index_1    : index of a residue in a given chain in a Biopython.PDB structure
                chain_id_2 : String (usually in ['A','B', 'C' ...]. The number of chains
                            depends on the specific protein and the resulting structure)
                index_2    : index of a residue in a given chain in a Biopython.PDB structure

                chain_id_1 and index_1 describe precisely one residue in a PDB structure,
                chain_id_2 and index_2 describe the second residue.

            Return:
                Return the C-alpha (!) distance between the two residues, described by
                chain_id_1/index_1 and chain_id_2/index_2. Round the returned value via int().

            The reason for using two different chains as an input is that also the distance
            between residues of different chains can be interesting.
            Different chains in a PDB structure can either occur between two different proteins
            (Heterodimers) or between different copies of the same protein (Homodimers).
        """
        ca_distance = 0
        for model in self.structure:
            chain1 = model[chain_id_1]
            chain2 = model[chain_id_2]
            residue1 = chain1[index_1]
            residue2 = chain2[index_2]

        dist1 = residue1["CA"].get_coord()
        dist2 = residue2["CA"].get_coord()
        ca_distance = np.linalg.norm(dist1 - dist2)

        return int(ca_distance)

    # 2.12 Contact Map  
    def get_contact_map(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return a complete contact map (see description in exercise sheet)
                for a given chain in a Biopython.PDB structure as numpy array.
                The values in the matrix describe the c-alpha distance between all residues
                in a chain of a Biopython.PDB structure.
                Only integer values of the distance have to be given (see below).
        """

        length = 0
        for model in self.structure:
            chains = model[chain_id]
        res = chains.get_residues()
        residues = []
        for i in iter(res):
            residues.append(i)

        residue = []
        for res in residues:
            try:
                residue.append(res["CA"])
            except KeyError:
                continue

        length = len(residue)
        contact_map = np.zeros((length, length), dtype=np.float32)

        for i in range(len(residue)):
            for j in range(len(residue)):
                try:
                    contact_map[i][j] = np.linalg.norm(residues[i]["CA"].get_coord() - residues[j]["CA"].get_coord())
                except KeyError:
                    continue

        return contact_map.astype(np.int64)  # return rounded (integer) values

    def remove_hetero_water(self,chain_id):

        for model in self.structure:
            chains = model[chain_id]

        #remove hetero atoms from the residue
        for residue in chains:
            if (residue.id[0] != ' '):
                chains.detach_child(residue.id)

        residues_lst = []
        for residue in chains:
            residues_lst.append(residue)

        residue_lst = []
        for i in range(len(residues_lst)):
            if(residues_lst[i].id[0]!="W"):
                residue_lst.append(residues_lst[i])
        return residue_lst

    def bfactor_per_residue(self,residue):
        # for atom in residue_lst[0]:
        #     atom_bfactor_lst.append(atom.get_bfactor())
        # arr = np.array(atom_bfactor_lst)
        # mean = np.mean(arr)

        b_factor_per_atom_lst = []
        for atom in residue:
            b_factor_per_atom_lst.append(atom.get_bfactor())
        bfactor = np.array(b_factor_per_atom_lst)
        mean = np.mean(bfactor)
        return mean


    # 2.13 B-Factors    
    def get_bfactors(self, chain_id):
        """
            Input:
                self: Use Biopython.PDB structure which has been stored in an object variable
                chain_id  : String (usually in ['A','B', 'C' ...]. The number of chains
                        depends on the specific protein and the resulting structure)
            Return:
                Return the B-Factors for all residues in a chain of a Biopython.PDB structure.
                The B-Factors describe the mobility of an atom or a residue.
                In a Biopython.PDB structure B-Factors are given for each atom in a residue.
                Calculate the mean B-Factor for a residue by averaging over the B-Factor
                of all atoms in a residue.
                Sometimes B-Factors are not available for a certain residue;
                (e.g. the residue was not resolved); insert np.nan for those cases.

                Finally normalize your B-Factors using Standard scores (zero mean, unit variance).
                You have to use np.nanmean, np.nanvar etc. if you have nan values in your array.
                The returned data structure has to be a numpy array rounded again to integer.
        """

        residue_lst = self.remove_hetero_water(chain_id)
        b_factor_lst = []
        atom_bfactor_lst = []

        for i in range(len(residue_lst)):
            atom_bfactor_lst.append(self.bfactor_per_residue(residue_lst[i]))

        b_factors = np.array(atom_bfactor_lst)
        mean = np.mean(b_factors)
        dev = np.std(b_factors)
        b_factors = (b_factors - mean)/dev

        return b_factors.astype(np.int64)  # return rounded (integer) values


def main():
    print('PDB parser class.')
    return None


if __name__ == '__main__':
    main()

# pathh = "C:/Users/ARYA/Desktop/Studies/SS2022/PP1/EX2/pp1cs22exercise2-ge35saf-main@3e89e4a00fb/tests/7ahl.cif"
# pdb = PDB_Parser(pathh)


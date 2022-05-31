# -*- coding: utf-8 -*-
from Bio import SeqIO  # Tip: This module might be useful for parsing...


############ Exercise 3: SwissProt ##########
class SwissProt_Parser:
    PARSER = SeqIO
    dict = {}
    sequence = ""
    dbx = []

    def __init__(self, path, frmt='uniprot-xml'):
        """
            Initialize every SwissProt_Parser with a path to a XML-formatted UniProt file.
            An example file is included in the repository (P09616.xml).
            Tip: Store the parsed XML entry in an object variable instead of parsing it
            again & again ...
        """
        self.sp_anno = SeqIO.parse(path, frmt)  # Parse the XML file once and re-use it in the functions below

    # 2.2 SwissProt Identification
    def get_dictionary(self):
        global dict
        global sequence
        global dbx
        for index, record in enumerate(self.sp_anno):
            dict = record.annotations
            sequence = record.seq
            dbx = record.dbxrefs
        return dict,sequence,dbx


    def get_sp_identification(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                identification: tuple consisting of in this order
                                    1. string: Unique SwissProt identifier for the given xml file
                                    2. string: Primary protein name
                                    3. string: Primary gene name
        """

        return_value = self.get_dictionary()
        dictionary = return_value[0]
        identifier = (dictionary["accessions"][0],dictionary["recommendedName_fullName"][0],dictionary["gene_name_primary"])
        return identifier


    # 2.3 SwissProt Sequence Information
    def get_sp_sequence_info(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                information: tuple consisting of in this order
                                    1. str: sequence of the UniProt entry
                                    2. int: sequence length of the UniProt entry
                                    3. int: sequence mass of the UniProt entry


                                    """

        rec = self.get_dictionary()
        rec_dict = rec[0]
        sequence = str(rec[1])
        sequence_info = (sequence,len(sequence),rec_dict["sequence_mass"])

        return sequence_info



    # 2.4 Organism 
    def get_organism(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the organsim as stated in the corresponding field
                of the XML data. Return value has to be a string.
        """


        dic = self.get_dictionary()[0]
        return dic["organism"]


    # 2.5 Localizations
    def get_localization(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Return the name of the subcellular localization as stated in the
                corresponding field.
                Return value has to be a list of strings.
        """


        d = self.get_dictionary()[0]
        return  d["comment_subcellularlocation_location"]



    # 2.6 Cross-references to PDB
    def get_pdb_support(self):
        """
            Input:
                self: Use XML entry which has been parsed & saved during object initialization
            Return:
                Returns a list of all PDB IDs which support the annotation of the
                given SwissProt XML file. Return the PDB IDs as list.
        """

        pdb_ids = []
        dbxrefs = self.get_dictionary()[2]
        for i in range(len(dbxrefs)):
            if(dbxrefs[i].startswith("PDB")):
                 pdb_ids.append(dbxrefs[i])

        pdb_values = []
        for pdb in pdb_ids:
            pdb_values.append(pdb.split(":")[1])

        pdb_ids = list(dict.fromkeys(pdb_values))
        return pdb_ids




def main():
    print('SwissProt XML Parser class')
    return None


if __name__ == '__main__':
    main()

#
# pathh = "C:/Users/ARYA/Desktop/Studies/SS2022/PP1/EX2/pp1cs22exercise2-ge35saf-main@3e89e4a00fb/tests/P09616.xml"
# sp = SwissProt_Parser(pathh)
# #
# identification = sp.get_sp_identification()
# print(identification)
#
# sequence = sp.get_sp_sequence_info()
# print(sequence)
#
# organism = sp.get_organism()
# print(organism)
#
# localization = sp.get_localization()
# print(localization)
#
# pdb = sp.get_pdb_support()
# print(pdb)



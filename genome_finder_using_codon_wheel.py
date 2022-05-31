class Genome:

    def __init__(self, genome):
        """
        Initialize the Genome class with the provided genome sequence.

        :param genome: String with the genome sequence.
        """
        self.genome = genome

    def get_at_content(self):
        """
        Return the AT content of the genome sequence, i.e. the combined
        fraction of 'A' and 'T' in the entire genome sequence.

        :return: AT content (float, rounded to 6 digits)
        """

        A_counter = 0
        T_counter = 0
        for i in range(0, len(self.genome)):
            if (self.genome[i] == "A"):
                A_counter += 1

        for i in range(0, len(self.genome)):
            if (self.genome[i] == "T"):
                T_counter += 1

        AT_counter = A_counter / len(self.genome) + T_counter / len(self.genome)
        return round(AT_counter, 6)

    def A_freq(self,genome):
        A_counter = 0
        for i in range(len(self.genome)):
            if (self.genome[i] == "A"):
                A_counter += 1
        return A_counter / len(self.genome)

    def T_freq(self,genome):
        T_counter = 0
        for i in range(len(self.genome)):
            if (self.genome[i] == "T"):
                T_counter += 1
        return T_counter / len(self.genome)

    def G_freq(self,genome):
        G_counter = 0
        for i in range(len(self.genome)):
            if (self.genome[i] == "G"):
                G_counter += 1
        return G_counter / len(self.genome)

    def C_freq(self,genome):
        C_counter = 0
        for i in range(len(self.genome)):
            if (self.genome[i] == "C"):
                C_counter += 1
        return C_counter / len(self.genome)

    def get_codon_frequency(self,nt1, nt2, nt3, A_frequency, T_frequency, G_frequency, C_frequency):
        if (nt1 + nt2 + nt3 == "AAA"):
            return round(A_frequency * A_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "AAT"):
            return round(A_frequency * A_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "AAG"):
            return round(A_frequency * A_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "AAC"):
            return round(A_frequency * A_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "ATA"):
            return round(A_frequency * T_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "ATT"):
            return round(A_frequency * T_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "ATG"):
            return round(A_frequency * T_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "ATC"):
            return round(A_frequency * T_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "AGA"):
            return round(A_frequency * G_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "AGT"):
            return round(A_frequency * G_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "AGG"):
            return round(A_frequency * G_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "AGC"):
            return round(A_frequency * G_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "ACA"):
            return round(A_frequency * C_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "ACT"):
            return round(A_frequency * C_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "ACG"):
            return round(A_frequency * C_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "ACC"):
            return round(A_frequency * C_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "TAA"):
            return round(T_frequency * A_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "TAT"):
            return round(T_frequency * A_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "TAG"):
            return round(T_frequency * A_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "TAC"):
            return round(T_frequency * A_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "TTA"):
            return round(T_frequency * T_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "TTT"):
            return round(T_frequency * T_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "TTG"):
            return round(T_frequency * T_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "TTC"):
            return round(T_frequency * T_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "TGA"):
            return round(T_frequency * G_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "TGT"):
            return round(T_frequency * G_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "TGG"):
            return round(T_frequency * G_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "TGC"):
            return round(T_frequency * G_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "TCA"):
            return round(T_frequency * C_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "TCT"):
            return round(T_frequency * C_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "TCG"):
            return round(T_frequency * C_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "TCC"):
            return round(T_frequency * C_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "GAA"):
            return round(G_frequency * A_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "GAT"):
            return round(G_frequency * A_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "GAG"):
            return round(G_frequency * A_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "GAC"):
            return round(G_frequency * A_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "GTA"):
            return round(G_frequency * T_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "GTT"):
            return round(G_frequency * T_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "GTG"):
            return round(G_frequency * T_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "GTC"):
            return round(G_frequency * T_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "GGA"):
            return round(G_frequency * G_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "GGT"):
            return round(G_frequency * G_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "GGG"):
            return round(G_frequency * G_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "GGC"):
            return round(G_frequency * G_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "GCA"):
            return round(G_frequency * C_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "GCT"):
            return round(G_frequency * C_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "GCG"):
            return round(G_frequency * C_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "GCC"):
            return round(G_frequency * C_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "CAA"):
            return round(C_frequency * A_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "CAT"):
            return round(C_frequency * A_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "CAG"):
            return round(C_frequency * A_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "CAC"):
            return round(C_frequency * A_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "CTA"):
            return round(C_frequency * T_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "CTT"):
            return round(C_frequency * T_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "CTG"):
            return round(C_frequency * T_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "CTC"):
            return round(C_frequency * T_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "CGA"):
            return round(C_frequency * G_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "CGT"):
            return round(C_frequency * G_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "CGG"):
            return round(C_frequency * G_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "CGC"):
            return round(C_frequency * G_frequency * C_frequency, 6)

        if (nt1 + nt2 + nt3 == "CCA"):
            return round(C_frequency * C_frequency * A_frequency, 6)
        if (nt1 + nt2 + nt3 == "CCT"):
            return round(C_frequency * C_frequency * T_frequency, 6)
        if (nt1 + nt2 + nt3 == "CCG"):
            return round(C_frequency * C_frequency * G_frequency, 6)
        if (nt1 + nt2 + nt3 == "CCC"):
            return round(C_frequency * C_frequency * C_frequency, 6)



    def get_codon_dist(self):
        """
        Return the expected codon distribution (fractions) based on the
        distribution (fractions) of the four different nucleotides (ATGC).

        :return: Tree-like structure made out of nested dictionaries. The nodes
                 represent the different nucleotides and the path down the tree
                 forms the corresponding codons. The leafs contain the expected
                 codon frequencies (rounded to 6 digits).
        """

        A_frequency = self.A_freq(self.genome)
        T_frequency = self.T_freq(self.genome)
        G_frequency = self.G_freq(self.genome)
        C_frequency = self.C_freq(self.genome)

        codon_dict = {
            "A": {
                "A": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "T": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "G": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "C": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                }
            },
            "T": {
                "A": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "T": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "G": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "C": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                }
            },
            "G": {
                "A": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "T": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "G": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "C": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                }
            },
            "C": {
                "A": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "T": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "G": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                },
                "C": {
                    "A": 0.0,
                    "T": 0.0,
                    "G": 0.0,
                    "C": 0.0
                }
            }
        }

        for i in range(0, len(self.genome), 3):
            if (i < len(self.genome) - 2):
                codon_dict[self.genome[i]][self.genome[i + 1]][self.genome[i + 2]] = self.get_codon_frequency(self.genome[i],self.genome[i + 1],self.genome[i + 2],A_frequency,T_frequency,G_frequency,C_frequency)
        for i in range(1, len(self.genome), 3):
            if (i < len(self.genome) - 2):
                codon_dict[self.genome[i]][self.genome[i + 1]][self.genome[i + 2]] = self.get_codon_frequency(self.genome[i],self.genome[i + 1],self.genome[i + 2],A_frequency,T_frequency,G_frequency,C_frequency)
        for i in range(2, len(self.genome), 3):
            if (i < len(self.genome) - 2):
                codon_dict[self.genome[i]][self.genome[i + 1]][self.genome[i + 2]] = self.get_codon_frequency(self.genome[i],self.genome[i + 1],self.genome[i + 2],A_frequency,T_frequency, G_frequency,C_frequency)

        return codon_dict

    def get_codon_wheel(self):
        codon_wheel = {
            "A": {
                "A": {
                    "A": "K",
                    "T": "N",
                    "G": "K",
                    "C": "N"
                },
                "T": {
                    "A": "I",
                    "T": "I",
                    "G": "M",
                    "C": "I"
                },
                "G": {
                    "A": "R",
                    "T": "S",
                    "G": "R",
                    "C": "S"
                },
                "C": {
                    "A": "T",
                    "T": "T",
                    "G": "T",
                    "C": "T"
                }
            },
            "T": {
                "A": {
                    "A": "TAA",
                    "T": "Y",
                    "G": "TAG",
                    "C": "Y"
                },
                "T": {
                    "A": "L",
                    "T": "F",
                    "G": "L",
                    "C": "F"
                },
                "G": {
                    "A": "TGA",
                    "T": "C",
                    "G": "W",
                    "C": "C"
                },
                "C": {
                    "A": "S",
                    "T": "S",
                    "G": "S",
                    "C": "S"
                }
            },
            "G": {
                "A": {
                    "A": "E",
                    "T": "D",
                    "G": "E",
                    "C": "D"
                },
                "T": {
                    "A": "V",
                    "T": "V",
                    "G": "V",
                    "C": "V"
                },
                "G": {
                    "A": "G",
                    "T": "G",
                    "G": "G",
                    "C": "G"
                },
                "C": {
                    "A": "A",
                    "T": "A",
                    "G": "A",
                    "C": "A"
                }
            },
            "C": {
                "A": {
                    "A": "Q",
                    "T": "H",
                    "G": "Q",
                    "C": "H"
                },
                "T": {
                    "A": "L",
                    "T": "L",
                    "G": "L",
                    "C": "L"
                },
                "G": {
                    "A": "R",
                    "T": "R",
                    "G": "R",
                    "C": "R"
                },
                "C": {
                    "A": "P",
                    "T": "P",
                    "G": "P",
                    "C": "P"
                }
            }
        }
        return codon_wheel

    def amino_acid_freq(self,sequence1, sequence2, sequence3, codon_dict):
        if (sequence1 + sequence2 + sequence3 not in ("TAA", "TAG", "TGA")):
            numerator = codon_dict[sequence1][sequence2][sequence3]
            denomenator = 1 - (codon_dict["T"]["A"]["A"] + codon_dict["T"]["A"]["G"] + codon_dict["T"]["G"]["A"])
            res = numerator / denomenator
            return res
        else:
            return codon_dict[sequence1][sequence2][sequence3]





    def calculate_amino_freq(self,sequence, codon_wheel, codon_dict):
        amino_acid_dist = {
            "L": 0.0,
            "S": 0.0,
            "R": 0.0,
            "V": 0.0,
            "P": 0.0,
            "T": 0.0,
            "A": 0.0,
            "G": 0.0,
            "I": 0.0,
            "F": 0.0,
            "Y": 0.0,
            "H": 0.0,
            "Q": 0.0,
            "N": 0.0,
            "K": 0.0,
            "D": 0.0,
            "E": 0.0,
            "C": 0.0,
            "M": 0.0,
            "W": 0.0
        }
        TAC_flag = False
        TAT_flag = False
        TGT_flag = False
        TGC_flag = False
        TGG_flag = False
        CAG_flag = False
        CAA_flag = False
        AAC_flag = False
        AAT_flag = False
        GAT_flag = False
        GAC_flag = False
        GAA_flag = False
        GAG_flag = False
        TTT_flag = False
        TTC_flag = False
        AAA_flag = False
        AAG_flag = False
        CAC_flag = False
        CAT_flag = False
        GGA_flag = False
        GGT_flag = False
        GGG_flag = False
        GGC_flag = False
        GCA_flag = False
        GCT_flag = False
        GCG_flag = False
        GCC_flag = False
        GTA_flag = False
        GTC_flag = False
        GTG_flag = False
        GTT_flag = False
        CCA_flag = False
        CCT_flag = False
        CCG_flag = False
        CCC_flag = False
        ACA_flag = False
        ACG_flag = False
        ACT_flag = False
        ACC_flag = False
        ATA_flag = False
        ATT_flag = False
        ATG_flag = False
        ATC_flag = False
        TCA_flag = False
        TCT_flag = False
        TCG_flag = False
        TCC_flag = False
        AGC_flag = False
        AGT_flag = False
        AGA_flag = False
        AGG_flag = False
        CGA_flag = False
        CGT_flag = False
        CGG_flag = False
        CGC_flag = False
        TTA_flag = False
        TTG_flag = False
        CTA_flag = False
        CTG_flag = False
        CTT_flag = False
        CTC_flag = False
        Y = 0.0
        C = 0.0
        W = 0.0
        M = 0.0
        Q = 0.0
        N = 0.0
        D = 0.0
        E = 0.0
        F = 0.0
        K = 0.0
        H = 0.0
        G = 0.0
        A = 0.0
        V = 0.0
        P = 0.0
        T = 0.0
        I = 0.0
        S = 0.0
        R = 0.0
        L = 0.0
        sequences = sequence
        for i in range(0, len(sequences), 3):
            if (i < len(sequences) - 2):
                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TAC", "TAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TAC" and TAC_flag == False):
                        TAC_flag = True
                        Y += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TAT" and TAT_flag == False):
                        TAT_flag = True
                        Y += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(TAC_flag == True and TAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = Y

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TGT", "TGC")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGT" and TGT_flag == False):
                        TGT_flag = True
                        C += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGC" and TGC_flag == False):
                        TGC_flag = True
                        C += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(TGC_flag == True and TGT_flag==True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = C

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TGG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGG" and TGG_flag == False):
                        TGG_flag = True
                        W += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)
                    if(TGG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = W

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ATG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATG" and ATG_flag == False):
                        ATG_flag = True
                        M+= self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(ATG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = M

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CAG", "CAA")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAG" and CAG_flag == False):
                        CAG_flag = True
                        Q += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAA" and CAA_flag == False):
                        CAA_flag = True
                        Q += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(CAG_flag==True and CAA_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = Q

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AAC", "AAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAC" and AAC_flag == False):
                        AAC_flag = True
                        N += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAT" and AAT_flag == False):
                        AAT_flag = True
                        N += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(AAC_flag == True and AAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = N

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GAT", "GAC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAT" and GAT_flag == False):
                        GAT_flag = True
                        D += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAC" and GAC_flag == False):
                        GAC_flag = True
                        D += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(GAT_flag == True and GAC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = D

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GAA", "GAG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAA" and GAA_flag == False):
                        GAA_flag = True
                        E += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAG" and GAG_flag == False):
                        GAG_flag = True
                        E += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(GAA_flag == True and GAG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = E

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TTT", "TTC")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTT" and TTT_flag == False):
                        TTT_flag = True
                        F += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTC" and TTC_flag == False):
                        TTC_flag = True
                        F += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(TTT_flag == True and TTC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = F

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AAA", "AAG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAA" and AAA_flag == False):
                        AAA_flag = True
                        K += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAG" and AAG_flag == False):
                        AAG_flag = True
                        K += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(AAA_flag == True and AAG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = K

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CAC", "CAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAC" and CAC_flag == False):
                        CAC_flag = True
                        H += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAT" and CAT_flag == False):
                        CAT_flag = True
                        H += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(CAC_flag == True and CAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = H

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GGA", "GGT", "GGG", "GGC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGA" and GGA_flag == False):
                        GGA_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGT" and GGT_flag == False):
                        GGT_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGG" and GGG_flag == False):
                        GGG_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGC" and GGC_flag == False):
                        GGC_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(GGA_flag == True and GGT_flag == True and GGG_flag== True and GGC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = G

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GCA", "GCT", "GCG", "GCC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCA" and GCA_flag == False):
                        GCA_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCT" and GCT_flag == False):
                        GCT_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCG" and GCG_flag == False):
                        GCG_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCC" and GCC_flag == False):
                        GCC_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(GCA_flag == True and GCT_flag  ==True and GCG_flag == True and GCC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = A

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GTA", "GTT", "GTG", "GTC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTA" and GTA_flag == False):
                        GTA_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTT" and GTT_flag == False):
                        GTT_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTG" and GTG_flag == False):
                        GTG_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTC" and GTC_flag == False):
                        GTC_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(GTA_flag == True and GTT_flag == True and GTC_flag==True and GTG_flag==True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]]=V

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CCA", "CCT", "CCG", "CCC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCA" and CCA_flag == False):
                        CCA_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCT" and CCT_flag == False):
                        CCT_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCG" and CCG_flag == False):
                        CCG_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCC" and CCC_flag == False):
                        CCC_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(CCA_flag== True and CCT_flag == True and CCG_flag == True and CCC_flag==True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = P

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ACA", "ACT", "ACG", "ACC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACA" and ACA_flag == False):
                        ACA_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACT" and ACT_flag == False):
                        ACT_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACC" and ACC_flag == False):
                        ACC_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACG" and ACG_flag == False):
                        ACG_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(ACA_flag == True and ACT_flag == True and ACG_flag==True and ACC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = T

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ATA", "ATT", "ATC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATA" and ATA_flag == False):
                        ATA_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATT" and ATT_flag == False):
                        ATT_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATC" and ATC_flag == False):
                        ATC_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)


                    if(ATA_flag==True and ATT_flag==True and ATC_flag==True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = I

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TCA", "TCT", "TCG", "TCC", "AGC", "AGT")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCA" and TCA_flag == False):
                        TCA_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCT" and TCT_flag == False):
                        TCT_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCG" and TCG_flag == False):
                        TCG_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCC" and TCC_flag == False):
                        TCC_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGC" and AGC_flag == False):
                        AGC_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGT" and AGT_flag == False):
                        AGT_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(TCA_flag == True and TCT_flag == True and TCC_flag==True and TCG_flag==True and AGC_flag==True and AGT_flag==True ):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = S

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AGA", "AGG", "CGA", "CGT", "CGG", "CGC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGA" and AGA_flag == False):
                        AGA_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGG" and AGG_flag == False):
                        AGG_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGA" and CGA_flag == False):
                        CGA_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGT" and CGT_flag == False):
                        CGT_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGG" and CGG_flag == False):
                        CGG_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGC" and CGC_flag == False):
                        CGC_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(AGA_flag == True and AGG_flag==True and CGA_flag == True and CGT_flag==True and CGC_flag == True and CGG_flag==True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = R

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TTA", "TTG", "CTA", "CTT", "CTG", "CTC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTA" and TTA_flag == False):
                        TTA_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTG" and TTG_flag == False):
                        TTG_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTA" and CTA_flag == False):
                        CTA_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTT" and CTT_flag == False):
                        CTT_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTG" and CTG_flag == False):
                        CTG_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTC" and CTC_flag == False):
                        CTC_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if(TTA_flag == True and TTG_flag==True and CTA_flag == True and CTG_flag == True and CTC_flag == True and CTT_flag==True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = L
        for i in range(1, len(sequences), 3):
            if (i < len(sequences) - 2):
                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TAC", "TAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TAC" and TAC_flag == False):
                        TAC_flag = True
                        Y += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TAT" and TAT_flag == False):
                        TAT_flag = True
                        Y += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TAC_flag == True and TAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = Y

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TGT", "TGC")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGT" and TGT_flag == False):
                        TGT_flag = True
                        C += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGC" and TGC_flag == False):
                        TGC_flag = True
                        C += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TGC_flag == True and TGT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = C

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TGG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGG" and TGG_flag == False):
                        TGG_flag = True
                        W += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)
                    if (TGG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = W

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ATG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATG" and ATG_flag == False):
                        ATG_flag = True
                        M += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (ATG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = M

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CAG", "CAA")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAG" and CAG_flag == False):
                        CAG_flag = True
                        Q += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAA" and CAA_flag == False):
                        CAA_flag = True
                        Q += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (CAG_flag == True and CAA_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = Q

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AAC", "AAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAC" and AAC_flag == False):
                        AAC_flag = True
                        N += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAT" and AAT_flag == False):
                        AAT_flag = True
                        N += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (AAC_flag == True and AAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = N

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GAT", "GAC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAT" and GAT_flag == False):
                        GAT_flag = True
                        D += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAC" and GAC_flag == False):
                        GAC_flag = True
                        D += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GAT_flag == True and GAC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = D

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GAA", "GAG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAA" and GAA_flag == False):
                        GAA_flag = True
                        E += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAG" and GAG_flag == False):
                        GAG_flag = True
                        E += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GAA_flag == True and GAG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = E

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TTT", "TTC")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTT" and TTT_flag == False):
                        TTT_flag = True
                        F += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTC" and TTC_flag == False):
                        TTC_flag = True
                        F += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TTT_flag == True and TTC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = F

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AAA", "AAG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAA" and AAA_flag == False):
                        AAA_flag = True
                        K += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAG" and AAG_flag == False):
                        AAG_flag = True
                        K += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (AAA_flag == True and AAG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = K

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CAC", "CAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAC" and CAC_flag == False):
                        CAC_flag = True
                        H += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAT" and CAT_flag == False):
                        CAT_flag = True
                        H += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (CAC_flag == True and CAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = H

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GGA", "GGT", "GGG", "GGC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGA" and GGA_flag == False):
                        GGA_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGT" and GGT_flag == False):
                        GGT_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGG" and GGG_flag == False):
                        GGG_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGC" and GGC_flag == False):
                        GGC_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GGA_flag == True and GGT_flag == True and GGG_flag == True and GGC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = G

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GCA", "GCT", "GCG", "GCC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCA" and GCA_flag == False):
                        GCA_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCT" and GCT_flag == False):
                        GCT_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCG" and GCG_flag == False):
                        GCG_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCC" and GCC_flag == False):
                        GCC_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GCA_flag == True and GCT_flag == True and GCG_flag == True and GCC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = A

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GTA", "GTT", "GTG", "GTC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTA" and GTA_flag == False):
                        GTA_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTT" and GTT_flag == False):
                        GTT_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTG" and GTG_flag == False):
                        GTG_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTC" and GTC_flag == False):
                        GTC_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GTA_flag == True and GTT_flag == True and GTC_flag == True and GTG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = V

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CCA", "CCT", "CCG", "CCC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCA" and CCA_flag == False):
                        CCA_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCT" and CCT_flag == False):
                        CCT_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCG" and CCG_flag == False):
                        CCG_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCC" and CCC_flag == False):
                        CCC_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (CCA_flag == True and CCT_flag == True and CCG_flag == True and CCC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = P

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ACA", "ACT", "ACG", "ACC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACA" and ACA_flag == False):
                        ACA_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACT" and ACT_flag == False):
                        ACT_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACC" and ACC_flag == False):
                        ACC_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACG" and ACG_flag == False):
                        ACG_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (ACA_flag == True and ACT_flag == True and ACG_flag == True and ACC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = T

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ATA", "ATT", "ATC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATA" and ATA_flag == False):
                        ATA_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATT" and ATT_flag == False):
                        ATT_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATC" and ATC_flag == False):
                        ATC_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)


                    if (ATA_flag == True and ATT_flag == True and ATC_flag == True ):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = I

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TCA", "TCT", "TCG", "TCC", "AGC", "AGT")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCA" and TCA_flag == False):
                        TCA_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCT" and TCT_flag == False):
                        TCT_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCG" and TCG_flag == False):
                        TCG_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCC" and TCC_flag == False):
                        TCC_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGC" and AGC_flag == False):
                        AGC_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGT" and AGT_flag == False):
                        AGT_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (
                            TCA_flag == True and TCT_flag == True and TCC_flag == True and TCG_flag == True and AGC_flag == True and AGT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = S

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AGA", "AGG", "CGA", "CGT", "CGG", "CGC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGA" and AGA_flag == False):
                        AGA_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGG" and AGG_flag == False):
                        AGG_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGA" and CGA_flag == False):
                        CGA_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGT" and CGT_flag == False):
                        CGT_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGG" and CGG_flag == False):
                        CGG_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGC" and CGC_flag == False):
                        CGC_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (AGA_flag == True and AGG_flag == True and CGA_flag == True and CGT_flag == True and CGC_flag == True and CGG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = R

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TTA", "TTG", "CTA", "CTT", "CTG", "CTC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTA" and TTA_flag == False):
                        TTA_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTG" and TTG_flag == False):
                        TTG_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTA" and CTA_flag == False):
                        CTA_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTT" and CTT_flag == False):
                        CTT_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTG" and CTG_flag == False):
                        CTG_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTC" and CTC_flag == False):
                        CTC_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TTA_flag == True and TTG_flag == True and CTA_flag == True and CTG_flag == True and CTC_flag == True and CTT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = L
        for i in range(2, len(sequences), 3):
            if (i < len(sequences) - 2):
                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TAC", "TAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TAC" and TAC_flag == False):
                        TAC_flag = True
                        Y += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TAT" and TAT_flag == False):
                        TAT_flag = True
                        Y += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TAC_flag == True and TAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = Y

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TGT", "TGC")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGT" and TGT_flag == False):
                        TGT_flag = True
                        C += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGC" and TGC_flag == False):
                        TGC_flag = True
                        C += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TGC_flag == True and TGT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = C

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TGG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TGG" and TGG_flag == False):
                        TGG_flag = True
                        W += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)
                    if (TGG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = W

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ATG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATG" and ATG_flag == False):
                        ATG_flag = True
                        M += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (ATG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = M

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CAG", "CAA")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAG" and CAG_flag == False):
                        CAG_flag = True
                        Q += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAA" and CAA_flag == False):
                        CAA_flag = True
                        Q += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (CAG_flag == True and CAA_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = Q

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AAC", "AAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAC" and AAC_flag == False):
                        AAC_flag = True
                        N += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAT" and AAT_flag == False):
                        AAT_flag = True
                        N += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (AAC_flag == True and AAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = N

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GAT", "GAC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAT" and GAT_flag == False):
                        GAT_flag = True
                        D += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAC" and GAC_flag == False):
                        GAC_flag = True
                        D += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GAT_flag == True and GAC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = D

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GAA", "GAG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAA" and GAA_flag == False):
                        GAA_flag = True
                        E += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GAG" and GAG_flag == False):
                        GAG_flag = True
                        E += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GAA_flag == True and GAG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = E

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TTT", "TTC")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTT" and TTT_flag == False):
                        TTT_flag = True
                        F += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTC" and TTC_flag == False):
                        TTC_flag = True
                        F += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TTT_flag == True and TTC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = F

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AAA", "AAG")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAA" and AAA_flag == False):
                        AAA_flag = True
                        K += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AAG" and AAG_flag == False):
                        AAG_flag = True
                        K += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (AAA_flag == True and AAG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = K

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CAC", "CAT")):
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAC" and CAC_flag == False):
                        CAC_flag = True
                        H += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CAT" and CAT_flag == False):
                        CAT_flag = True
                        H += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (CAC_flag == True and CAT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = H

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GGA", "GGT", "GGG", "GGC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGA" and GGA_flag == False):
                        GGA_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGT" and GGT_flag == False):
                        GGT_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)
                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGG" and GGG_flag == False):
                        GGG_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GGC" and GGC_flag == False):
                        GGC_flag = True
                        G += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GGA_flag == True and GGT_flag == True and GGG_flag == True and GGC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = G

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GCA", "GCT", "GCG", "GCC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCA" and GCA_flag == False):
                        GCA_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCT" and GCT_flag == False):
                        GCT_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCG" and GCG_flag == False):
                        GCG_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GCC" and GCC_flag == False):
                        GCC_flag = True
                        A += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GCA_flag == True and GCT_flag == True and GCG_flag == True and GCC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = A

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("GTA", "GTT", "GTG", "GTC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTA" and GTA_flag == False):
                        GTA_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTT" and GTT_flag == False):
                        GTT_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTG" and GTG_flag == False):
                        GTG_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "GTC" and GTC_flag == False):
                        GTC_flag = True
                        V += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (GTA_flag == True and GTT_flag == True and GTC_flag == True and GTG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = V

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("CCA", "CCT", "CCG", "CCC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCA" and CCA_flag == False):
                        CCA_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCT" and CCT_flag == False):
                        CCT_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCG" and CCG_flag == False):
                        CCG_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CCC" and CCC_flag == False):
                        CCC_flag = True
                        P += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (CCA_flag == True and CCT_flag == True and CCG_flag == True and CCC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = P

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ACA", "ACT", "ACG", "ACC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACA" and ACA_flag == False):
                        ACA_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACT" and ACT_flag == False):
                        ACT_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACC" and ACC_flag == False):
                        ACC_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ACG" and ACG_flag == False):
                        ACG_flag = True
                        T += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (ACA_flag == True and ACT_flag == True and ACG_flag == True and ACC_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = T

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("ATA", "ATT", "ATC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATA" and ATA_flag == False):
                        ATA_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATT" and ATT_flag == False):
                        ATT_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "ATC" and ATC_flag == False):
                        ATC_flag = True
                        I += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)


                    if (ATA_flag == True and ATT_flag == True and ATC_flag == True ):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = I

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TCA", "TCT", "TCG", "TCC", "AGC", "AGT")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCA" and TCA_flag == False):
                        TCA_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCT" and TCT_flag == False):
                        TCT_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCG" and TCG_flag == False):
                        TCG_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TCC" and TCC_flag == False):
                        TCC_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGC" and AGC_flag == False):
                        AGC_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGT" and AGT_flag == False):
                        AGT_flag = True
                        S += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (
                            TCA_flag == True and TCT_flag == True and TCC_flag == True and TCG_flag == True and AGC_flag == True and AGT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = S

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("AGA", "AGG", "CGA", "CGT", "CGG", "CGC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGA" and AGA_flag == False):
                        AGA_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "AGG" and AGG_flag == False):
                        AGG_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGA" and CGA_flag == False):
                        CGA_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGT" and CGT_flag == False):
                        CGT_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGG" and CGG_flag == False):
                        CGG_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CGC" and CGC_flag == False):
                        CGC_flag = True
                        R += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (AGA_flag == True and AGG_flag == True and CGA_flag == True and CGT_flag == True and CGC_flag == True and CGG_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = R

                if (sequences[i] + sequences[i + 1] + sequences[i + 2] in ("TTA", "TTG", "CTA", "CTT", "CTG", "CTC")):

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTA" and TTA_flag == False):
                        TTA_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "TTG" and TTG_flag == False):
                        TTG_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTA" and CTA_flag == False):
                        CTA_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTT" and CTT_flag == False):
                        CTT_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTG" and CTG_flag == False):
                        CTG_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (sequences[i] + sequences[i + 1] + sequences[i + 2] == "CTC" and CTC_flag == False):
                        CTC_flag = True
                        L += self.amino_acid_freq(sequences[i], sequences[i + 1], sequences[i + 2], codon_dict)

                    if (TTA_flag == True and TTG_flag == True and CTA_flag == True and CTG_flag == True and CTC_flag == True and CTT_flag == True):
                        amino_acid_dist[codon_wheel[sequences[i]][sequences[i + 1]][sequences[i + 2]]] = L



        return amino_acid_dist




    def get_amino_acid_dist(self):
        """
        Return the expected amino acid distribution (fractions) based on the
        expected distribution (fractions) of the different codons.

        :return: Dictionary that contains the expected amino acid distribution.
                 The keys are the 20 different amino acids, the values are the
                 corresponding frequencies (rounded to 6 digits).
        """

        codon_dict = self.get_codon_dist()
        codon_wheel = self.get_codon_wheel()
        ANN = self.calculate_amino_freq(self.genome, codon_wheel, codon_dict)
        for aa in ANN.keys():
            ANN[aa] = round(ANN[aa], 6)
        return ANN



# gc = Genome("AGCCGTCCTGCCGGACTGCTGGAGGCGGCCACAGCGCCATGTTGGATGCTCTGCTCGTTGAGTGAAGAAAATCCACCGGCATCGCCTGAGCCCCGCTACCGAGAAGGGCGCCGCTTCCTCCGGGGAGGGGGATAAAGATCCCCCGCCGCCGGCCCATGAGGATATTGCCGTGAAAGGCACAGCGACTGCAGCAGGAACCGGACCCGGCACCGGAGCGGCGGCGGCGGCGGCAGCAGCGGTACCGCCTCCTCACCCGGCGGCGGCAGCAGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCAGCGGTCCCCCCTCCTCACCCGAACATCAGGGCCCTCCAGACTCAGGCGCCCCAACAAATTCCTAGAGGACCTGTGCAACAACCTCTTGAGGATCGAATCTTCACTCCCGCTGTCTCAGCAGTCTACAGCACGGTAACACAAGTGGCAAGACAGCCGGGAACCCCTACCCCATCCCCTTATTCAGCACATGAAATAAACAAGGGGCATCCAAATCTTGCGGCAACGCCCCCGGGACATGCATCGTCCCCTGGACTCTCTCAAACCCCTTATCCCTCTGGACAGAATGCAGGTCCAACCACGCTGGTATACCCTCAAACCCCTCAGACAATGAATTCACAACCTCAAACCCGTTCTCCGTTTTTCCAGAGGCCTCAAATACAGCCTCCTAGAGCTACCATCCCGAACAGCAGTCCTTCCATTCGTCCTGGTGCACAGACACCCACTGCAGTGTACCAGGCTAATCAGCACATCATGATGGTTAACCATCTGCCCATGCCGTACCCAGTGCCCCAGGGGCCTCAGTACTGTATACCACAGTACCGTCATAGTGGCCCTCCTTATGTTGGGCCCCCCCAACAATATCCAGTTCAACCACCGGGGCCAGGTCCTTTTTATCCTGGACCAGGACCTGGGGACTTCCCCAATGCTTATGGAACGCCTTTTTACCCAAGTCAGCCGGTGTATCAGTCAGCACCTATCATAGTGCCTACGCAGCAACAGCCGCCTCCAGCCAAGAGAGAGAAAAAAACTATAAGAATTCGGGATCCAAACCAGGGAGGTAAAGACATAACAGAGGAGATTATGTCTGGAGGTGGCAGCAGAAATCCTACTCCACCCATAGGAAGACCCACGTCCACACCTACTCCTCCTCAGCTGCCCAGCCAGGTCCCCGAGCACAGCCCTGTGGTTTATGGGACTGTGGAGAGCGCTCATCTTGCTGCCAGCACCCCTGTCACTGCAGCTAGCGACCAGAAGCAAGAGGAGAAGCCAAAACCAGATCCAGTGTTAAAGTCTCCTTCCCCAGTCCTTAGGCTAGTCCTCAGTGGAGAGAAGAAAGAACAAGAAGGCCAGACATCTGAAACTACTGCAATAGTATCCATAGCAGAGCTTCCTCTGCCTCCATCACCTACCACTGTTTCTTCTGTTGCTCGAAGTACAATTGCAGCCCCCACCTCTTCTGCTCTTAGTAGCCAACCAATATTCACCACTGCTATAGATGACAGATGTGAACTCTCATCCCCAAGAGAAGACACAATTCCTATACCCAGCCTCACATCTTGCACAGAAACATCAGACCCTTTACCAACAAATGAAAATGATGATGATATATGCAAGAAACCCTGTAGTGTAGCACCTAATGATATTCCACTGGTTTCTAGTACTAACCTAATTAATGAAATAAATGGAGTTAGCGAAAAATTATCAGCCACGGAGAGCATTGTGGAAATAGTAAAACAGGAAGTATTGCCATTGACTCTTGAATTGGAGATTCTCGAAAATCCCCCAGAAGAAATGAAACTGGAGTGTATCCCAGCTCCCATCACCCCTTCCACAGTTCCTTCCTTTCCTCCAACTCCTCCAACTCCTCCAGCTTCTCCTCCTCACACTCCAGTCATTGTTCCTGCTGCTGCCACTACTGTTAGTTCTCCGAGTGCTGCCATCACAGTCCAGAGAGTCCTAGAGGAGGACGAGAGCATAAGAACTTGCCTTAGTGAAGATGCAAAAGAGATTCAGAACAAAATAGAGGTAGAAGCAGATGGGCAAACAGAAGAGATTTTGGATTCTCAAAACTTAAATTCAAGAAGGAGCCCTGTCCCAGCTCAAATAGCTATAACTGTACCAAAGACATGGAAGAAACCAAAAGATCGGACCCGAACCACTGAAGAGATGTTAGAGGCAGAATTGGAGCTTAAAGCTGAAGAGGAGCTTTCCATTGACAAAGTACTTGAATCTGAACAAGATAAAATGAGCCAGGGGTTTCATCCTGAAAGAGACCCCTCTGACCTAAAAAAAGTGAAAGCTGTGGAAGAAAATGGAGAAGAAGCTGAGCCAGTACGTAATGGTGCTGAGAGTGTTTCTGAGGGTGAAGGAATAGATGCTAATTCAGGCTCCACAGATAGTTCTGGTGATGGGGTTACATTTCCATTTAAACCAGAATCCTGGAAGCCTACTGATACTGAAGGTAAGAAGCAGTATGACAGGGAGTTTCTGCTGGACTTCCAGTTCATGCCTGCCTGTATACAAAAACCAGAGGGCCTGCCTCCTATCAGTGATGTGGTTCTTGACAAGATCAACCAACCCAAATTGCCAATGCGAACTCTGGATCCTCGAATTTTGCCTCGAGGACCAGACTTTACACCAGCCTTTGCTGATTTTGGAAGGCAGACACCTGGTGGAAGAGGCGTACCTTTGTTGAATGTTGGGTCACGAAGATCTCAACCTGGCCAAAGAAGAGAACCCAGAAAGATCATCACAGTTTCTGTAAAAGAAGATGTACACCTGAAAAAGGCAGAAAATGCCTGGAAGCCAAGCCAAAAACGAGACAGCCAAGCCGATGATCCCGAAAACATTAAAACCCAGGAGCTTTTTAGAAAAGTTCGAAGTATCTTAAATAAATTGACACCACAGATGTTCAATCAACTGATGAAGCAAGTGTCAGGACTTACTGTTGACACAGAGGAGCGGCTGAAAGGAGTTATTGACCTGGTCTTTGAGAAGGCTATTGATGAACCCAGTTTCTCTGTGGCTTACGCAAACATGTGTCGATGTCTAGTAACGCTGAAAGTACCCATGGCAGACAAGCCTGGTAACACAGTGAATTTCCGGAAGCTGCTACTGAACCGTTGCCAGAAGGAGTTTGAAAAAGATAAAGCAGATGATGATGTCTTTGAGAAGAAGCAGAAAGAACTTGAGGCTGCCAGTGCTCCAGAGGAGAGGACAAGGCTTCATGATGAACTGGAAGAAGCCAAGGACAAAGCCCGGCGGAGATCCATTGGCAACATCAAGTTTATTGGAGAACTCTTTAAACTCAAAATGCTGACTGAAGCCATCATGCATGACTGTGTGGTGAAGCTGCTAAAGAACCATGATGAAGAATCCCTGGAGTGCCTGTGTCGCCTGCTCACCACCATTGGCAAAGACTTGGACTTTGAAAAAGCAAAGCCACGTATGGACCAGTACTTTAATCAGATGGAGAAAATTGTGAAAGAAAGAAAAACCTCATCTAGGATTCGGTTCATGCTTCAAGATGTTATAGACCTAAGGCTGTGCAATTGGGTATCTCGAAGAGCAGATCAAGGGCCTAAAACTATCGAACAGATTCACAAAGAGGCTAAAATAGAAGAACAAGAAGAGCAAAGGAAGGTCCAGCAACTCATGACCAAAGAGAAGAGAAGACCAGGTGTCCAGAGAGTGGACGAAGGTGGGTGGAACACTGTACAAGGGGCCAAGAACAGTCGGGTACTGGACCCCTCAAAATTCCTAAAAATCACTAAGCCTACAATTGATGAAAAAATTCAGCTGGTACCTAAAGCACAGCTAGGCAGCTGGGGAAAAGGCAGCAGTGGTGGAGCAAAGGCAAGTGAGACTGATGCCTTACGGTCAAGTGCTTCCAGTTTAAACAGATTCTCTGCCCTGCAACCTCCAGCACCCTCAGGGTCCACGCCATCCACGCCTGTAGAGTTTGATTCCCGAAGGACCTTAACTAGTCGTGGAAGTATGGGCAGGGAGAAGAATGACAAGCCCCTTCCATCTGCAACAGCTCGGCCAAATACTTTCATGAGGGGTGGCAGCAGTAAAGACCTGCTAGACAATCAGTCTCAAGAAGAGCAGCGGAGAGAGATGCTGGAGACCGTGAAGCAGCTCACAGGAGGTGTGGATGTGGAGAGGAACAGCACTGAGGCTGAGCGAAATAAAACAAGGGAGTCAGCAAAACCAGAAATTTCAGCAATGTCAGCTCATGACAAGGCTGCATTATCAGAAGAGGAACTGGAGAGGAAGTCGAAATCTATCATTGATGAATTTCTACACATTAATGATTTTAAGGAAGCCATGCAGTGTGTGGAAGAGCTGAATGCCCAGGGCCTACTACATGTTTTTGTGAGAGTGGGAGTGGAGTCCACCCTGGAAAGGAGCCAGATCACCAGGGATCACATGGGCCAATTACTCTATCAGCTGGTACAGTCAGAAAAACTCAGCAAACAGGACTTTTTCAAAGGTTTTTCAGAAACTTTGGAATTGGCAGATGACATGGCCATTGATATTCCCCATATTTGGTTGTACCTTGCTGAACTGGTGACCCCCATGTTAAAAGAAGGTGGAATCTCCATGAGAGAACTTACCATAGAATTTAGCAAACCTTTACTTCCTGTTGGAAGAGCTGGGGTCTTGCTATCTGAAATATTGCACCTACTATGCAAACAAATGAGCCATAAGAAAGTGGGAGCCTTATGGAGGGAGGCTGACCTCAGCTGGAAGGACTTTTTACCAGAAGGAGAAGATGTACATAATTTTCTTTTGGAGCAGAAGTTGGACTTCATAGAGTCTGACAGTCCCTGTTCCTCTGAAGCACTTTCAAAGAAAGAACTGTCTGCCGAAGAGCTGTATAAGCGACTCGAGAAACTCATTATTGAGGACAAAGCGAATGATGAACAGATCTTTGACTGGGTAGAGGCTAATCTAGACGAAATCCAGATGAGTTCACCTACATTCCTTAGAGCTTTAATGACTGCTGTTTGTAAAGCAGCTATTATAGATTTGCTGCGGATGTTTTTTGATTGTCTATATGACGAGGAGGTGATCTCCGAGGATGCCTTCTACAAATGGGAGAGCAGCAAGGACCCTGCAGAGCAGAATGGGAAGGGCGTGGCTCTGAAATCTGTCACGGCATTCTTCACGTGGCTGCGGGAAGCAGAAGAGGAGTCTGAGGATAACTAAAACTTCAAATACACAAAATGAAACAAAAGAAACAATTTAAGTATTTTTTTAAAAAGTTTCACGTCTTCGCCAATCACAGTGCAGCAAGGCCAATTCTCGCAGAAACCCCCACGTGTGCACGAGTGGGAGAGGGGAAAGAGAAAAAAAGGTGATCATGGAGGAAAAAGGTACTGGATAAAAGTAAACTTCAAACCTTAGGGCGGGAGCACTAAAACCAAAATACATGTATTATTTATAGAAAATATTTTCTGTTTTAATCTTTTCTTTTTAAACAAGGACTCATACTTAAAAAAATGTTTAGCAAAAAAAAAAAAAGTTGAGAACTTTTAATTTATTTTAAGGACTGCAAATGCCAGTGTAATTTTTTAATTTGCAGTTTCTGTAAACAACTTGTATAATAGAAAAGCAGAGAAATAAATTTCCCTCCCCTTCAAGATGCACCTCATGTTTGTTTTAAGGTATAGCATTTAGTCCAGATTTGAGAAAGTTTGGGGTGAACAAGGTAAGAAAGATTTTTTTTTTTTTGGCATCAAATCTTTCTGCCTGCCTCTCAGCTTGCTTCAGAAAATTTAAAAAATCACAATAGTAATCAAAACATACATAACATTGAAACAGAAGGAAATGCTGTGGACCACAGAACTCCAAGAATTGTTTAAAAAAAAAAAAGTGCTACCCTGAGAAAAGTACTCTTAATACTCTTGAAATCTTTAGAGCAACTTTAAGGCTTGTAAATACATAGAACAAATATTTAAAAAAACAAAAAGAAATTGACTCAGTACTATTTCTTTTCACTTTGAAAATATAAAGAACAAAATAAAGACAAACATTGCAAGTTTAAAA")
# # # gc.get_amino_acid_dist()
# codon_dist1 = Genome.get_codon_dist(gc)
# print(codon_dist1)
# amino_acid = Genome.get_amino_acid_dist(gc)
# print(amino_acid)


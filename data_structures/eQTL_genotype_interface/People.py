######################################################################
# People.py
# Data type encoding people, storing their genotypes at a list of
# loci of interest, and storing the alleles for which each person
# is heterozygous.
#
# Loci are assumed to be ordered in some fashion, represented as
# integers.
#
# People are also assumed to be ordered.
######################################################################

import numpy as np


def main():
    # test 1 - simple example
    maternal_genotypes = np.array([[0, 1, 0], [1, 1, 0]])
    paternal_genotypes = np.array([[0, 1, 1], [0, 1, 1]])

    People(maternal_genotypes, paternal_genotypes)

    # test 2 - provided sample data
    maternal_file = "Xm1.txt"
    paternal_file = "Xp1.txt"

    people = People(np.loadtxt(maternal_file), np.loadtxt(paternal_file))

    for i in range(people.maternal_genotypes.shape[0]):
        print(people.maternal_genotypes[i, :])
        print(people.paternal_genotypes[i, :])
        print(sorted(list(people.get_heterozygote_loci()[i])))
        print('*'*40)


class People():

    def __init__(self, maternal_genotypes, paternal_genotypes):
        """
        Initiate a set of people given the maternal and paternal 
        genotype (read from a genotype file).

        Genotype files are text files of space separated integers 
        representing genotypes. The loci are assumed to be of the 
        same order for maternal and paternal files. 
        """
        self.maternal_genotypes = maternal_genotypes
        self.paternal_genotypes = paternal_genotypes
        assert self.maternal_genotypes.shape == self.paternal_genotypes.shape, \
            "ERROR: maternal and paternal genotypes must be of the same dimensions!"
        self.heterozygotes = [set()
                              for _ in range(self.maternal_genotypes.shape[0])]
        for i, loc in zip(*np.where(self.maternal_genotypes != self.paternal_genotypes)):
            self.heterozygotes[i].add(loc)

    def get_heterozygote_loci(self):
        """
        Returns a list of sets where the i-th set contains the loci  
        for which person i is heterozygous. 
        """
        return self.heterozygotes


class Person():

    def __init__(self, maternal_genotype, paternal_genotype):
        """
        Store the maternal genotype, paternal genotype, and locations 
        where the person is heterozygous. 
        """
        self.maternal_genotype = maternal_genotype
        self.paternal_genotype = paternal_genotype
        self.locations_heterozygous = set()
        for index, (i, j) in enumerate(zip(maternal_genotype,
                                           paternal_genotype)):
            if i != j:
                self.locations_heterozygous.add(index)

        def get_maternal(self):
            return self.maternal_genotype

        def get_paternal(self):
            return self.paternal_genotype

        def get_het_loci(self):
            return self.locations_heterozygous


if __name__ == '__main__':
    main()

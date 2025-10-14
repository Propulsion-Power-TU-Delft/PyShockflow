import numpy as np
import json
import os

with open(os.path.dirname(os.path.realpath(__file__)) + "/group_valency.json") as file:
    group_valencies = json.load(file)

linear_molecules = ['linear']
ring_molecules = ['aromatic', 'cyclohexane', 'cyclopentane']
alkyl_groups = ['CH3', 'CH2', ')CH', ')CH(']
branch_groups = ['C_arom', 'CH_hex', 'CH_pent']
polar_groups = ['CH=0', ')C=O', 'OCH2', 'OCH3', 'COO']
associating_groups = ['OH', 'NH2']
linear_doublebond_groups = ['=CH2', '=CH', '=C(']


class MolecularStructure:

    def __init__(self, molecule_type, groups, group_occurrences):

        """
        :param molecule_type:
            'alkane'
            'alkene'
            'aromatic' (possibly with branch)
            'cyclohexane' (possibly with branch)
            'cyclopentane' (possibly with branch)

        """

        if isinstance(molecule_type, str):
            molecule_type = [molecule_type]

        for moltype in molecule_type:
            if moltype not in linear_molecules and moltype not in ring_molecules:
                raise Exception('Molecule Type not recognized')
        self.MoleculeType = molecule_type

        if np.shape(groups) != np.shape(group_occurrences):
            raise Exception('shape of groups and group occurrences is not compatible')

        self.Groups = groups
        self.GroupOccurrences = group_occurrences

        shape = np.shape(self.Groups)
        self.NCmp = shape[0]
        self.NGroups = shape[1]

    def get_open_ends(self):
        valencies = np.empty([self.NCmp])
        for i in range(self.NCmp):
            molecule_valency = 0
            for j in range(self.NGroups):
                molecule_valency += self.GroupOccurrences[i, j]*(2 - group_valencies[self.Groups[i, j]])
            valencies[i] = molecule_valency
        return valencies

    def check_octet_rule(self):
        open_ends = self.get_open_ends()
        octet = np.empty(self.NCmp)
        for i in range(self.NCmp):
            octet[i] = open_ends[i]
            if self.MoleculeType[i] in linear_molecules:
                octet[i] -= 2
        return octet

    def check_disjunct_molecule(self):
        check = np.empty([self.NCmp, self.NGroups])
        for i in range(self.NCmp):
            sumoccurrences = np.sum(self.GroupOccurrences[i, :][np.where(self.Groups[i, :] != 'None')])
            if self.MoleculeType[i] in linear_molecules:
                rhs = 2
                if any(True for ldbg in linear_doublebond_groups if (ldbg in self.Groups[i, :] and self.GroupOccurrences[i, :][np.where(self.Groups[i, :] == ldbg)] > 0)):
                    sumoccurrences -= 1 #considering the double bond structure such as 2 =CH2 as a single group.
            else:
                rhs = 0
            for j in range(self.NGroups):
                if self.Groups[i, j] == 'None':
                    check[i, j] = -1
                else:
                    check[i, j] = rhs - (sumoccurrences - self.GroupOccurrences[i, j]*(group_valencies[self.Groups[i, j]] - 1))  #negative if not violated
        return check

    def check_feasible_branching(self):
        feasible_branching = np.empty(self.NCmp)
        for i in range(self.NCmp):
            if self.MoleculeType[i] in linear_molecules:
                feasible_branching[i] = -1
            else:
                if any(True for alk in alkyl_groups if (alk in self.Groups[i, :] and self.GroupOccurrences[i, :][np.where(self.Groups[i, :] == alk)] > 0)):
                    if self.MoleculeType[i] == 'aromatic':
                        if 'C_arom' in self.Groups[i, :]:
                            feasible_branching[i] = 1 - self.GroupOccurrences[i, :][np.where(self.Groups[i, :] == 'C_arom')]
                        else:
                            feasible_branching[i] = 1
                    elif self.MoleculeType[i] == 'cyclohexane':
                        if 'CH_hex' in self.Groups[i, :]:
                            feasible_branching[i] = 1 - self.GroupOccurrences[i, :][np.where(self.Groups[i, :] == 'CH_hex')]
                        else:
                            feasible_branching[i] = 1
                    elif self.MoleculeType[i] == 'cyclopentane':
                        if 'CH_pent' in self.Groups[i, :]:
                            feasible_branching[i] = 1 - self.GroupOccurrences[i, :][np.where(self.Groups[i, :] == 'CH_pent')]
                        else:
                            feasible_branching[i] = 1
                else:
                    feasible_branching[i] = -1

        return feasible_branching

    def get_fitting_constraints(self):
        #this prevents extrapolation from the database used for fitting the GC parameters
        constr = np.empty([self.NCmp, 2])
        for i in range(self.NCmp):
            if self.MoleculeType[i] in linear_molecules:
                constr[i, 0] = -1
            else:
                constr[i, 0] = 0.0
                count = 0
                for group in self.Groups[i, :]:
                    if any([group in polar_groups, group in associating_groups, group in linear_doublebond_groups]):
                        constr[i, 0] += self.GroupOccurrences[i, count]
                    count += 1
            ngroups = 0
            for group in self.Groups[i, :]:
                if any([group in polar_groups, group in associating_groups, group in linear_doublebond_groups]):
                    ngroups += 1
            if ngroups > 1:
                constr[i, 1] = ngroups - 1
            else:
                constr[i, 1] = -1

        return constr

    def complete_open_ends(self, completion_groups):
        octet_rule = self.check_octet_rule()
        groups_completion = np.full([self.NCmp, 1], 'None',  dtype=self.Groups.dtype)
        occurrences_completion = np.full([self.NCmp, 1], -1.0)
        append = False
        for i in range(self.NCmp):
            if octet_rule[i] < 0:
                if completion_groups[i] not in self.Groups[i, :]:
                    append = True
                    groups_completion[i, 0] = completion_groups[i]
                    occurrences_completion[i, 0] = -octet_rule[i]/(2 - group_valencies[completion_groups[i]])
                else:

                    index = np.where(self.Groups[i, :] == completion_groups[i])
                    self.GroupOccurrences[i, index] -= octet_rule[i]/(2 - group_valencies[completion_groups[i]])

        if append:
            self.Groups = np.append(self.Groups, groups_completion, axis = 1)
            self.GroupOccurrences = np.append(self.GroupOccurrences, occurrences_completion, axis = 1)
            self.NGroups += 1

    def check_max_groups(self, max_groups):
        constr = np.empty(self.NCmp)
        for i in range(self.NCmp):
            tmp = self.GroupOccurrences[i, :][np.where(self.Groups[i] != 'None')]
            constr[i] = np.sum(tmp) - max_groups

        return constr

    def check_min_groups(self, min_groups):
        constr = np.empty(self.NCmp)
        for i in range(self.NCmp):
            tmp = self.GroupOccurrences[i, :][np.where(self.Groups[i] != 'None')]
            constr[i] = min_groups - np.sum(tmp)

        return constr













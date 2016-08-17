"""
A package for keeping track of energy flow and so forth
"""
import numpy

try:
    from . import Theory
    from . import Utilities
except:
    import Theory
    import Utilities


class Feynman:
    """
    Please remember that this is essentially the subset sum problem which
    is NP-Complete.
    """
    def __init__(self, theory: Theory.Theory):
        interactions = theory.getInt()
        fieldNames = set()
        for inter in interactions:
            for f in inter.fields:
                fieldNames.add(f.name)
        # Gives the ordered list of fields
        self.indexf = sorted(fieldNames)

        # Gives the indices of fields
        self.findex = {self.indexf[i]: 2*i for i in range(len(self.indexf))}
        self.inters = []
        for inter in interactions:
            outvecinter = [0 for _ in range(2*len(self.indexf))]
            for f in inter.fields:
                outvecinter[self.findex[f.name] + (1 if f.anti else 0)] += 1
            addedset = set()
            for i, o in Utilities.genParticlePerms(outvecinter):
                minusi = -i
                normo = o
                checka = (inter.coupling, tuple(minusi), tuple(normo))
                if checka not in addedset:
                    addedset.add(checka)
                    self.inters.append((inter.coupling,
                                        minusi,
                                        normo))
                # antiminusi = Utilities.swapAntis(-i)
                # antinormo = Utilities.swapAntis(o)
                # checkb = (inter.coupling, tuple(antiminusi), tuple(antinormo))
                # if checkb not in addedset:
                #     addedset.add(checkb)
                #     self.inters.append((inter.coupling,
                #                         antiminusi,
                #                         antinormo))

    def convertDictToVec(self, d):
        vec = [0 for _ in range(2 * len(self.indexf))]
        for name,anti in d:
            if not anti:
                vec[self.findex[name]] = d[(name,anti)]
            else:
                vec[self.findex[name]+1] = d[(name,anti)]
        return numpy.array(vec)

    def listPaths(self, startstate: dict, endstate: dict, maxorder=4):
        """
        Lists each path to get to given endstate (endstate * state > 0)
        :param startstate: Start state  as dict{(fieldname, anti?): num_particles}
        :param endstate: State of interest as dict{(fieldname, anti?): num_particles}
        :param maxorder: maximum number of nodes allowed
        :return: [paths]
        """
        startvec = self.convertDictToVec(startstate)
        endvec = self.convertDictToVec(endstate)
        paths = Utilities.vectorPath(startvec,endvec,[(v[1],v[2]) for v in self.inters],maxorder)
        return [[self.inters[j] for j in path] for path in paths]

    def listTrees(self, startstate: dict, endstate: dict, maxorder=4):
        """
        Lists each path to a given endstate without loops
        :param startstate: Start state  as dict{field: num_particles}
        :param endstate: State of interest as dict{field: num_particles}
        :param maxorder: maximum number of nodes allowed
        :return: [paths]
        """
        pass

    def amplitude(self, startstate: dict, endstate: dict):
        """
        Calculates the total amplitude of the endstate
        :param startstate: Start state  as dict{field: num_particles}
        :param endstate: State of interest as dict{field: num_particles}
        :return: A
        """
        pass

    def generateLatexForPaths(self,paths):
        pass
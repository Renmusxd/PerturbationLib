"""
A package for keeping track of energy flow and so forth
"""
try:
    from . import Theory
except:
    import Theory


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
                fieldNames.update(f.name)
        self.indexf = sorted(fieldNames) # Memory is cheap
        self.findex = {self.indexf[i]: i for i in range(len(self.indexf))}
        self.inters = {}
        for inter in interactions:
            vecinter = [0 for _ in range(len(self.indexf))]
            for f in inter.fields:
                vecinter[self.findex[f.name]] += -1 if f.anti else 1
            self.inters[inter.coupling] = vecinter

    def listPaths(self, startstate: dict, endstate: dict):
        """
        Lists each path to get to given endstate (endstate * state > 0)
        :param startstate: Start state  as dict{field: num_particles}
        :param endstate: State of interest as dict{field: num_particles}
        :return: [paths]
        """
        startvec = [0 for _ in range(len(self.indexf))]
        for k in startstate:
            startvec[self.findex[k]] = startstate[k]

    def listTrees(self, startstate: dict, endstate: dict):
        """
        Lists each path to a given endstate without loops
        :param startstate: Start state  as dict{field: num_particles}
        :param endstate: State of interest as dict{field: num_particles}
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

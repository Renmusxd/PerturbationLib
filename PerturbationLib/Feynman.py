"""
A package for keeping track of energy flow and so forth
"""
from PerturbationLib.Theory import *
from PerturbationLib.Utilities import *
from typing import Generator, Mapping


class Feynman:
    """
    Please remember that this is essentially the subset sum problem which
    is NP-Complete.
    """
    def __init__(self, theory: Theory):
        self.interactions = theory.getInt()

        # All fields and antifields
        fieldnames = set(f for inter in self.interactions
                         for field in inter.getFields()
                         for f in [field, field.antifield()])
        # Gives the ordered list of fields and antifields, weird key is to place all matter fields before antimatter.
        self.fieldlist = list(sorted(fieldnames, key=lambda x: str(int(x.anti)) + ":" + repr(x)))

        # Gives the indices of fields
        self.findex = {field: i for i, field in enumerate(self.fieldlist)}
        swapvec = [self.findex[field.antifield()] for field in self.fieldlist]
        self.inters = []
        for inter in self.interactions:
            outvecinter = [0 for _ in range(len(self.fieldlist))]
            for f in inter.fields:
                outvecinter[self.findex[f]] += 1

            addedset = set()
            for i, o in Utilities.genInOutPairs(outvecinter, swapvec=swapvec):
                taken, given = -i, o
                checka = (inter, tuple(taken), tuple(given))
                if checka not in addedset:
                    addedset.add(checka)
                    self.inters.append((inter, taken, given))

    def convertDictToVec(self, d: Mapping[Field, int]) -> Sequence[int]:
        vec = [0 for _ in range(len(self.fieldlist))]
        for field in d:
            vec[self.findex[field]] = d[field]
        return numpy.array(vec)

    def listPaths(self, startstate: Mapping[Field, int], endstate: Mapping[Field, int],
                  maxorder=4) -> Generator['Diagram', None, None]:
        """
        Lists each path to get to given endstate (endstate * state > 0)
        :param startstate: Start state  as dict{field: num_particles}
        :param endstate: State of interest as dict{field: num_particles}
        :param maxorder: maximum number of nodes allowed
        :return: [paths]
        """
        startvec = self.convertDictToVec(startstate)
        endvec = self.convertDictToVec(endstate)
        for path in Utilities.genVectorPaths(startvec, endvec, {v[0]: (v[1], v[2]) for v in self.inters}, maxorder):
            diagrams = [Diagram(*[f for f in startstate for _ in range(startstate[f])])]
            for inter, req, prod in path:
                removed_fields = {self.fieldlist[i]: req[i] for i in range(len(self.fieldlist)) if req[i] > 0}
                added_fields = {self.fieldlist[i]: prod[i] for i in range(len(self.fieldlist)) if prod[i] > 0}


class Vertex:
    def __init__(self, inter: Interaction):
        self.fields = tuple(sorted(inter.getRawFields()))
        self.free_fields = {field: self.fields.count(field) for field in self.fields}
        self.field_paths = {}
        self.vertex_paths = {}

    def addPath(self, field: Field, vertex: 'Vertex'):
        field = field.matterField()
        if field not in self.free_fields:
            raise ValueError("Field not in vertex.")
        if self.free_fields[field] == 0:
            raise ValueError("No field paths remaining.")
        self.free_fields[field] -= 1

        if field not in self.field_paths:
            self.field_paths[field] = []
        if vertex not in self.vertex_paths:
            self.vertex_paths[vertex] = []
        self.field_paths[field].append(vertex)
        self.vertex_paths[vertex].append(field)

    def __hash__(self):
        return hash(self.fields)

    def __eq__(self, other):
        if not isinstance(other, Vertex):
            raise ValueError("Cannot compare non-like types.")
        return self.fields == other.fields

    def __lt__(self, other):
        if not isinstance(other, Vertex):
            raise ValueError("Cannot compare non-like types.")
        return self.fields < other.fields

    def __le__(self, other):
        if not isinstance(other, Vertex):
            raise ValueError("Cannot compare non-like types.")
        return self.fields <= other.fields

    def __gt__(self, other):
        if not isinstance(other, Vertex):
            raise ValueError("Cannot compare non-like types.")
        return self.fields > other.fields

    def __ge__(self, other):
        if not isinstance(other, Vertex):
            raise ValueError("Cannot compare non-like types.")
        return self.fields >= other.fields


class Diagram:
    def __init__(self, *fields: Field):
        self.fields = tuple(sorted(fields))

    def applyInteraction(self, interaction: Interaction,
                         remove_particles: Mapping[Field, int],
                         add_particles: Mapping[Field, int]) -> Generator['Diagram', None, None]:
        pass

    def serialize(self) -> str:
        pass

    def __eq__(self, other):
        if not isinstance(other, Diagram):
            raise ValueError("Cannot compare Diagram and {}".format(type(other)))
        return self.serialize() == other.serialize()

    def __lt__(self, other):
        if not isinstance(other, Diagram):
            raise ValueError("Cannot compare Diagram and {}".format(type(other)))
        return self.serialize() < other.serialize()

    def __le__(self, other):
        if not isinstance(other, Diagram):
            raise ValueError("Cannot compare Diagram and {}".format(type(other)))
        return self.serialize() <= other.serialize()

    def __gt__(self, other):
        if not isinstance(other, Diagram):
            raise ValueError("Cannot compare Diagram and {}".format(type(other)))
        return self.serialize() > other.serialize()

    def __ge__(self, other):
        if not isinstance(other, Diagram):
            raise ValueError("Cannot compare Diagram and {}".format(type(other)))
        return self.serialize() >= other.serialize()

    def __hash__(self):
        return hash(self.serialize())

    def __repr__(self):
        pass


from PerturbationLib import Utilities
from PerturbationLib.Symmetries import Symmetry
from typing import Iterable, MutableMapping, List, Sequence


class Field:
    """
    A class containing information about a given field,
    maintains a list of multiplets for each symmetry as well
    as other important information.
    """
    def __init__(self, name: str, *symmetries: Symmetry, anti: bool = False):
        self.syms = list(symmetries)
        self.name = name
        self.anti = anti

    def combineWithSyms(self, syms: Iterable[Symmetry]) -> List[Symmetry]:
        '''
        Returns results of [syms] * [field.syms] as list of syms (each may be sum of reprs)
        :param syms:
        :return: [syms]
        '''
        combinedsyms = []

        for os in syms:
            sym = None
            for s in self.syms:
                if os.name == s.name:
                    sym = s
                    break
            if sym is None:
                raise Exception("Fields do not share symmetry: "+os.name)
            combinedsyms.append(os.combine(sym))
        return combinedsyms

    def antifield(self) -> 'Field':
        '''
        :return: antifield with correct symmetry multiplets
        '''
        antisyms = [s.inverse() for s in self.syms]
        return Field(self.name, *antisyms, anti=(not self.anti))

    def copy(self) -> 'Field':
        return Field(self.name, *self.syms, anti=self.anti)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return type(self) == type(other) and repr(self) == repr(other)

    def __repr__(self):
        return "\\bar{"+self.name+"}" if self.anti else self.name

    def _latex(self, *args):
        return self.__repr__()

    def _repr_latex_(self):
        return "${}$".format(self._latex())


class Interaction:
    """
    Represents an interaction term, notably a set of fields and
    an associated coupling constant
    """
    def __init__(self, fields: Iterable[Field], coupling: str):
        self.fields = fields
        self.coupling = coupling
        self.intRepr = []
        for f in fields:
            newlist = []
            added = False
            for fname, fdel, fcount in self.intRepr:
                if f.name == fname:
                    newlist.append((fname, fdel + (-1 if f.anti else 1), fcount+1))
                    added = True
                else:
                    newlist.append((fname, fdel, fcount))
            if not added:
                newlist.append((f.name, (-1 if f.anti else 1), 1))
            self.intRepr = newlist
        self.intRepr = tuple([(a, abs(b), c) for (a, b, c) in self.intRepr])

    def __eq__(self, other):
        return (isinstance(other,self.__class__)) \
               and other.intRepr == self.intRepr

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.intRepr)

    def __repr__(self):
        return self.coupling + "".join([repr(f) for f in self.fields])

    def _latex(self, *args):
        return self.__repr__()

    def _repr_latex_(self):
        return "${}$".format(self._latex())


class Theory:
    """
    Maintains a list of applied symmetries and added fields,
    calculated allowed terms up to a truncation level.
    """
    def __init__(self, *symmetries: Symmetry, fields: Sequence[Field] = None,
                 gaugefields: MutableMapping[Symmetry, Field] = None,
                 trunc: int = 4):
        self.syms = symmetries
        self.fields = []
        self.trunc = trunc
        self.Lk = None
        self.Lint = None
        if gaugefields is None:
            self.gaugeFields = {}
        else:
            self.gaugeFields = gaugefields

        if fields is not None:
            for field in fields:
                self.addField(field)

    def addField(self, field: Field, gaugeforsym: Symmetry = None):
        if not gaugeforsym and self.fieldAllowed(field):
            self.fields.append(field)
        elif gaugeforsym:
            # gaugeforsym is the name of the symmetry to
            # which the field lends locality
            self.gaugeFields[gaugeforsym] = field
        else:
            raise Exception("Field does not have required symmetries")
        self.Lk = None
        self.Lint = None

    def getL(self) -> List[Interaction]:
        """
        :return: whole lagrangian
        """
        return self.getK() + self.getInt()

    def getK(self) -> List[Interaction]:
        """
        :return: return kinetic terms
        """
        if (self.Lk is None) or (self.Lint is None):
            self.calculateL()
        return self.Lk

    def getInt(self) -> List[Interaction]:
        """
        :return: interactions terms
        """
        if (self.Lk is None) or (self.Lint is None):
            self.calculateL()
        return self.Lint

    def calculateL(self) -> None:
        """
        Populate with all allowed interactions and kinetic terms
        """
        self.Lk = []
        self.Lint = []
        # TODO calculate kinetic terms

        # Make all combinations of trunc #fields and trunc #antifields
        intnum = 0
        Lintset = set()
        lkset = set()
        fieldweights = [(f, self.trunc) for f in self.fields] + [(f.antifield(), self.trunc) for f in self.fields]
        for encoding in Utilities.truncCombinations(fieldweights, self.trunc):
            accsyms = [s.singlet() for s in self.syms]  # Start with all singlets
            for field in encoding:
                accsyms = field.combineWithSyms(accsyms)
            # Make sure there's a singlet in each symmetry
            allhavesinglets = True
            for s in accsyms:
                allhavesinglets = allhavesinglets and s.containsSinglet()
            if allhavesinglets:
                if len(encoding) != 2 or encoding[0].name != encoding[1].name:
                    couple = Interaction(encoding, "g_{"+str(intnum)+"}")
                    if couple not in Lintset:
                        Lintset.add(couple)
                        intnum += 1
                else:
                    lkset.add(Interaction(encoding, "m_{"+str(encoding[0].name)+"}"))
        self.Lk = sorted(list(lkset), key=lambda k: len(k.fields))
        self.Lint = sorted(list(Lintset), key=lambda k: len(k.fields))

    def fieldAllowed(self, field: Field) -> bool:
        '''
        Returns whether a field may be added to the theory, precisely whether
        it contains the appropriate symmetries
        '''
        for f in self.syms:
            n = f.name
            found = False
            for fp in field.syms:
                if fp.name == n:
                    found = True
                    break
            if not found:
                return False
        return True

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if (self.Lk is None) or (self.Lint is None):
            self.calculateL()
        s = []
        if self.Lk is None:
            s.append("\mathcal{L}_k")
        else:
            s.append(" + ".join([repr(k) for k in self.Lk]))
        if self.Lint is None:
            s.append("\mathcal{L}_{int}")
        else:
            s.append(" + ".join([repr(i) for i in self.Lint]))
        return " + ".join(s)

    def _latex(self, *args):
        return self.__repr__()

    def _repr_latex_(self):
        return "${}$".format(self._latex())


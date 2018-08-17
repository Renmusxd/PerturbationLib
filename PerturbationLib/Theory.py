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
        """
        Returns results of [syms] * [field.syms] as list of syms (each may be sum of reprs)
        :param syms:
        :return: [syms]
        """
        combinedsyms = []

        for osym in syms:
            sym = None
            for s in self.syms:
                if osym.name == s.name:
                    sym = s
                    break
            if sym is None:
                raise Exception("Fields do not share symmetry: "+osym.name)
            combinedsyms.append(osym.combine(sym))
        return combinedsyms

    def antifield(self) -> 'Field':
        """
        :return: antifield with correct symmetry multiplets
        """
        antisyms = [s.inverse() for s in self.syms]
        return Field(self.name, *antisyms, anti=(not self.anti))

    def matterField(self) -> 'Field':
        return self.antifield() if self.anti else self.copy()

    def antiMatterField(self) -> 'Field':
        return self.antifield() if not self.anti else self.copy()

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
        self.fields = list(fields)
        self.coupling = coupling
        self.intRepr = []
        for f in self.fields:
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

    def getFields(self) -> List[Field]:
        return self.fields

    def getRawFields(self) -> List[Field]:
        return [field.matterField() for field in self.fields]

    def __eq__(self, other):
        return (isinstance(other, self.__class__)) \
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

    def getK(self, filter_anti_dups=True) -> List[Interaction]:
        """
        :return: return kinetic terms
        """
        if (self.Lk is None) or (self.Lint is None):
            self.calculateL(filter_anti_dups=filter_anti_dups)
        return self.Lk

    def getInt(self, filter_anti_dups=True) -> List[Interaction]:
        """
        :return: interactions terms
        """
        if (self.Lk is None) or (self.Lint is None):
            self.calculateL(filter_anti_dups=filter_anti_dups)
        return self.Lint

    def calculateL(self, filter_anti_dups=True) -> None:
        """
        Populate with all allowed interactions and kinetic terms
        """
        self.Lk = []
        self.Lint = []

        # Make all combinations of trunc #fields and trunc #antifields
        intnum = 0
        fieldweights = [(f, self.trunc) for f in self.fields] + [(f.antifield(), self.trunc) for f in self.fields]
        seen_l = set()
        for encoding in sorted(Utilities.truncCombinations(fieldweights, self.trunc), key=lambda k: len(k)):
            field_name_tuple = tuple(sorted([field.name for field in encoding]))
            if filter_anti_dups and field_name_tuple in seen_l:
                continue

            accsyms = [s.singlet() for s in self.syms]  # Start with all singlets
            for field in encoding:
                accsyms = field.combineWithSyms(accsyms)

            # Make sure there's a singlet in each symmetry
            allhavesinglets = True
            for s in accsyms:
                allhavesinglets = allhavesinglets and s.containsSinglet()

            if allhavesinglets:
                seen_l.add(field_name_tuple)
                if len(encoding) != 2 or encoding[0].name != encoding[1].name:
                    self.Lint.append(Interaction(encoding, "g_{"+str(intnum)+"}"))
                    intnum += 1
                else:
                    self.Lk.append(Interaction(encoding, "m_{"+str(encoding[0].name)+"}"))

    def fieldAllowed(self, field: Field) -> bool:
        """
        Returns whether a field may be added to the theory, precisely whether
        it contains the appropriate symmetries
        """
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


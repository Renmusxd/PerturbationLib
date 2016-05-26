from Utilities import truncCombinations

class Theory:
    """
    Maintains a list of applied symmetries and added fields,
    calculated allowed terms up to a truncation level.
    """

    def __init__(self,symmetries,fields=None,gaugeFields=None,trunc=4):
        if fields is None:
            self.fields = []
        else:
            self.fields = fields
        if gaugeFields is None:
            self.gaugeFields = []
        else:
            self.gaugeFields = gaugeFields
        self.syms = symmetries
        self.trunc = trunc
        self.Lk = None
        self.Lint = None


    def addField(self, field, gaugeField=False):
        if not gaugeField and self.fieldAllowed(field):
            self.fields.append(field)
        elif gaugeField:
            # TODO check
            self.gaugeFields.append(gaugeField)
        else:
            raise Exception("Field does not have required symmetries")
        self.Lk = None
        self.Lint = None

    def getL(self):
        """
        :return: whole lagrangian
        """
        return self.getK() + self.getInt()

    def getK(self):
        """
        :return: return kinetic terms
        """
        if (self.Lk is None) or (self.Lint is None):
            self.calculateL()
        return self.Lk

    def getInt(self):
        """
        :return: interactions terms
        """
        if (self.Lk is None) or (self.Lint is None):
            self.calculateL()
        return self.Lint

    def calculateL(self):
        """
        Populate with all allowed interactions and kinetic terms
        """
        self.Lk = []
        self.Lint = []
        # TODO calculate kinetic terms

        # Make all combinations of trunc #fields and trunc #antifields
        intNum = 0
        LintSet = set()
        fieldWeights = [(f,self.trunc) for f in self.fields] + [(f.antifield(),self.trunc) for f in self.fields]
        for encoding in truncCombinations(fieldWeights,self.trunc):
            accSyms = [s.singlet() for s in self.syms] # Start with all singlets
            for field in encoding:
                accSyms = field.precombineWithSyms(accSyms)
            # Make sure there's a singlet in each symmetry
            allHaveSinglets = True
            for s in accSyms:
                allHaveSinglets = allHaveSinglets and s.containsSinglet()
            if allHaveSinglets:
                if len(encoding)!=2 or encoding[0].name != encoding[1].name:
                    couple = Interaction(encoding,"g_{"+str(intNum)+"}")
                    if not couple in LintSet:
                        LintSet.add(couple)
                        intNum += 1
                else:
                    LintSet.add(Interaction(encoding,"m_{"+str(encoding[0].name)+"}"))
        self.Lint = sorted(list(LintSet), key=lambda k: len(k.fields))

    def fieldAllowed(self,field):
        '''
        Returns whether a field may be added to the theory, precicely whether
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
        s = ""
        if self.Lk is None:
            s += "\mathcal{L}_k"
        else:
            s += " + ".join([repr(k) for k in self.Lk])
        if self.Lint is None:
            s += "\mathcal{L}_{int}"
        else:
            s += " + ".join([repr(i) for i in self.Lint])
        return s

class Field:
    """
    A class containing information about a given field,
    maintains a list of multiplets for each symmetry as well
    as other important information.
    """
    def __init__(self,name,symmetries=None,anti=False):
        if symmetries is None:
            self.syms = []
        else:
            self.syms = symmetries
        self.name = name
        self.anti = anti

    def precombineWithSyms(self,syms):
        '''
        Returns results of [syms] * [field.syms] as list of syms (each may be sum of reprs)
        :param syms:
        :return: [syms]
        '''
        combinedSyms = []
        for s in self.syms:
            osym = None
            for os in syms:
                if os.name == s.name:
                    osym = os
                    break
            if osym is None:
                raise Exception("Fields do not share symmetry: "+s.name)
            combinedSyms.append(s.precombine(osym))
        return combinedSyms

    def antifield(self):
        '''
        :return: antifield with correct symmetry multiplets
        '''
        antisyms = [s.inverse() for s in self.syms]
        return Field(self.name,antisyms,anti = not self.anti)

    def copy(self):
        return Field(self.name,self.syms.copy(),self.anti)

    def __repr__(self):
        return "\\bar{"+self.name+"}" if self.anti else self.name


class Interaction:
    """
    Represents an interaction term, notably a set of fields and
    an associated coupling constant
    """
    def __init__(self, fields, coupling):
        self.fields = fields
        self.coupling = coupling
        self.intRepr = []
        for f in fields:
            newlist = []
            added = False
            for (fname,fdel,fcount) in self.intRepr:
                if (f.name == fname):
                    newlist.append( (fname,fdel + (-1 if f.anti else 1), fcount+1))
                    added = True
                else:
                    newlist.append((fname,fdel,fcount))
            if not added:
                newlist.append((f.name,(-1 if f.anti else 1),1))
            self.intRepr = newlist
        self.intRepr = tuple([(a,abs(b),c) for (a,b,c) in self.intRepr])

    def __eq__(self, other):
        return (isinstance(other,self.__class__)) \
               and other.intRepr == self.intRepr

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.intRepr)

    def __repr__(self):
        return self.coupling + "".join([repr(f) for f in self.fields])


from PerturbationLib import Symmetries


class SU(Symmetries.Symmetry):
    """
    SU(N)
    """
    def __init__(self, N, name=None, multiplets=None):
        if multiplets is None:
            multiplets = []
        super().__init__(N,name or 'SU('+str(N)+')', multiplets)
        for m in multiplets:
            if (type(m)==tuple or type(m)==list):
                if len(m)!=(self.N-1):
                    raise Exception("Multiplet of incorrect size")
            elif type(m)==int:
                # TODO
                raise Exception("Cannot yet convert integers to multiplets")
            else:
                raise Exception("Incorrect multiplet type: "+str(type(m))+":"+str(m))

    def containsSinglet(self):
        """
        :return: true if repr contains a singlet
        """
        return tuple((0 for _ in range(self.N-1))) in self.multiplets

    def singlet(self):
        """
        Gives a singlet of the given group
        :return:
        """
        return SU(self.N, self.name, [tuple((0 for _ in range(self.N-1)))])

    def inverse(self):
        """
        Gives the conjugate representation
        :return:
        """
        return SU(self.N,self.name,[m[::-1] for m in self.multiplets])

    def combineRepr(self, r1, r2):
        """
        Combines two representations of the symmetry.
        :return: [combined repr]
        """
        t1 = Tableau(rep=r1)
        t2 = Tableau(rep=r2)
        ts = t1.combine(t2)
        return [t.getMultiplet() for t in ts]


class Tableau:
    """
    A tableau object for SU(N) combinations
    """
    def __init__(self, rep=None, rowshape=None, lettercols=None, letterrows=None):
        self.rep = rep
        self.rows = tuple(rowshape) if rowshape else None
        if self.rep is None and self.rows is None:
            raise Exception("Multiplet or Tableau shape must be defined")
        self.N = len(rep)+1 if rep else len(rowshape)
        if lettercols is None:
            self.lettercols = []
        else:
            self.lettercols = lettercols
        if letterrows is None:
            self.letterrows = [{} for _ in range(self.N)]
        else:
            self.letterrows = letterrows

    def combine(self, other) -> list:
        """
        Combine this tableau with a single other
        :param other: another Tableau instance
        :return: [Tableau]
        """
        structure = self.getRows()
        letters = list(other.getRows())
        takeList = [self]
        putList = []
        for i in range(len(letters)):
            while letters[i]>0:
                # Take a letter 'i'
                letters[i] -= 1
                for r in range(i,len(structure)):
                    for table in takeList:
                        rowlen = table.rowLen(r)
                        if table.canPlaceLetter(i,r,rowlen):
                            clone = table.clone()
                            clone.placeLetter(i,r,rowlen)
                            putList.append(clone)
                takeList = putList
                putList = []
        for obj in takeList:
            obj.solidify()
        return takeList

    def rowLen(self, r: int) -> int:
        """
        Returns length of row plus letters
        :param r: row index
        :return: length of row
        """
        return self.getRows()[r] + sum([self.letterrows[r][k] for k in self.letterrows[r]])

    def canPlaceLetter(self, i: int, r: int, c: int) -> bool:
        # Must be adjacent to existing structure
        if self.rowLen(r) != c:
            # Must be placed directly to the right of structure
            return False
        if r > 0 and self.rowLen(r-1) <= c:
            # Cannot exceed length of row above
            return False

        # Must not have letter already in column
        if c < len(self.lettercols) and i in self.lettercols[c]:
            return False

        # Check if enough previous letters up to this line
        if i > 0:
            prevtot = 0
            for j in range(r):
                prv = 0
                cur = 0
                if (i-1) in self.letterrows[j]:
                    prv = self.letterrows[j][i-1]
                if i in self.letterrows[j]:
                    cur = self.letterrows[j][i]
                prevtot += prv - cur
                if prevtot < 0:
                    raise Exception("More {0} than {1} by row {2}".format(i,i-1,j))
            if prevtot == 0:
                return False

        return True

    def placeLetter(self, i, r, c):
        if self.canPlaceLetter(i,r,c):
            if i in self.letterrows[r]:
                self.letterrows[r][i] += 1
            else:
                self.letterrows[r][i] = 1
            if c >= len(self.lettercols):
                self.lettercols += [[]] * (c+1 - len(self.lettercols))
            self.lettercols[c].append(i)
            # We know it's not yet in this column
        else:
            raise Exception("Cannot place letter")

    def solidify(self):
        """
        Convert letters into rep and cut down to size
        """
        self.rows = tuple([self.rowLen(r) for r in range(self.N)])
        self.rep = None
        self.lettercols = []
        self.letterrows = [{} for _ in range(self.N)]

    def getMultiplet(self) -> tuple:
        if (self.rep is None) and (self.rows is not None):
            rep = []
            for i in range(self.N-1):
                rep.append(self.rows[i] - self.rows[i+1])
            self.rep = tuple(rep)
        return self.rep

    def getRows(self) -> tuple:
        if (self.rows is None) and (self.rep is not None):
            rar = [0 for _ in range(self.N)]
            for i in range(self.N - 2, -1, -1):
                for j in range(i + 1):
                    rar[j] += self.rep[i]
            m = min(rar)
            for i in range(self.N):
                rar[i] -= m
            m = max(rar)
            if m == 0:
                rar = [1 for _ in range(self.N)]
            self.rows = tuple(rar)
        return self.rows

    def clone(self):
        return Tableau(rowshape=self.getRows(),
                       lettercols=[l.copy() for l in self.lettercols],
                       letterrows=[d.copy() for d in self.letterrows])

    def __eq__(self, other) -> bool:
        if type(other)==type(self):
            return str(other)==str(self)
        return False

    def __hash__(self) -> int:
        return hash(str(self))

    def __str__(self) -> str:
        rows = self.getRows()
        strlist = []
        for r in range(self.N):
            if self.rowLen(r) - rows[r] > 0:
                strlist.append(str(rows[r])+"{"+str(self.rowLen(r)-rows[r])+"}")
            else:
                strlist.append(str(rows[r]))
        return "("+ ", ".join(strlist)+")"

    def __repr__(self) -> str:
        return self.__str__()


if __name__ == "__main__":
    s1 = SU(3, multiplets=[(1,1)])
    s2 = SU(3, multiplets=[(1,1)])

    s3 = s1.combine(s2)

    print("{0} * {1} = {2}".format(s1.multiplets, s2.multiplets, s3.multiplets))
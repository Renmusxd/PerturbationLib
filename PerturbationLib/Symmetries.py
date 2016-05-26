
class Symmetry:
    """
    A class containing information as to how to
    combine and calculate representations.
    Multiplets are stored in (a,b,c...) notation to
    remove ambiguity. A conversion method will be implemented for
    each class.
    """
    def __init__(self, name, N, multiplets):
        """
        Create a symmetry
        :param name: unique identifier for symmetry
        :param multiplets: sum of multiplets (i.e. (1) + (3) + ...)
        """
        self.name = name
        self.N = N
        self.multiplets = list(multiplets)

    def combine(self, sym):
        """
        Combine a given U(N) symmetry with another U(N) symmetry
        :param sym: other U(N) with same N
        :return: new U(N) with sum of multiplets
        """
        if self.name==sym.name:
            mapRepr = set()
            for r1 in self.multiplets:
                for r2 in sym.multiplets:
                    mapRepr.add(self.combineRepr(r1, r2))
            return self.__class__(self.name, self.N, mapRepr)
        raise Exception("Symmetries do not match: "+self.name+", "+sym.name)

    def precombine(self, sym):
        """
        Combine a given U(N) symmetry with another U(N) symmetry
        :param sym: other U(N) with same N
        :return: new U(N) with sum of multiplets
        """
        if self.name==sym.name:
            mapRepr = set()
            for r1 in self.multiplets:
                for r2 in sym.multiplets:
                    mapRepr.add(self.combineRepr(r2, r1))
            return self.__class__(self.name, self.N, mapRepr)
        raise Exception("Symmetries do not match: "+self.name+", "+sym.name)

    def combineRepr(self, r1, r2):
        raise Exception("Method was not overridden")

    def containsSinglet(self):
        raise Exception("Method was not overridden")

    def singlet(self):
        '''
        Gives a singlet of the given group
        :return:
        '''
        raise Exception("Method was not overridden")

    def inverse(self):
        '''
        Gives the conjugate representation
        :return:
        '''
        raise Exception("Method was not overridden")


class U(Symmetry):
    """
    U(N) symmetry
    """
    def __init__(self, name, N, multiplets):
        super().__init__(name, N, multiplets)
        if N!=1:
            raise Exception("U(N>1) not yet implemented")

    def combineRepr(self, r1, r2):
        """
        Combines the multiplets for two U(N) symmetries, returns new multiplet
        :param r2:
        :return:
        """
        if self.N==1:
            # Return a tuple
            return r1[0] + r2[0],
        else:
            raise Exception("U(N>1) not yet implemented")

    def containsSinglet(self):
        if self.N==1:
            for m in self.multiplets:
                if m[0] == 0:
                    return True
        else:
            raise Exception("U(N>1) not yet implemented")

    def singlet(self):
        if self.N == 1:
            return U(self.name, self.N, [(0,)])
        else:
            raise Exception("U(N>1) not yet implemented")

    def inverse(self):
        if self.N == 1:
            mults = [(-x[0],) for x in self.multiplets]
            return U(self.name, self.N, mults)
        else:
            raise Exception("U(N>1) not yet implemented")

    def __repr__(self):
        return "U("+str(self.N)+"){"+(" + ".join([str(x) for x in self.multiplets]))+"}"

class SU(Symmetry):
    """
    SU(N)
    """
    def __init__(self, name, N, multiplets):
        super().__init__(name, N, multiplets)
        raise Exception("Not yet implemented")


class SO(Symmetry):
    """
    SO(N)
    """
    def __init__(self, name, N, multiplets):
        super().__init__(name, N, multiplets)
        raise Exception("Not yet implemented")

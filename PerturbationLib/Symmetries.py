from abc import ABCMeta, abstractmethod


class Symmetry(metaclass=ABCMeta):
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
        Combine a given symmetry with another symmetry
        :param sym: other U(N) with same N
        :return: new U(N) with sum of multiplets
        """
        if self.name==sym.name and self.N ==sym.N:
            mapRepr = set()
            for r1 in self.multiplets:
                for r2 in sym.multiplets:
                    mapRepr.update(self.combineRepr(r1, r2))
            return self.__class__(self.name, self.N, mapRepr)
        raise Exception("Symmetries do not match: "+self.name+", "+sym.name)

    @abstractmethod
    def combineRepr(self, r1, r2):
        '''
        Combines two representations of the symmetry.
        :return: [combined repr]
        '''
        raise Exception("Method was not overridden")

    @abstractmethod
    def containsSinglet(self):
        '''
        :return: true if repr contains a singlet
        '''
        raise Exception("Method was not overridden")

    @abstractmethod
    def singlet(self):
        '''
        Gives a singlet of the given group
        :return:
        '''
        raise Exception("Method was not overridden")

    @abstractmethod
    def inverse(self):
        '''
        Gives the conjugate representation
        :return:
        '''
        raise Exception("Method was not overridden")

    def __repr__(self):
        return "\{self.name\}"
    def _latex(self,*args):
        return self.__repr__()
    def _repr_latex(self):
        return "$"+self._latex()+"$"

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
        :param r1:
        :param r2:
        :return:
        """
        if self.N==1:
            return [(r1[0] + r2[0],)]
        else:
            raise Exception("U(N>1) not yet implemented")

    def containsSinglet(self):
        if self.N==1:
            return (0,) in self.multiplets
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
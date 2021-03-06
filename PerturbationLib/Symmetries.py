from abc import ABCMeta, abstractmethod
from typing import Sequence, Tuple


class Symmetry(metaclass=ABCMeta):
    """
    A class containing information as to how to
    combine and calculate representations.
    Multiplets are stored in (a,b,c...) notation to
    remove ambiguity. A conversion method will be implemented for
    each class.
    """
    def __init__(self, name: str, multiplets: Sequence[Tuple[int, ...]]):
        """
        Create a symmetry
        :param name: unique identifier for symmetry
        :param multiplets: sum of multiplets (i.e. (1) + (3) + ...)
        """
        self.name = name
        if multiplets is None:
            self.multiplets = [self.singletRepr()]
        else:
            self.multiplets = list(multiplets)

    def combine(self, sym: 'Symmetry') -> 'Symmetry':
        """
        Combine a given symmetry with another symmetry
        :param sym: other U(N) with same N
        :return: new U(N) with sum of multiplets
        """
        if self.matchesSymmetry(sym):
            mapRepr = []
            for r1 in self.multiplets:
                for r2 in sym.multiplets:
                    mapRepr += self.combineRepr(r1, r2)
            return self.constructWithNewRepr(mapRepr)
        raise Exception("Symmetries do not match: "+self.name+", "+sym.name)

    def singlet(self) -> 'Symmetry':
        """
        Gives a singlet of the given group
        :return:
        """
        return self.constructWithNewRepr([self.singletRepr()])

    def inverse(self) -> 'Symmetry':
        """
        Gives the conjugate representation
        :return:
        """
        return self.constructWithNewRepr(self.inverseRepr())

    def __call__(self, *args: Tuple[int, ...]) -> 'Symmetry':
        return self.constructWithNewRepr(args)

    @abstractmethod
    def constructWithNewRepr(self, newrepr: Sequence[Tuple[int, ...]]) -> 'Symmetry':
        """
        Create a new version of this symmetry with new representations.
        :param newrepr: new representations
        :return: Symmetry with new representations
        """
        raise Exception("Method was not overridden")

    @abstractmethod
    def combineRepr(self, r1: Tuple[int, ...], r2: Tuple[int, ...]) -> Tuple[int, ...]:
        """
        Combines two representations of the symmetry.
        :return: [combined repr]
        """
        raise Exception("Method was not overridden")

    @abstractmethod
    def containsSinglet(self) -> bool:
        """
        :return: true if repr contains a singlet
        """
        raise Exception("Method was not overridden")

    @abstractmethod
    def singletRepr(self) -> Tuple[int, ...]:
        """
        Gives a singlet representation of the given group
        :return:
        """
        raise Exception("Method was not overridden")

    @abstractmethod
    def inverseRepr(self) -> Sequence[Tuple[int, ...]]:
        """
        Gives the conjugate representation representation/multiplets
        :return:
        """
        raise Exception("Method was not overridden")

    @abstractmethod
    def matchesSymmetry(self, sym: 'Symmetry') -> bool:
        """
        True/False if symmetries may be combined.
        :param sym:
        :return:
        """
        return self.name == sym.name

    def __repr__(self):
        return self.name

    def _latex(self, *args):
        return self.__repr__()
    
    def _repr_latex(self):
        return "$"+self._latex()+"$"


class U(Symmetry):
    """
    U(N) symmetry
    """
    def __init__(self, N: int, name: str = None, multiplets: Sequence[Tuple[int, ...]] = None):
        self.N = N
        super().__init__(name or 'U('+str(N)+')', multiplets)
        if N != 1:
            raise Exception("U(N>1) not yet implemented")

    def constructWithNewRepr(self, newrepr: Sequence[Tuple[int, ...]]) -> 'U':
        """
        Create a new version of this symmetry with new representations.
        :param newrepr: new representations
        :return: Symmetry with new representations
        """
        return U(self.N, self.name, newrepr)

    def combineRepr(self, r1, r2) -> Sequence[Tuple[int, ...]]:
        """
        Combines the multiplets for two U(N) symmetries, returns new multiplet
        :param r1:
        :param r2:
        :return:
        """
        if self.N == 1:
            return [(r1[0] + r2[0],)]
        else:
            raise Exception("U(N>1) not yet implemented")

    def containsSinglet(self) -> bool:
        if self.N == 1:
            return (0,) in self.multiplets
        else:
            raise Exception("U(N>1) not yet implemented")

    def singletRepr(self) -> Tuple[int, ...]:
        if self.N == 1:
            return 0,
        raise Exception("U(N>1) not yet implemented")

    def inverseRepr(self) -> 'Sequence[Tuple[int]]':
        if self.N == 1:
            return [(-x[0],) for x in self.multiplets]
        else:
            raise Exception("U(N>1) not yet implemented")

    def matchesSymmetry(self, sym: 'Symmetry') -> bool:
        if isinstance(sym, U):
            return self.name == sym.name and self.N == sym.N
        return False

    def __repr__(self):
        return "U("+str(self.N)+"){"+(" + ".join([str(x) for x in self.multiplets]))+"}"

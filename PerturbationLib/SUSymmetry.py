from . import Symmetries


class SU(Symmetries.Symmetry):
    """
    SU(N)
    """
    def __init__(self, name, N, multiplets):
        super().__init__(name, N, multiplets)
        raise Exception("Not yet implemented")

    def combineRepr(self, r1, r2):
        """
        Combines two representations of the symmetry.
        :return: combined repr (may be a sum of repr)
        """
        raise Exception("Method was not overridden")

    def containsSinglet(self):
        """
        :return: true if repr contains a singlet
        """
        raise Exception("Method was not overridden")

    def singlet(self):
        """
        Gives a singlet of the given group
        :return:
        """
        raise Exception("Method was not overridden")

    def inverse(self):
        """
        Gives the conjugate representation
        :return:
        """
        raise Exception("Method was not overridden")
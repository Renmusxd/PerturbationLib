import numpy

from PerturbationLib.Feynman import Feynman
from PerturbationLib.Theory import Theory, Field
from PerturbationLib.Symmetries import U
from PerturbationLib.SUSymmetry import SU

if __name__ == "__main__":
    usym = U(1)
    susym = SU(2)

    t = Theory(usym, susym, trunc=4)
    f1 = Field("\phi", usym((1,)), susym((0,)))
    t.addField(f1)
    f2 = Field("\psi", usym((2,)), susym((0,)))
    t.addField(f2)
    f3 = Field("\zeta", usym((0,)), susym((0,)))
    t.addField(f3)

    print(t.getInt())

    startd = {("\psi", False): 1}
    endd = {("\phi", False): 2}

    feyn = Feynman(t)
    print(feyn)
    # paths = feyn.listPaths(startd, endd, maxorder=4)
    #
    # endvec = feyn.convertDictToVec(endd)
    #
    # for path in paths:
    #     start = feyn.convertDictToVec(startd)
    #     print(path)
    #     for inter in path:
    #         start += inter[1] + inter[2]
    #         print(inter[1], inter[2], start)
    #     assert(numpy.array_equal(start, endvec))

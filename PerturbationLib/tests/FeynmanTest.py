import numpy

from PerturbationLib.Feynman import Feynman
from PerturbationLib.Theory import Theory, Field
from PerturbationLib.Symmetries import U
from PerturbationLib.SUSymmetry import SU

if __name__ == "__main__":
    t = Theory([U(1),SU(2)], trunc=4)
    f1 = Field("\phi",[U(1, multiplets=[(1,)]), SU(2, multiplets=[(0,)])])
    t.addField(f1)
    f2 = Field("\psi",[U(1, multiplets=[(2,)]), SU(2, multiplets=[(0,)])])
    t.addField(f2)
    f3 = Field("\zeta",[U(1, multiplets=[(0,)]), SU(2, multiplets=[(0,)])])
    t.addField(f3)

    print(t.getInt())

    startd = {("\psi",False):1}
    endd = {("\phi",False):2}

    feyn = Feynman(t)
    paths = feyn.listPaths(startd,endd,maxorder=4)

    endvec = feyn.convertDictToVec(endd)

    for path in paths:
        start = feyn.convertDictToVec(startd)
        print(path)
        for inter in path:
            start += inter[1] + inter[2]
            print(inter[1],inter[2],start)
        assert(numpy.array_equal(start,endvec))
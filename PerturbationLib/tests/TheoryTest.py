from PerturbationLib.Symmetries import U
from PerturbationLib.Theory import Theory, Field
from PerturbationLib.Utilities import *


if __name__ == "__main__":
    usym = U(1)
    t = Theory([usym], trunc=4)
    f = Field("\phi", usym((1,)))
    t.addField(f)
    f = Field("\psi", usym((-3,)))
    t.addField(f)
    print(t)

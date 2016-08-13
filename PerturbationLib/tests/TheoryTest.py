from PerturbationLib.Symmetries import U
from PerturbationLib.Theory import Theory, Field
from PerturbationLib.Utilities import *

if __name__ == "__main__":
    t = Theory([U(1)], trunc=4)
    f = Field("\phi",[U(1,multiplets=[(1,)])])
    t.addField(f)
    f = Field("\psi",[U(1,multiplets=[(-3,)])])
    t.addField(f)
    print(t)
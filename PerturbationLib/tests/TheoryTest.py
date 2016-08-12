try:
    from ..Symmetries import U
    from ..Theory import Theory, Field
    from ..Utilities import *
except:
    from Symmetries import U
    from Theory import Theory, Field
    from Utilities import *

print(weightedTruncatedPowerset([("a",2),("b",2)],4))

if __name__ == "__main__":
    t = Theory([U("charge",1,[])], trunc=4)
    f = Field("\phi",[U("charge",1,[(1,)])])
    t.addField(f)
    f = Field("\psi",[U("charge",1,[(-3,)])])
    t.addField(f)
    print(t)
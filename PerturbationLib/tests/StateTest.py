try:
    from ..ParticleState import *
except:
    from ParticleState import *

def A():
    a = SingleState([1,0])
    b = SingleState([2,0])
    c = SingleState([0,1])
    d = a+b+c+c
    e = d.transpose() * c
    assert e==2
    f = ( AOp(1) * c)
    g = f.transpose() * SingleState([0,0])
    assert g==1

def B():
    d = ( XOp(1) * SingleState([1,0]) )
    e = ( AOp(1) * d ) + ( COp(1) * d )

def C():
    a_, a = COp(0), AOp(0)
    b = SingleState([0,0])
    assert 2 == (a_ * a_ * b).transpose() * (a_ * (a_ + a) * b)
from PerturbationLib.ParticleState import *


def A():
    a = SingleState([1,0])
    b = SingleState([2,0])
    c = SingleState([0,1])
    d = a+b+c+c
    e = d.transpose() * c
    f = ( XOp(1) * c)

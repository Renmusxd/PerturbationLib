"""
A library for SHO level perturbation theory. Mostly for keeping track of constants.
"""

import itertools
import sympy
from sympy import latex
from typing import Iterable


class Braket:
    BRA = True
    KET = False

class SingleState:
    def __init__(self, particles: Iterable[int], multiplier: float = 1.0, bk=Braket.KET):
        self.particles = list(particles)
        self.mult = multiplier
        self.bk = bk

    def __eq__(self, other):
        '''
        Same state, disregard multiplier
        :param other: other SingleState
        :return: True/False
        '''
        if other.__class__==self.__class__:
            return canAddSingleStates(self,other)
        return False

    def __ne__(self, other):
        '''
        Not __eq__
        :param other: other SingleState
        :return: False/True
        '''
        return not self.__eq__(other)

    def __rmul__(self, other):
        '''
        :param other: other SingleState or float/int
        :return:
        '''
        if other.__class__==SingleState:
            if other.bk==Braket.BRA and self.bk==Braket.KET:
                matching = True
                for i in range(len(self.particles)):
                    if self.particles[i]!=other.particles[i]:
                        matching = False
                        break
                if matching:
                    return other.mult * self.mult
            else:
                raise Exception("We do not yet support tensor multiplication")
        elif other.__class__==State:
            return sum( [self * state for state in other.states ] )
        else:
            try:
                s = self.copy()
                s.mult *= other
                return s
            except:
                return other.__mul__(self)

    def __mul__(self, other):
        '''
        Inner product of two states
        :param other:
        :return: int
        '''
        if other.__class__==SingleState:
            if self.bk==Braket.BRA and other.bk==Braket.KET:
                matching = True
                for i in range(len(self.particles)):
                    if self.particles[i]!=other.particles[i]:
                        matching = False
                        break
                if matching:
                    return other.mult * self.mult
                return 0
            else:
                raise Exception("We do not yet support tensor multiplication")
        elif other.__class__==State:
            return sum( [self * state for state in other.states ] )
        else:
            s = self.copy()
            s.mult *= other
            return s
        raise TypeError("unsupported operand type(s)")

    def __add__(self, other):
        return addStates(self, other)

    def copy(self):
        return SingleState(self.particles, self.mult)

    def transpose(self):
        '''
        Takes transpose of the state
        :return: conjugate of multiplier and ket/bra of state
        '''
        s = self.copy()
        s.bk = not s.bk
        s.mult = sympy.conjugate(s.mult)
        return s

    def __str__(self):
        mstr = "" if self.mult==1 else latex(self.mult)
        if self.bk == Braket.KET:
            return mstr+"|"+",".join([str(i) for i in self.particles])+">"
        elif self.bk == Braket.BRA:
            return mstr+"<"+",".join([str(i) for i in self.particles])+"|"
    def __repr__(self):
        return self.__str__()
    def _repr_latex_(self):
        return "$"+self._latex()+"$"
    def _latex(self, *args):
        return "$"+str(self)+"$"

class State:
    def __init__(self, states):
        '''
        Creates a new state
        init([1,0,0]) -> |1,0,0>
        :param numparticles: number of particles available
        :param mult: multiplier
        :param states: array of SingleStates: |a> + |b> + ...
        '''
        self.states = [state.copy() for state in states if abs(state.mult)>0]

    @classmethod
    def initFromList(cls, particleStates, bk=Braket.KET):
        '''
        State constructor
        :param particleStates: Array of particle states
        :return:
        '''
        c = cls([])
        c.states = cls( [SingleState(ps, bk=bk) for ps in particleStates] )
        return c

    def __add__(self, state):
        return addStates(self, state)

    def __rmul__(self, other):
        if other.__class__==SingleState or other.__class__==State:
            total = 0
            # For each of local states
            for m_s in self.states:
                # For each other state
                total = total + (other * self)
            return total
        else:
            s = self.copy()
            for s_s in s.states:
                s_s.mult *= other
            return s


    def __mul__(self, other):
        '''
        Applied with: self * state
        Not: state * self
        :param state: state to which to be applied
        :return: inner product with state
        '''
        if other.__class__==SingleState or other.__class__==State:
            total = 0
            # For each of local states
            for m_s in self.states:
                # For each other state
                total = total + (m_s * other)
            return total
        else:
            s = self.copy()
            for s_s in s.states:
                s_s.mult *= other
            return s

    def __eq__(self, other):
        if other.__class__==State:
            for i in range(len(self.states)):
                if self.states[i]!=other.states[i]:
                    return False
            return True
        return False

    def __neg__(self):
        newstate = self.copy()
        for i in range(len(self.states)):
            newstate.states[i].mult = -self.states[i].mult
        return newstate

    def transpose(self):
        '''
        Takes transpose of each state
        :return: conjugate of multipliers and ket/bra of states
        '''
        newstate = self.copy()
        for i in range(len(self.states)):
            newstate.states[i] = newstate.states[i].transpose()
        return newstate

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        strlist = [str(self.states[i]) for i in range(len(self.states))]
        return " + ".join(strlist)
    def _repr_latex_(self):
        return "$"+self._latex()+"$"
    def _latex(self, *args):
        return "$"+str(self)+"$"

    def copy(self):
        return State(self.states)

class Operator:
    mult = 1

class XOp(Operator):
    '''
    x = c(a + a')
    '''
    X_CONST = 1
    def __init__(self, index,mult=1):
        self.index = index
        self.mult = XOp.X_CONST * mult

    def __call__(self, state):
        return self.__mul__(state)

    def __mul__(self, other):
        if other.__class__==SingleState or other.__class__==State:
            return (COp(self.index) * other) + (AOp(self.index) * other)
        elif issubclass(other.__class__,Operator):
            # a * O -> [a O]
            return OpProduct([self, other])
        new = XOp(self.index)
        new.mult *= other
        return new

    def __add__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpSum:
                return OpSum([self] + other.ops)
            else:
                return OpSum([self] + [other])
        else:
            raise TypeError("Cannot add operator to anything except operators")

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        mstr = "" if (self.mult==1) else latex(self.mult)
        return mstr + "x_" + str(self.index)
    def _repr_latex_(self):
        return "$"+self._latex()+"$"
    def _latex(self, *args):
        mstr = "" if (self.mult==1) else latex(self.mult)
        return mstr + "x_{" + str(self.index)+"}"
    def copy(self):
        return XOp(self.index,self.mult)

class COp(Operator):
    def __init__(self, index,mult=1):
        self.index = index
        self.mult = mult

    def __call__(self, state):
        return self.__mul__(state)

    def __mul__(self, other):
        if other.__class__==State:
            other = other.copy()
            for i in range(len(other.states)):
                other.states[i].particles[self.index] += 1
                other.states[i].mult *= sympy.sqrt(other.states[i].particles[self.index])
            other.states = [s for s in other.states if abs(s.mult) > 0]
            return other
        elif other.__class__==SingleState:
            other = other.copy()
            other.particles[self.index] += 1
            other.mult *= sympy.sqrt(other.particles[self.index])
            return other
        elif issubclass(other.__class__,Operator):
            # a * O -> [a O]
            return OpProduct([self, other])
        new = self.copy()
        new.mult *= other
        return new

    def __add__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpSum:
                return OpSum([self] + other.ops)
            else:
                return OpSum([self] + [other])
        else:
            raise TypeError("Cannot add operator to anything except operators")

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        mstr = "" if (self.mult==1) else latex(self.mult)
        return mstr+"a'_"+str(self.index)
    def _repr_latex_(self):
        return "$"+self._latex()+"$"
    def _latex(self, *args):
        mstr = "" if (self.mult==1) else latex(self.mult)
        return mstr+"a^{\dag}_{"+str(self.index)+"}"
    def copy(self):
        return COp(self.index,self.mult)

class AOp(Operator):
    def __init__(self, index, mult=1):
        self.index = index
        self.mult = mult

    def __call__(self, state):
        return self.__mul__(state)

    def __mul__(self, other):
        if other.__class__==State:
            other = other.copy()
            for i in range(len(other.states)):
                if other.states[i].particles[self.index]<=0:
                    other.states[i].mult = 0
                else:
                    other.states[i].mult *= sympy.sqrt(other.states[i].particles[self.index])
                    other.states[i].particles[self.index] -= 1
            other.states = [s for s in other.states if abs(s.mult) > 0]
            return other
        elif other.__class__==SingleState:
            other = other.copy()
            if other.particles[self.index]>0:
                other.mult *= sympy.sqrt(other.particles[self.index])
                other.particles[self.index] -= 1
                return other
            else:
                return State([])
        elif issubclass(other.__class__,Operator):
            # a * O -> [a O]
            return OpProduct([self, other])
        new = self.copy()
        new.mult *= other
        return new

    def __add__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpSum:
                return OpSum([self] + other.ops)
            else:
                return OpSum([self] + [other])
        else:
            raise TypeError("Cannot add operator to anything except operators")

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        mstr = "" if (self.mult==1) else latex(self.mult)
        return mstr+"a_"+str(self.index)
    def _repr_latex_(self):
        return "$"+self._latex()+"$"
    def _latex(self, *args):
        mstr = "" if (self.mult==1) else latex(self.mult)
        return mstr+"a_{"+str(self.index)+"}"
    def copy(self):
        return AOp(self.index,self.mult)


class OpProduct(Operator):
    def __init__(self, ops, mult=1):
        self.ops = [op.copy() for op in ops]
        self.mult = mult
        for op in self.ops:
            self.mult *= op.mult
            op.mult = 1
    def __rmul__(self, other):
        if issubclass(other.__class__, Operator):
            if other.__class__==OpProduct:
                return OpProduct(other.ops + self.ops)
            else:
                return OpProduct([other] + self.ops)
        else:
            new = self.copy()
            new.mult *= other
            return new
    def __mul__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpProduct:
                return OpProduct(self.ops + other.ops)
            else:
                return OpProduct(self.ops + [other])
        elif other.__class__ == SingleState or other.__class__ == State:
            acc = other
            for op in reversed(self.ops):
                acc = op * acc
            return acc
        else:
            new = self.copy()
            new.mult *= other
            return new

    def __add__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpSum:
                return OpSum([self] + other.ops)
            else:
                return OpSum([self] + [other])
        else:
            raise TypeError("Cannot add operator to anything except operators")

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        mstr = "" if (self.mult == 1) else latex(self.mult)
        return mstr+"(" + (" ".join([str(op) for op in self.ops])) + ")"
    def _repr_latex_(self):
        return "$"+self._latex()+"$"
    def _latex(self, *args):
        mstr = "" if (self.mult == 1) else latex(self.mult)
        return mstr+"(" + (" ".join([op._latex() for op in self.ops])) + ")"
    def copy(self):
        return OpProduct(self.ops, self.mult)

class OpSum(Operator):
    def __init__(self, ops, mult=1):
        self.ops = [op.copy() for op in ops]
        self.mult = mult
    def __rmul__(self, other):
        if issubclass(other.__class__, Operator):
            if other.__class__==OpProduct:
                return OpProduct(other.ops + [self])
            else:
                return OpProduct([other] + [self])
        else:
            new = self.copy()
            new.mult *= other
            return new
    def __mul__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpProduct:
                return OpProduct([self] + other.ops)
            else:
                return OpProduct([self] + [other])
        elif other.__class__ == SingleState or other.__class__ == State:
            acc = State([])
            for op in self.ops:
                acc += op * other * self.mult
            return acc
        else:
            new = self.copy()
            new.mult *= other
            return new

    def __add__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpSum:
                # Add mult to each child
                if self.mult!=1:
                    for op in self.ops:
                        op.mult *= self.mult
                    self.mult = 1
                return OpSum(self.ops + other.ops)
            else:
                return OpSum(self.ops + [other])
        else:
            raise TypeError("Cannot add operator to anything except operators")

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        mstr = "" if (self.mult == 1) else latex(self.mult)
        return mstr + "(" + (" + ".join([str(op) for op in self.ops])) + ")"
    def _repr_latex_(self):
        return "$"+self._latex()+"$"
    def _latex(self, *args):
        mstr = "" if (self.mult == 1) else latex(self.mult)
        return mstr + "(" + (" + ".join([op._latex() for op in self.ops])) + ")"

    def copy(self):
        return OpSum(self.ops, self.mult)

class DeltaH(Operator):
    def __init__(self, ladders, mult=1):
        '''
        Creates a perturbation operator
        POperator([1,2]) -> x.y^2
        :param ladders: position operators for perturbation
        :param mult: multiplier of perturbation 'lambda'
        '''
        self.pert = ladders[:]

    def __call__(self, state):
        '''
        Corresponds to applying operator to state
        :param state: state to which to be applied
        :return: set of state outputs
        '''
        return self.__mul__(state)

    def __mul__(self, state):
        '''
        Corresponds to applying operator to state
        :param state: state to which to be applied
        :return: set of state outputs
        '''
        # Creation:         a'|n> = sqrt(n+1).|n+1>
        # Annihilation:      a|n> = sqrt(n).|n-1>
        output_states = []
        # For each particle
        for l_indx in range(len(self.pert)):
            # For each possible number of creation (implicit and annihilation) operators
            for ladder_value in range(self.pert[l_indx]):
                # Create list of operators
                used_ops = [COp(l_indx)]*ladder_value + [AOp(l_indx)]*(self.pert[l_indx]-(ladder_value+1))
                print(used_ops)
                perms = itertools.permutations(used_ops)
                # For every possible arrangement of the operators
                acc_state = state.copy()
                for perm in perms:
                    sel_state = state.copy()
                    for op in perm:
                        sel_state = op * sel_state
                    print(sel_state)
                    acc_state.mult += sel_state.mult
                output_states.append(acc_state)
        return output_states


# ===================================== #
# Helper functions for above operations # =====================================
# ===================================== #


def canAddSingleStates(state1, state2):
    if len(state1.particles)==len(state2.particles):
        if state1.bk!=state2.bk:
            return False
        for i in range(len(state1.particles)):
            if state1.particles[i]!=state2.particles[i]:
                return False
        return True
    return False


def addSingleStates(state1, state2):
    if canAddSingleStates(state1,state2):
        return SingleState(state1.particles,
                           state1.mult + state2.mult,
                           state1.bk)
    else:
        return State([state1, state2])


def addSingleStateToState(state, singlestate):
    newstates = [s.copy() for s in state.states]
    added = False
    for i in range(len(newstates)):
        if canAddSingleStates(singlestate,newstates[i]):
            newstates[i] = addSingleStates(singlestate, newstates[i])
            added = True
            break
    if not added:
        newstates.append(singlestate.copy())
    return State(newstates)


def addStateToState(state1, state2):
    newstate = state1.copy()
    for s in (s.copy() for s in state2.states):
        newstate = newstate + s
    return newstate


def addStates(state1, state2):
    if state1.__class__ == State:
        if state2.__class__ == State:
            return addStateToState(state1, state2)
        else:
            return addSingleStateToState(state1, state2)
    elif state2.__class__ == State:
        if state1.__class__ == State:
            return addStateToState(state2, state1)
        else:
            return addSingleStateToState(state2, state1)
    elif state1.__class__ == SingleState and state2.__class__ == SingleState:
        return addSingleStates(state1, state2)
    else:
        TypeError("Cannot add "+str(state1)+" and "+str(state2))
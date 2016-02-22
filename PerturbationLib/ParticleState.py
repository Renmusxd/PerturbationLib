import math
import itertools
import sympy
from sympy import Symbol

class braket:
    BRA = True
    KET = False

class SingleState:
    def __init__(self, particles, mult=1, bk=braket.KET):
        self.particles = particles[:]
        self.mult = mult
        self.bk = bk

    def __eq__(self, other):
        '''
        Same state, disregard multiplier
        :param other: other SingleState
        :return: True/False
        '''
        if other.__class__==self.__class__ and self.bk==other.bk:
            for i in xrange(len(self.particles)):
                if self.particles[i]!=other.particles[i]:
                    return False
            return True
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
            if other.bk==braket.BRA and self.bk==braket.KET:
                matching = True
                for i in xrange(len(self.particles)):
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
            if self.bk==braket.BRA and other.bk==braket.KET:
                matching = True
                for i in xrange(len(self.particles)):
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
        if self==other:
            return SingleState(self.particles, mult=self.mult+other.mult, bk=self.bk)
        elif other.__class__==SingleState:
            s = other.copy()
            if s==self:
                s.mult += self.mult
                return s
            else:
                return State([self, s])
        elif other.__class__==State:
            s = other.copy()
            added = False
            for s_s in s.states:
                if s_s==self:
                    s_s.mult += self.mult
                    added = True
                    break
            if not added:
                s.states.append(self.copy())
            return s
        else:
            return other.__add__(self)


    def copy(self):
        return SingleState(self.particles, self.mult)

    def transpose(self):
        s = self.copy()
        s.bk = not s.bk
        return s

    def __str__(self):
        mstr = "" if self.mult==1 else str(self.mult)
        if self.bk == braket.KET:
            return mstr+"|"+",".join([str(i) for i in self.particles])+">"
        elif self.bk == braket.BRA:
            return mstr+"<"+",".join([str(i) for i in self.particles])+"|"
    def __repr__(self):
        return self.__str__()

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
    def initFromList(cls, particleStates, bk=braket.KET):
        '''
        State constructor
        :param particleStates: Array of particle states
        :return:
        '''
        c = cls([])
        c.states = cls( [SingleState(ps, bk=bk) for ps in particleStates] )
        return c

    def __add__(self, state):
        if state.__class__==SingleState:
            return state + self
        elif state.__class__==State:
            retstate = self.copy()
            for i in xrange(len(state.states)):
                i_was_added = False
                for j in xrange(len(retstate.states)):
                    if retstate.states[j]==state.states[i]:
                        retstate.states[j].mult += state.states[i].mult
                        i_was_added = True
                        break
                if not i_was_added:
                    retstate.states.append(state.states[i].copy())

            return retstate

    def __rmul__(self, other):
        try:
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
        except:
            return other.__mul__(self)


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
            for i in xrange(len(self.states)):
                if self.states[i]!=other.states[i]:
                    return False
            return True

        return False

    def __neg__(self):
        newstate = self.copy()
        for i in xrange(len(self.states)):
            newstate.states[i].mult = -self.states[i].mult
        return newstate

    def transpose(self):
        newstate = self.copy()
        for i in xrange(len(self.states)):
            newstate.states[i] = newstate.states[i].transpose()
        return newstate

    def __repr__(self):
        strlist = [str(self.states[i]) for i in xrange(len(self.states))]
        return " + ".join(strlist)

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
        mstr = "" if (self.mult==1) else str(self.mult)
        return mstr + "x_" + str(self.index)

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
            for i in xrange(len(other.states)):
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
        mstr = "" if (self.mult==1) else str(self.mult)
        return mstr+"a'_"+str(self.index)

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
            for i in xrange(len(other.states)):
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
        mstr = "" if (self.mult==1) else str(self.mult)
        return mstr+"a_"+str(self.index)

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
        mstr = "" if (self.mult==1) else str(self.mult)
        return mstr+"(" + (" ".join([str(op) for op in self.ops])) + ")"

    def copy(self):
        return OpProduct(self.ops,self.mult)

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
                acc += op * other
            return acc
        else:
            new = self.copy()
            new.mult *= other
            return new

    def __add__(self, other):
        if issubclass(other.__class__,Operator):
            if other.__class__ == OpSum:
                return OpSum(self.ops + other.ops)
            else:
                return OpSum(self.ops + [other])
        else:
            raise TypeError("Cannot add operator to anything except operators")

    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        return "(" + (" + ".join([str(op) for op in self.ops])) + ")"

    def copy(self):
        return OpSum(self.ops,self.mult)

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
        for l_indx in xrange(len(self.pert)):
            # For each possible number of creation (implicit and annihilation) operators
            for ladder_value in xrange(self.pert[l_indx]):
                # Create list of operators
                used_ops = [COp(l_indx)]*ladder_value + [AOp(l_indx)]*(self.pert[l_indx]-(ladder_value+1))
                print used_ops
                perms = itertools.permutations(used_ops)
                # For every possible arrangement of the operators
                acc_state = state.copy()
                for perm in perms:
                    sel_state = state.copy()
                    for op in perm:
                        sel_state = op(sel_state)
                    print sel_state
                    acc_state.mult += sel_state.mult
                output_states.append(acc_state)
        return output_states
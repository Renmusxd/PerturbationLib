import math
import numpy


def truncatedPowerset(values,trunc,depth=0):
    if len(values)==0 or depth>=trunc:
        return [[]]
    elif len(values)==1:
        h = values[0]
        return [[h],[]]
    elif len(values)>1:
        h, t = values[0],values[1:]
        woh = truncatedPowerset(t,trunc,depth)
        wh =[[h] + l for l in truncatedPowerset(t,trunc,depth+1)]
        return wh + woh


def truncCombinations(values,trunc):
    return weightedTruncatedPowerset(values,trunc)[:-1]


def weightedTruncatedPowerset(values,trunc,depth=0):
    '''
    Acts like the above, but groups values into equivalence classes of appropriate sizes
    and returns the powerset with that modulo
    :param values: tuples of (obj,weight)
    :param trunc: max size of set in powerset
    :return: powerset mod obj
    '''
    if len(values)==0 or depth>=trunc:
        return [[]]
    elif len(values)==1:
        (h,hw) = values[0]
        return [[h]*i for i in reversed(range(min(hw+1,trunc-depth+1)))]
    elif len(values)>1:
        (h,hw), t = values[0], values[1:]
        woh = weightedTruncatedPowerset(t, trunc, depth)

        wht = t
        if hw>1:
            wht = [(h,hw-1)] + t
        wh = [[h] + l for l in weightedTruncatedPowerset(wht,trunc,depth+1)]

        return wh + woh


def genInOutPairs(fieldvec):
    """
    Generates all pairs of vectors (i,o) such that i + o = fieldvec
    :param fieldvec:
    :yield: tuples (i,o)
    """
    fi = numpy.array([0 for _ in range(len(fieldvec))])
    fo = numpy.array(fieldvec)
    fb = [False for _ in range(len(fieldvec))]

    i = 0

    yield numpy.array(fi), numpy.array(fo)
    while i < len(fieldvec):
        if fo[i] > 0:
            fo[i] -= 1
            fi[i] += 1
            i = 0
            yield numpy.array(fi), numpy.array(fo)
        else:
            fo[i] = fi[i]
            fi[i] = 0
            i += 1


def genParticlePerms(fieldvec):
    s = len(fieldvec)
    if s%2 != 0:
        raise Exception("Must have field and anti-field slot for each field")
    perms = list(genInOutPairs(fieldvec))[1:-1]
    for fi, fo in perms:
        for i in range(0, s, 2):
            fo[i], fo[i + 1] = fo[i + 1], fo[i]
    return perms


def swapAntis(vec):
    vec = numpy.array(vec)
    for i in range(0, len(vec), 2):
        vec[i], vec[i + 1] = vec[i + 1], vec[i]
    return vec


def chooseNFromM(n,m):
    if n == 0:
        yield [0 for _ in range(m)]
    elif m == 1:
        yield [n]
    else:
        for i in range(n+1):
            for tail in chooseNFromM(n-i,m-1):
                yield [i] + tail


def vectorPath(start, end, nodes, trunc):
    result = []
    if trunc == 0:
        return result
    for i in range(len(nodes)):
        req, prod = nodes[i]
        s1 = numpy.array(start + req)
        for x in s1:
            # If a single element is negative, this is an invalid path
            if x < 0:
                break
        else:
            s = s1 + prod
            if numpy.array_equal(s, end):
                result.append([i])
            subpaths = vectorPath(s,end,nodes,trunc-1)
            result += [[i] + path for path in subpaths]
    return result


if __name__ == "__main__":
    l = list(chooseNFromM(3,3))
    print(len(l),l)
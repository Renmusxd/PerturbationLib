import math
import numpy
from functools import reduce
from typing import Sequence, TypeVar, List, Tuple, Generator, Mapping

T = TypeVar('T')


def truncatedPowerset(values: Sequence[T], trunc: int, depth: int = 0) -> List[List[T]]:
    """
    Get the powerset of values up to order trunc
    :param values: set of values from which to make powerset
    :param trunc: maximum size of powerset items
    :param depth: internal variable to subtract from trunc
    :return:
    """
    if len(values) == 0 or depth >= trunc:
        return [[]]
    elif len(values) == 1:
        h = values[0]
        return [[h], []]
    elif len(values) > 1:
        h, t = values[0], values[1:]
        woh = truncatedPowerset(t, trunc, depth)
        wh = [[h] + l for l in truncatedPowerset(t, trunc, depth+1)]
        return wh + woh


def truncCombinations(values: Sequence[Tuple[T, int]], trunc: int) -> List[List[T]]:
    return weightedTruncatedPowerset(values, trunc)[:-1]


def weightedTruncatedPowerset(values: Sequence[Tuple[T, int]], trunc: int, depth: int = 0) -> List[List[T]]:
    '''
    Acts like the above, but groups values into equivalence classes of appropriate sizes
    and returns the powerset with that modulo
    :param values: tuples of (obj,weight)
    :param trunc: max size of set in powerset
    :param depth: internally used variable to subtract from trunc
    :return: powerset mod obj
    '''
    if len(values) == 0 or depth >= trunc:
        return [[]]
    elif len(values) == 1:
        (h, hw) = values[0]
        return [[h]*i for i in reversed(range(min(hw+1, trunc-depth+1)))]
    elif len(values) > 1:
        (h, hw), t = values[0], values[1:]
        woh = weightedTruncatedPowerset(t, trunc, depth)

        wht = t
        if hw > 1:
            wht = [(h, hw-1)] + t
        wh = [[h] + l for l in weightedTruncatedPowerset(wht, trunc, depth+1)]

        return wh + woh


def genInOutPairs(fieldvec: Sequence[int],
                  swapvec: Sequence[int] = None) -> Generator[Tuple[numpy.ndarray, numpy.ndarray], None, None]:
    """
    Generates all pairs of vectors (i,o) such that i[k] + o[swapvec[k]] = fieldvec
    :param fieldvec: list of particles
    :param swapvec: list of places to swap
    :yield: tuples (i,o)
    """
    if swapvec is None:
        swapvec = list(range(len(fieldvec)))

    fi = numpy.array([0 for _ in range(len(fieldvec))])
    fo = numpy.array(fieldvec)

    i = 0
    yield fi, fo
    while i < len(fieldvec):
        if fo[swapvec[i]] > 0:
            fo[swapvec[i]] -= 1
            fi[i] += 1
            i = 0
            yield fi, fo
        else:
            fo[swapvec[i]] = fi[i]
            fi[i] = 0
            i += 1


def swapAntis(vec: Sequence[int]) -> Sequence[int]:
    vec = numpy.array(vec)
    for i in range(0, len(vec), 2):
        vec[i], vec[i + 1] = vec[i + 1], vec[i]
    return vec


def chooseNFromM(n: int, m: int) -> Generator[List[int], None, None]:
    if n == 0:
        yield [0 for _ in range(m)]
    elif m == 1:
        yield [n]
    else:
        for i in range(n+1):
            for tail in chooseNFromM(n-i, m-1):
                yield [i] + tail


def genVectorPaths(start: Sequence[int], end: Sequence[int],
                   nodes: Mapping[T, Tuple[Sequence[int], Sequence[int]]],
                   trunc: int) -> Generator[List[Tuple[T, Sequence[int], Sequence[int]]], None, None]:
    start, end = numpy.asarray(start), numpy.asarray(end)
    if trunc == 0:
        yield []

    for node in nodes:
        node_tuple = (node, nodes[node][0], nodes[node][1])
        required_particles, produced_particles = numpy.asarray(nodes[node][0]), numpy.asarray(nodes[node][1])
        comb = start + required_particles
        if not numpy.all(numpy.positive(comb)):
            continue
        after_inter = comb + produced_particles
        if numpy.array_equal(after_inter, end):
            yield [node_tuple]
        else:
            for subpath in genVectorPaths(start, end, nodes, trunc - 1):
                yield [node_tuple] + subpath


if __name__ == "__main__":
    for x in genInOutPairs([1,1,1], swapvec=[1,0,2]):
        print(x)

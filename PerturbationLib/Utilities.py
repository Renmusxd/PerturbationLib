import math

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
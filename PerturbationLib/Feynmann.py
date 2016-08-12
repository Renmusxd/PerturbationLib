"""
A package for keeping track of energy flow and so forth
"""

class Interaction:
    def __init__(self, theory):
        pass

    def amplitude(self, endstate):
        '''
        Calculates the total amplitude of the endstate
        :param endstate: State of interest
        :return: A
        '''
        pass

    def listPaths(self, endstate):
        '''
        Lists each path to get to given endstate (endstate * state > 0)
        :param endstate: State of interest
        :return: [paths]
        '''
        pass
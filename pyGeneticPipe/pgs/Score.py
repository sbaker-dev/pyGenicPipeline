from pyGeneticPipe.core.Input import Input


class Score(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_pgs(self):
        return
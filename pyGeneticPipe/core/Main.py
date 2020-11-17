from pyGeneticPipe.shell.ShellMaker import ShellMaker
from pyGeneticPipe.clean.Cleaner import Cleaner
from pyGeneticPipe.core.Input import Input


class Main(ShellMaker, Cleaner, Input):
    def __init__(self, args):
        """
        This Class inherits all other classes that can be used, and then execute the job via getattr
        :param args: Json args that have been set via GUI or manually
        """
        super().__init__(args)

        getattr(Main, self.operation)(self)



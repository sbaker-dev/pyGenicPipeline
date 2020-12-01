from pyGeneticPipe.support.ShellMaker import ShellMaker
from pyGeneticPipe.utils.misc import terminal_time
from pyGeneticPipe.pgs.Cleaner import Cleaner
from pyGeneticPipe.core.Input import Input
from colorama import init


class Main(ShellMaker, Cleaner, Input):
    def __init__(self, args):
        """
        This Class inherits all other classes that can be used, and then execute the job via getattr
        :param args: Json args that have been set via GUI or manually
        """
        init(autoreset=True)
        super().__init__(args)
        print(f"Starting {self.operation}: {terminal_time()}")
        getattr(Main, self.operation)(self)

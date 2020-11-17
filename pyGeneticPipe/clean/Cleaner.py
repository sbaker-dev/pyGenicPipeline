from pyGeneticPipe.utils.misc import directory_iterator
from pyGeneticPipe.core.Input import Input
from pathlib import Path


class Cleaner(Input):
    def __init__(self, args):
        super().__init__(args)

    def create_validation_group(self):
        print("Hello from cleaner")

        # Create the validation group
        validation_group = self.project_file.create_group(self.h5_validation)

        # Create a dataset of all the chromosomes we have to work with
        self._validation_chromosomes(validation_group)

    def _validation_chromosomes(self, validation_group):
        """
        This will create a dataset of all the chromosomes that we have to work with with our validation group in the
        h5py file
        """
        valid_chromosomes = []
        for file in directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:
                valid_chromosomes.append(int(Path(self.load_directory, file).stem.split("_")[-1]))
        valid_chromosomes.sort()

        validation_group.create_dataset(self.h5_valid_chromosome, data=valid_chromosomes)

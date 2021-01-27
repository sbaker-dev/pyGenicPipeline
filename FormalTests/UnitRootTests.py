from pyGenicPipeline.Processors import *

from miscSupports import load_yaml
from pathlib import Path
import subprocess
import unittest


class MyTestCase(unittest.TestCase):

    @staticmethod
    def _write_path():
        """Set the write path"""
        return Path(Path(__file__).parent, "TestDirectory")

    @staticmethod
    def _config_loader():
        """Some processes like plink need their paths/attributes set externally from a yaml config file"""
        return load_yaml(Path(Path(__file__).parent.absolute(), "config.yaml"))

    @staticmethod
    def _file_exists_cleanup(file_path):
        """Some tests we simply want to validate a file exists, what is within them """

        assert file_path.exists(), f"plink failed to produce {file_path.name}"
        print(f"Found {file_path.name}: Deleting")
        file_path.unlink()

    def test_plink(self):
        """
        Test Plink 1.9 installed and can run
        """
        plink_path = Path(self._config_loader()["plink_1_path"])
        write_path = self._write_path()
        assert (plink_path.exists() and write_path.exists())

        job = subprocess.run(f"{plink_path}"
                             r" --noweb"
                             r" --dummy 2 2"
                             r" --freq"
                             r" --make-bed"
                             f" --out {Path(write_path, 'toy_data')}")

        # Check for process
        assert not job.returncode, "Subprocess of plink1 failed"

        # Check the two test outputs have been created, delete them if they have
        self._file_exists_cleanup(Path(write_path, 'toy_data.frq'))
        self._file_exists_cleanup(Path(write_path, 'toy_data.log'))

    def test_split_bed_by_chromosome(self):
        self.assertEqual(1, 1)

        plink_path = Path(self._config_loader()["plink_1_path"])
        load_file = Path(Path(__file__).parent, "Data", "EUR")
        write_path = self._write_path()
        assert (plink_path.exists() and write_path.exists())

        # # This pyGeneticPipe job takes the following direct args:
        # Plink_Path = apps / plink / 2.00
        # Load_File = Path

        # For each chromosome within the file specified create a new file
        for chrome in range(1, 3):
            subprocess.run(f"{plink_path}"
                           r" --noweb"
                           f" --bfile {load_file}"
                           f" --chr {chrome}"
                           r" --make-bed"
                           f" --out {load_file}_{chrome}")

        # Check the two test outputs have been created, delete them if they have
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_1.bed'))
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_1.bim'))
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_1.fam'))
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_1.log'))
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_2.bed'))
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_2.bim'))
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_2.fam'))
        self._file_exists_cleanup(Path(Path(__file__).parent, "Data", 'EUR_2.log'))


if __name__ == '__main__':
    # Path to write data to
    unittest.main()

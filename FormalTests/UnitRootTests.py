from pyGenicPipeline import ArgMaker
from pyGenicPipeline import Main

from miscSupports import load_yaml
from csvObject import CsvObject
from pathlib import Path
import numpy as np
import subprocess
import unittest
import shutil


class MyTestCase(unittest.TestCase):

    @staticmethod
    def _write_path():
        """Set the write path"""
        return Path(Path(__file__).parent, "TestDirectory")

    @staticmethod
    def _data_directory():
        """Path to the data directory"""
        return Path(Path(__file__).parent, "Data")

    @staticmethod
    def _config_loader():
        """Some processes like plink need their paths/attributes set externally from a yaml config file"""
        return load_yaml(Path(Path(__file__).parent.absolute(), "config.yaml"))

    @staticmethod
    def _file_exists_cleanup(file_path):
        """Some tests we simply want to validate a file exists, what is within them """

        assert file_path.exists(), f"Process failed to produce {file_path.name}"
        print(f"Found {file_path.name}: Deleting")
        file_path.unlink()

    @staticmethod
    def _directory_exists_cleanup(directory_path):
        """Unlink directories of files"""

        assert directory_path.exists(), f"{directory_path} was not created"
        shutil.rmtree(directory_path)

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
        load_file = Path(self._data_directory(), "EUR")
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

    def test_pgs_summary_cleaner(self):
        """
        Validates that the summary clean works
        """
        # Load the args into a yaml, then run
        self._set_summary_cleaner()
        Main(f"{self._write_path()}/pgs_clean_summary_stats.yaml")

        # Check the args was set and the output written
        self._file_exists_cleanup(Path(self._write_path(), "pgs_clean_summary_stats.yaml"))

        # Assert the output file was written
        output = Path(self._write_path(), "PGS", "Cleaned", "Cleaned_2.csv")
        assert output.exists(), "Failed to created the cleaned output"

        # Load the output, validate some of the contents to check if the file has changed between versions
        cleaned = CsvObject(output, set_columns=True, column_types=[int, int, str, str, str, float, float, float])
        self.assertEqual(cleaned.row_length, 8)
        self.assertEqual(cleaned.column_length, 39570)
        self.assertEqual(np.mean(cleaned["Beta"]), 2.804870250211524e-06)

        # Remove directory structures
        self._directory_exists_cleanup(Path(self._write_path(), "PGS"))
        self._directory_exists_cleanup(Path(self._write_path(), "PySnpTools_Meta"))

    def _set_summary_cleaner(self):
        """Use arg maker to make the yaml args"""
        args_maker = ArgMaker(self._config_loader()["defaults"])
        args_dict = args_maker.pgs_clean_summary_stats

        args_dict["Mandatory"]["Working_Directory"] = self._write_path()
        args_dict["Mandatory"]["Load_Directory"] = self._data_directory()
        args_dict["Mandatory"]["Summary_Path"] = Path(self._data_directory(), "Height.QC.gz")
        args_dict["Mandatory"]["Summary_Sample_Size"] = 253288
        args_dict["Mandatory"]["Summary_Effect_Type"] = "OR"
        args_dict["Mandatory"]["HapMap3"] = False
        args_dict["Mandatory"]["Z_Scores"] = False

        args_maker.write_yaml_args(args_dict, self._write_path())

if __name__ == '__main__':
    # Path to write data to
    unittest.main()

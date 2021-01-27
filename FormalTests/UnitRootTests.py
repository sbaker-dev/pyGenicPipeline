from pyGenicPipeline.Processors import *

from miscSupports import load_yaml
from pathlib import Path
import subprocess
import unittest
import sys
import os


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

        assert file_path.exists(), "plink failed to produce toy_data.frq"
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


if __name__ == '__main__':
    # Path to write data to
    unittest.main()

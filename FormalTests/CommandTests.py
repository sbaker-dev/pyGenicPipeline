import unittest
import subprocess
from pathlib import Path
from pyGeneticPipe.clean.Cleaner import Cleaner
import h5py


class MyTestCase(unittest.TestCase):

    @staticmethod
    def test_plink():
        """
        Test Plink 1.9 installed and can run
        """
        testing_path = Path(__file__).parent.absolute()
        plink_path = Path(r"C:\Users\Samuel\plink.exe")
        assert plink_path.exists()

        job = subprocess.run(f"{plink_path}"
                             r" --dummy 2 2"
                             r" --freq"
                             r" --make-bed"
                             f" --out {Path(testing_path, 'TestDirectory', 'toy_data')}")

        # Check for process
        assert not job.returncode, "Subprocess of plink1 failed"

        # Create paths to these files
        freq_file = Path(testing_path, "TestDirectory", 'toy_data.frq')
        log_file = Path(testing_path, "TestDirectory", 'toy_data.log')

        # Check the two test outputs have been created, delete them if they have
        assert freq_file.exists(), "plink failed to produce toy_data.frq"
        print(f"Found {freq_file.name}: Deleting")
        freq_file.unlink()

        assert log_file.exists(), "plink failed to produce toy_data.log"
        print(f"Found {log_file.name}: Deleting")
        log_file.unlink()

    def test_cvg_bed(self):
        """Test the create_validation_group for .bed"""
        testing_path = Path(__file__).parent.absolute()

        # For .bed (within config)
        bed_handel = Cleaner(Path(testing_path, "create_validation_group.yaml"))
        bed_handel.load_type = ".bed"
        print("Validating for .bed")
        self._validate_cvg(bed_handel)

    def test_cvg_bgen(self):
        """Test the create_validation_group for .bgen"""
        testing_path = Path(__file__).parent.absolute()

        # For .bgen (within config)
        bgen_handel = Cleaner(Path(testing_path, "create_validation_group.yaml"))
        print("Validating for .bgen")
        self._validate_cvg(bgen_handel)

    def _validate_cvg(self, job_handel):
        """
        We found an error that .bgen files load [rs123, rs123] yet we where expect rs123. When this happended we ended
        up with length issues as well as the problems associate with the snps. This tests for both of these conditions
        and is generalised to test on .bgen and .bed
        """
        job_handel.create_validation_group()

        with h5py.File(Path(job_handel.working_dir, job_handel.project_name), "r") as f:

            self.assertEqual(f["Validation"]["Valid_Snps"].shape[0], 540528)

            # Check we found all chromosomes of sample file
            self.assertEqual(f["Validation"]["Valid_Chromosomes"].shape[0], 23)

            # Check that we only have the snp name once
            self.assertEqual(len(f["Validation"]["Valid_Snps"][0].decode().split(",")), 1)

        # Check the two test outputs have been created, delete them if they have
        created_file = Path(job_handel.working_dir, job_handel.project_name)
        assert created_file.exists(), "Failed to write file"
        print(f"Found {created_file.name}: Deleting")
        created_file.unlink()


if __name__ == '__main__':
    unittest.main()

import unittest
import subprocess
from pathlib import Path


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


if __name__ == '__main__':
    unittest.main()

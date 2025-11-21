"""
Integration test to parse all structures from a directory of .dbn files.

This test reads all .dbn files from a specified directory, parses each one,
and records which ones fail and why. Uses multiprocessing for parallel processing.
"""

import pytest
import multiprocessing as mp
from pathlib import Path
from typing import List, Tuple, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import traceback

from rna_secstruct.parser import Parser


def read_dbn_file(filepath: Path) -> Tuple[str, str]:
    """
    Read sequence and structure from a .dbn file.
    
    Args:
        filepath: Path to the .dbn file
        
    Returns:
        Tuple of (sequence, structure)
        
    Raises:
        ValueError: If file format is invalid
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Filter out comment lines (starting with #)
    non_comment_lines = [line.strip() for line in lines if line.strip() and not line.strip().startswith('#')]
    
    if len(non_comment_lines) < 2:
        raise ValueError(f"File {filepath} does not contain both sequence and structure lines")
    
    sequence = non_comment_lines[0]
    structure = non_comment_lines[1]
    
    return sequence, structure


def parse_single_file(filepath_str: str) -> Tuple[str, Optional[str], Optional[str]]:
    """
    Parse a single .dbn file and return the result.
    
    This function must be picklable for multiprocessing, so it takes a string path.
    
    Args:
        filepath_str: String path to the .dbn file
        
    Returns:
        Tuple of (filepath_str, error_message, traceback_string)
        If successful, returns (filepath_str, None, None)
        If failed, returns (filepath_str, error_message, traceback_string)
    """
    filepath = Path(filepath_str)
    try:
        sequence, structure = read_dbn_file(filepath)
        parser = Parser()
        parser.parse(sequence, structure)
        return (str(filepath), None, None)
    except Exception as e:
        error_msg = f"{type(e).__name__}: {str(e)}"
        tb_str = traceback.format_exc()
        return (str(filepath), error_msg, tb_str)


@pytest.mark.integration
def test_parse_all_dbn_files():
    """
    Test parsing all .dbn files from the specified directory.
    
    This test:
    1. Finds all .dbn files in the directory
    2. Parses each file in parallel using multiprocessing
    3. Records which files fail and why
    4. Reports statistics at the end
    """
    dbn_dir = Path("/Users/jyesselman2/Downloads/dbnFiles")
    
    if not dbn_dir.exists():
        pytest.skip(f"Directory {dbn_dir} does not exist")
    
    # Find all .dbn files
    dbn_files = list(dbn_dir.glob("*.dbn"))
    
    if not dbn_files:
        pytest.skip(f"No .dbn files found in {dbn_dir}")
    
    print(f"\nFound {len(dbn_files)} .dbn files to parse")
    
    # Parse files in parallel using multiprocessing
    # Using multiprocessing for CPU-bound parsing work
    # Convert Path objects to strings for pickling
    n_jobs = mp.cpu_count()  # Use all available CPU cores
    
    failures: List[Tuple[str, str, str]] = []
    successes = 0
    
    # Convert Path objects to strings for multiprocessing
    dbn_file_strings = [str(f) for f in dbn_files]
    
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # Submit all tasks
        future_to_file = {
            executor.submit(parse_single_file, filepath_str): filepath_str 
            for filepath_str in dbn_file_strings
        }
        
        # Process results as they complete
        for future in as_completed(future_to_file):
            filepath_str, error_msg, tb_str = future.result()
            
            if error_msg is None:
                successes += 1
                if successes % 1000 == 0:
                    print(f"Processed {successes} files successfully...")
            else:
                failures.append((filepath_str, error_msg, tb_str))
                if len(failures) <= 20:  # Print first 20 failures
                    filepath = Path(filepath_str)
                    print(f"FAILED: {filepath.name} - {error_msg}")
    
    # Print summary
    total = len(dbn_files)
    print(f"\n{'='*60}")
    print(f"Parsing Summary:")
    print(f"  Total files: {total}")
    print(f"  Successful: {successes}")
    print(f"  Failed: {len(failures)}")
    print(f"  Success rate: {successes/total*100:.2f}%")
    print(f"{'='*60}")
    
    # Group failures by error type
    if failures:
        error_types: Dict[str, int] = {}
        for _, error_msg, _ in failures:
            error_type = error_msg.split(':')[0]  # Get exception type
            error_types[error_type] = error_types.get(error_type, 0) + 1
        
        print(f"\nFailure breakdown by error type:")
        for error_type, count in sorted(error_types.items(), key=lambda x: -x[1]):
            print(f"  {error_type}: {count}")
        
        # Save detailed failure report
        report_path = Path(__file__).parent / "parse_failures_report.txt"
        with open(report_path, 'w') as f:
            f.write("Detailed Failure Report\n")
            f.write("=" * 60 + "\n\n")
            for filepath_str, error_msg, tb_str in failures:
                filepath = Path(filepath_str)
                f.write(f"File: {filepath.name}\n")
                f.write(f"Path: {filepath_str}\n")
                f.write(f"Error: {error_msg}\n")
                if tb_str:
                    f.write(f"Traceback:\n{tb_str}\n")
                f.write("-" * 60 + "\n\n")
        
        print(f"\nDetailed failure report saved to: {report_path}")
    
    # Assert that we processed all files
    assert successes + len(failures) == total, \
        f"Expected to process {total} files, but got {successes} successes and {len(failures)} failures"
    
    # Optionally fail the test if there are too many failures
    # Uncomment the line below if you want the test to fail when failures exceed a threshold
    # failure_rate = len(failures) / total
    # assert failure_rate < 0.1, f"Failure rate {failure_rate*100:.2f}% exceeds 10% threshold"


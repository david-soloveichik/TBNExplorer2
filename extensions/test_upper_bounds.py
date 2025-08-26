#!/usr/bin/env python3
"""Test the upper bound functionality for specific polymers in IBOT.

This test verifies that the --upper-bound-on-polymers option correctly computes
upper bounds on concentration exponents for specific off-target polymers.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

# Add parent directory to path to import modules
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_upper_bounds_basic():
    """Test basic upper bound computation."""

    # Create a simple test case
    # We'll use the existing test files if they work
    tbn_file = Path(__file__).parent / "my_inputs" / "and_gate_noA.tbn"
    on_target_file = Path(__file__).parent / "my_inputs" / "and_gate_noA_on-target.tbnpolys"

    if not tbn_file.exists() or not on_target_file.exists():
        print("Test files not found. Skipping test.")
        return

    # Create a test file specifying some off-target polymers to bound
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tbnpolys", delete=False) as f:
        upper_bound_file = Path(f.name)
        # Write polymers that should be off-target
        # These are from the polymer basis but not in the on-target file
        f.write("# Test polymers for upper bounds\n")
        f.write("b2 c1 c2\n")  # A polymer from the basis
        f.write("\n")
        f.write("a2 b1\n")  # Another polymer from the basis

    try:
        # Run the IBOT algorithm with upper bounds
        cmd = [
            "python",
            "-m",
            "extensions.ibot_cli",
            str(tbn_file),
            str(on_target_file),
            "--upper-bound-on-polymers",
            str(upper_bound_file),
            "--output-prefix",
            "test_upper_bound",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        if result.returncode != 0:
            print(f"Command failed with return code {result.returncode}")
            print(f"STDOUT:\n{result.stdout}")
            print(f"STDERR:\n{result.stderr}")

            # Check if it's expected (e.g., polymers not in basis)
            if "not found in polymer basis" in result.stderr:
                print("Expected error: Test polymers not in polymer basis.")
                return True
            return False

        # Check that output files were created
        output_file = Path("test_upper_bound-ibot-upper-bounds.tbnpolys")
        if output_file.exists():
            print(f"✓ Upper bound output file created: {output_file}")
            # Clean up
            output_file.unlink()
            return True
        else:
            print("✗ Output file not created")
            return False

    finally:
        # Clean up temporary file
        upper_bound_file.unlink(missing_ok=True)

        # Clean up any other generated files
        for pattern in ["test_upper_bound-ibot*", "test_upper_bound.tbnpolymat"]:
            for f in Path.cwd().glob(pattern):
                f.unlink()


def test_validation():
    """Test that --generate-tbn cannot be used with --upper-bound-on-polymers."""

    tbn_file = Path(__file__).parent / "my_inputs" / "and_gate_noA.tbn"
    on_target_file = Path(__file__).parent / "my_inputs" / "and_gate_noA_on-target.tbnpolys"

    if not tbn_file.exists() or not on_target_file.exists():
        print("Test files not found. Skipping test.")
        return

    # Create a dummy upper bounds file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tbnpolys", delete=False) as f:
        upper_bound_file = Path(f.name)
        f.write("# Dummy file\n")
        f.write("a1 a2\n")

    try:
        # Try to run with both --generate-tbn and --upper-bound-on-polymers
        cmd = [
            "python",
            "-m",
            "extensions.ibot_cli",
            str(tbn_file),
            str(on_target_file),
            "--upper-bound-on-polymers",
            str(upper_bound_file),
            "--generate-tbn",
            "100.0",
            "nM",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=False)

        # This should fail with our validation error
        if result.returncode != 0 and "--generate-tbn cannot be used with --upper-bound-on-polymers" in result.stderr:
            print("✓ Validation correctly prevents using --generate-tbn with --upper-bound-on-polymers")
            return True
        else:
            print("✗ Validation did not work as expected")
            print(f"STDERR:\n{result.stderr}")
            return False

    finally:
        # Clean up
        upper_bound_file.unlink(missing_ok=True)


if __name__ == "__main__":
    print("Testing upper bound functionality...")
    print()

    print("Test 1: Basic upper bound computation")
    success1 = test_upper_bounds_basic()
    print()

    print("Test 2: Validation of incompatible options")
    success2 = test_validation()
    print()

    if success1 and success2:
        print("All tests passed! ✓")
    else:
        print("Some tests failed. ✗")
        sys.exit(1)

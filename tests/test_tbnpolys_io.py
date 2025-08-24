"""Tests for the tbnpolys_io module."""

import tempfile
from pathlib import Path

import pytest

from tbnexplorer2.model import TBN, BindingSite, Monomer
from tbnexplorer2.tbnpolys_io import TbnpolysParser, TbnpolysWriter


def create_test_tbn():
    """Create a simple TBN for testing."""
    binding_sites = {"a": 0, "b": 1, "c": 2}

    monomers = [
        Monomer(
            binding_sites=[
                BindingSite("a", False),
                BindingSite("a", True),
                BindingSite("b", False),
            ],
            concentration=100.0,
            name="monomer1",
            original_line="monomer1: a a* b, 100.0",
        ),
        Monomer(
            binding_sites=[
                BindingSite("b", False),
                BindingSite("c", True),
            ],
            concentration=50.0,
            name=None,  # Unnamed monomer
            original_line="b c*, 50.0",
        ),
        Monomer(
            binding_sites=[
                BindingSite("c", False),
                BindingSite("c", False),
            ],
            concentration=75.0,
            name="C",
            original_line="c c >C, 75.0",
        ),
    ]

    return TBN(monomers, binding_sites, concentration_units="nM")


class TestTbnpolysWriter:
    """Test TbnpolysWriter class."""

    def test_write_simple_polymer(self):
        """Test writing a simple polymer."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)

        # Create a polymer with counts [2, 1, 0]
        polymers = [[2, 1, 0]]

        content = writer.format_polymers(polymers)

        assert "2 | monomer1" in content
        assert "b c*" in content  # Unnamed monomer shows binding sites
        assert "C" not in content  # Monomer with count 0 not shown

    def test_write_multiple_polymers(self):
        """Test writing multiple polymers."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)

        polymers = [
            [1, 0, 1],  # monomer1 and C
            [0, 2, 0],  # Only unnamed monomer
        ]

        content = writer.format_polymers(polymers)
        lines = content.split("\n")

        # First polymer
        assert "monomer1" in content
        assert "C" in content

        # Second polymer
        assert "2 | b c*" in content

        # Check for empty line between polymers
        empty_line_count = sum(1 for line in lines if line.strip() == "")
        assert empty_line_count >= 1

    def test_write_with_concentrations(self):
        """Test writing polymers with concentrations."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)

        polymers = [[1, 1, 0]]
        concentrations = [25.5]

        content = writer.format_polymers(polymers, concentrations, units="nM")

        assert "# Concentration: 25.50 nM" in content

    def test_write_with_header_comment(self):
        """Test writing with a header comment."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)

        polymers = [[1, 0, 0]]
        header = "Test polymer basis\nWith multiple lines"

        content = writer.format_polymers(polymers, header_comment=header)

        assert "# Test polymer basis" in content
        assert "# With multiple lines" in content

    def test_format_concentration_values(self):
        """Test concentration formatting for various values."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)

        # Test various concentration values
        test_cases = [
            (0, "0"),
            (0.00001, "1.00e-05"),
            (0.005, "5.00e-03"),  # < 0.01 uses scientific notation
            (0.01, "0.0100"),  # >= 0.01 uses 4 decimal places
            (0.5, "0.5000"),  # >= 0.01 uses 4 decimal places
            (1.5, "1.500"),  # >= 1 uses 3 decimal places
            (5.5, "5.500"),  # >= 1 uses 3 decimal places
            (15.5, "15.50"),  # >= 10 uses 2 decimal places
            (55.5, "55.50"),  # >= 10 uses 2 decimal places
            (155.5, "155.5"),  # >= 100 uses 1 decimal place
            (555.5, "555.5"),  # >= 100 uses 1 decimal place
            (5555.5, "5.6e+03"),  # >= 1000 uses scientific notation
        ]

        for value, expected_prefix in test_cases:
            formatted = writer._format_concentration(value, "nM")
            assert formatted.startswith(expected_prefix), f"Failed for {value}: got {formatted}"

    def test_write_to_file(self):
        """Test writing polymers to a file."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)

        polymers = [[1, 1, 1]]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbnpolys", delete=False) as f:
            file_path = Path(f.name)

        try:
            writer.write_polymers(polymers, file_path)

            # Read back and verify
            with open(file_path) as f:
                content = f.read()

            assert "monomer1" in content
            assert "b c*" in content
            assert "C" in content
        finally:
            file_path.unlink()


class TestTbnpolysParser:
    """Test TbnpolysParser class."""

    def test_parse_simple_polymer(self):
        """Test parsing a simple polymer."""
        content = """# Comment line
monomer1
b c*
"""
        parser = TbnpolysParser()
        polymers = parser.parse_content(content)

        assert len(polymers) == 1
        assert len(polymers[0]) == 2
        assert polymers[0][0] == (1, "monomer1")
        assert polymers[0][1] == (1, "b c*")

    def test_parse_polymer_with_multiplicity(self):
        """Test parsing polymers with multiplicity prefixes."""
        content = """2 | monomer1
3 | a b* c
monomer2
"""
        parser = TbnpolysParser()
        polymers = parser.parse_content(content)

        assert len(polymers) == 1
        assert polymers[0][0] == (2, "monomer1")
        assert polymers[0][1] == (3, "a b* c")
        assert polymers[0][2] == (1, "monomer2")

    def test_parse_multiple_polymers(self):
        """Test parsing multiple polymers separated by empty lines."""
        content = """# Polymer 1
monomer1
2 | monomer2

# Polymer 2
3 | a b*
c c*
"""
        parser = TbnpolysParser()
        polymers = parser.parse_content(content)

        assert len(polymers) == 2

        # First polymer
        assert len(polymers[0]) == 2
        assert polymers[0][0] == (1, "monomer1")
        assert polymers[0][1] == (2, "monomer2")

        # Second polymer
        assert len(polymers[1]) == 2
        assert polymers[1][0] == (3, "a b*")
        assert polymers[1][1] == (1, "c c*")

    def test_parse_with_comments(self):
        """Test parsing with various comment styles."""
        content = """# Header comment
monomer1  # inline comment
# Comment between monomers
2 | monomer2

# Next polymer
monomer3
"""
        parser = TbnpolysParser()
        polymers = parser.parse_content(content)

        assert len(polymers) == 2
        assert polymers[0][0] == (1, "monomer1")
        assert polymers[0][1] == (2, "monomer2")
        assert polymers[1][0] == (1, "monomer3")

    def test_parse_empty_content(self):
        """Test parsing empty content."""
        parser = TbnpolysParser()

        # Empty string
        polymers = parser.parse_content("")
        assert len(polymers) == 0

        # Only comments
        polymers = parser.parse_content("# Just comments\n# More comments")
        assert len(polymers) == 0

        # Only empty lines
        polymers = parser.parse_content("\n\n\n")
        assert len(polymers) == 0

    def test_parse_with_tbn_context(self):
        """Test parsing with TBN context to resolve monomers."""
        tbn = create_test_tbn()
        parser = TbnpolysParser(tbn)

        content = """monomer1
b c*
C
"""
        polymers = parser.parse_content(content)

        assert len(polymers) == 1
        assert len(polymers[0]) == 3

        # With TBN context, monomers are resolved
        assert polymers[0][0][0] == 1  # multiplicity
        assert polymers[0][0][1] == tbn.monomers[0]  # monomer1
        assert polymers[0][1][1] == tbn.monomers[1]  # unnamed monomer
        assert polymers[0][2][1] == tbn.monomers[2]  # C

    def test_parse_binding_sites_different_order(self):
        """Test that binding sites can be in different order."""
        tbn = create_test_tbn()
        parser = TbnpolysParser(tbn)

        # The unnamed monomer has binding sites "b c*"
        # Test with reversed order
        content = """c* b"""

        polymers = parser.parse_content(content)
        assert len(polymers) == 1
        assert polymers[0][0][1] == tbn.monomers[1]

    def test_parse_file(self):
        """Test parsing from a file."""
        content = """monomer1
2 | monomer2
"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbnpolys", delete=False) as f:
            f.write(content)
            file_path = Path(f.name)

        try:
            parser = TbnpolysParser()
            polymers = parser.parse_file(file_path)

            assert len(polymers) == 1
            assert polymers[0][0] == (1, "monomer1")
            assert polymers[0][1] == (2, "monomer2")
        finally:
            file_path.unlink()

    def test_invalid_monomer_resolution(self):
        """Test error handling for invalid monomer resolution."""
        tbn = create_test_tbn()
        parser = TbnpolysParser(tbn)

        # Try to resolve a non-existent monomer
        content = """nonexistent_monomer"""

        with pytest.raises(ValueError, match="Cannot resolve monomer"):
            parser.parse_content(content)


class TestRoundTrip:
    """Test round-trip writing and parsing."""

    def test_round_trip_simple(self):
        """Test that we can write and read back the same data."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)
        parser = TbnpolysParser(tbn)

        # Original polymers
        original_polymers = [
            [1, 0, 1],
            [2, 1, 0],
            [0, 3, 0],
        ]

        # Write to string
        content = writer.format_polymers(original_polymers)

        # Parse back
        parsed_polymers = parser.parse_content(content)

        # Verify we get the same polymers back
        assert len(parsed_polymers) == len(original_polymers)

        for orig_poly, parsed_poly in zip(original_polymers, parsed_polymers):
            # Reconstruct the polymer vector from parsed data
            reconstructed = [0] * len(tbn.monomers)
            for multiplicity, monomer in parsed_poly:
                # Find monomer index
                monomer_idx = tbn.monomers.index(monomer)
                reconstructed[monomer_idx] = multiplicity

            assert reconstructed == orig_poly

    def test_round_trip_with_concentrations(self):
        """Test round-trip with concentrations."""
        tbn = create_test_tbn()
        writer = TbnpolysWriter(tbn)

        polymers = [[1, 1, 0]]
        concentrations = [42.5]

        # Write with concentrations
        content = writer.format_polymers(polymers, concentrations, units="nM")

        # The parser doesn't parse concentrations (they're just comments)
        # but we can verify the format is preserved
        assert "# Concentration: 42.50 nM" in content

import os
import tempfile

import pytest

from tbnexplorer2.parser import TBNParser
from tbnexplorer2.polymat_io import PolymatData, PolymatReader, PolymatWriter


class TestParametrizedTBN:
    """Test parametrized TBN file parsing with {{var}} template syntax."""

    def test_simple_parametrized_concentration(self):
        """Test parsing a simple parametrized concentration."""
        content = """\\UNITS: nM
monomer1: a b, {{conc1}}
monomer2: c d, {{conc2}}"""

        variables = {"conc1": 100.0, "conc2": 50.5}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            monomers, binding_sites, units, used_vars = TBNParser.parse_file(f.name, variables=variables)

        os.unlink(f.name)

        assert len(monomers) == 2
        assert monomers[0].concentration == 100.0
        assert monomers[1].concentration == 50.5
        assert used_vars == {"conc1": 100.0, "conc2": 50.5}

    def test_missing_variable_error(self):
        """Test that missing template variables raise an error."""
        content = """\\UNITS: nM
monomer1: a b, {{missing_var}}"""

        variables = {"other_var": 100.0}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            with pytest.raises(ValueError, match="Template variable 'missing_var' not provided"):
                TBNParser.parse_file(f.name, variables=variables)

        os.unlink(f.name)

    def test_mixed_parametrized_and_literal(self):
        """Test mixing parametrized and literal concentrations."""
        content = """\\UNITS: nM
monomer1: a b, {{var1}}
monomer2: c d, 75.5
monomer3: e f, {{var2}}"""

        variables = {"var1": 100.0, "var2": 25.0}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            monomers, binding_sites, units, used_vars = TBNParser.parse_file(f.name, variables=variables)

        os.unlink(f.name)

        assert len(monomers) == 3
        assert monomers[0].concentration == 100.0
        assert monomers[1].concentration == 75.5
        assert monomers[2].concentration == 25.0
        assert used_vars == {"var1": 100.0, "var2": 25.0}

    def test_unused_variables(self):
        """Test that unused variables don't appear in used_vars."""
        content = """\\UNITS: nM
monomer1: a b, {{var1}}"""

        variables = {"var1": 100.0, "var2": 50.0, "unused": 25.0}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            monomers, binding_sites, units, used_vars = TBNParser.parse_file(f.name, variables=variables)

        os.unlink(f.name)

        assert monomers[0].concentration == 100.0
        assert used_vars == {"var1": 100.0}
        assert "var2" not in used_vars
        assert "unused" not in used_vars

    def test_template_with_spaces(self):
        """Test that template syntax works with spaces inside braces."""
        content = """\\UNITS: nM
monomer1: a b, {{ var1 }}
monomer2: c d, {{var2}}"""

        variables = {"var1": 100.0, "var2": 50.0}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            monomers, binding_sites, units, used_vars = TBNParser.parse_file(f.name, variables=variables)

        os.unlink(f.name)

        assert monomers[0].concentration == 100.0
        assert monomers[1].concentration == 50.0

    def test_parameters_in_polymat_file(self):
        """Test that parameters are saved and loaded from .tbnpolymat files."""
        # Create a polymat data with parameters
        import numpy as np

        data = PolymatData(
            polymers=[np.array([1, 0, 0]), np.array([0, 1, 0])],
            n_monomers=3,
            n_polymers=2,
            parameters={"conc1": 100.0, "conc2": 50.5},
            matrix_hash="test_hash",
        )

        with tempfile.NamedTemporaryFile(suffix=".tbnpolymat", delete=False) as f:
            writer = PolymatWriter(f.name)
            writer.write(data)

            # Read it back
            reader = PolymatReader(f.name)
            loaded_data = reader.read()

        os.unlink(f.name)

        assert loaded_data.parameters == {"conc1": 100.0, "conc2": 50.5}
        assert loaded_data.matrix_hash == "test_hash"

    def test_no_parameters_in_polymat_file(self):
        """Test that polymat files without parameters work correctly."""
        import numpy as np

        data = PolymatData(
            polymers=[np.array([1, 0, 0])],
            n_monomers=3,
            n_polymers=1,
            matrix_hash="test_hash",
        )

        with tempfile.NamedTemporaryFile(suffix=".tbnpolymat", delete=False) as f:
            writer = PolymatWriter(f.name)
            writer.write(data)

            # Read it back
            reader = PolymatReader(f.name)
            loaded_data = reader.read()

        os.unlink(f.name)

        assert loaded_data.parameters is None
        assert loaded_data.matrix_hash == "test_hash"

    def test_negative_concentration_with_units(self):
        """Test that negative concentrations are allowed with parametrized values when units are present."""
        content = """\\UNITS: nM
monomer1: a b, {{positive}}
monomer2: c d, {{negative}}"""

        variables = {"positive": 100.0, "negative": -50.0}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            monomers, binding_sites, units, used_vars = TBNParser.parse_file(f.name, variables=variables)

        os.unlink(f.name)

        # This should work fine since units are present
        assert len(monomers) == 2
        assert monomers[0].concentration == 100.0
        assert monomers[1].concentration == -50.0

    def test_invalid_template_syntax(self):
        """Test that invalid template syntax is treated as literal."""
        content = """\\UNITS: nM
monomer1: a b, 100.0"""

        variables = {}

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            monomers, binding_sites, units, used_vars = TBNParser.parse_file(f.name, variables=variables)

        os.unlink(f.name)

        assert monomers[0].concentration == 100.0
        assert len(used_vars) == 0

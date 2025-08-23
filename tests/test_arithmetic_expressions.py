"""Test arithmetic expressions in parametrized TBN files."""

import tempfile
from pathlib import Path

import pytest

from tbnexplorer2.parser import TBNParser


class TestArithmeticExpressions:
    """Test arithmetic expression evaluation in template syntax."""

    def test_simple_addition(self):
        """Test simple addition in template."""
        content = """\\UNITS: nM
monomer1: a b, {{x + y}}
monomer2: c d, {{a + 10}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 50, "y": 30, "a": 90}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert len(monomers) == 2
            assert monomers[0].concentration == 80  # 50 + 30
            assert monomers[1].concentration == 100  # 90 + 10
            assert "x" in used_vars
            assert "y" in used_vars
            assert "a" in used_vars

            Path(f.name).unlink()

    def test_subtraction(self):
        """Test subtraction in template."""
        content = """\\UNITS: nM
monomer1: a b, {{base - offset}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"base": 100, "offset": 25}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 75  # 100 - 25
            assert "base" in used_vars
            assert "offset" in used_vars

            Path(f.name).unlink()

    def test_multiplication(self):
        """Test multiplication in template."""
        content = """\\UNITS: nM
monomer1: a b, {{conc * factor}}
monomer2: c d, {{base * 2.5}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"conc": 20, "factor": 3, "base": 40}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 60  # 20 * 3
            assert monomers[1].concentration == 100  # 40 * 2.5

            Path(f.name).unlink()

    def test_division(self):
        """Test division in template."""
        content = """\\UNITS: nM
monomer1: a b, {{total / parts}}
monomer2: c d, {{value / 2}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"total": 150, "parts": 3, "value": 90}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 50  # 150 / 3
            assert monomers[1].concentration == 45  # 90 / 2

            Path(f.name).unlink()

    def test_power_operation(self):
        """Test power/exponentiation in template."""
        content = """\\UNITS: nM
monomer1: a b, {{base ** exp}}
monomer2: c d, {{2 ** n}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"base": 5, "exp": 2, "n": 3}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 25  # 5 ** 2
            assert monomers[1].concentration == 8  # 2 ** 3

            Path(f.name).unlink()

    def test_parentheses_and_order(self):
        """Test parentheses and order of operations."""
        content = """\\UNITS: nM
monomer1: a b, {{(a + b) * c}}
monomer2: c d, {{a + b * c}}
monomer3: e f, {{((x + y) * 2) / z}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"a": 10, "b": 20, "c": 3, "x": 15, "y": 5, "z": 4}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 90  # (10 + 20) * 3
            assert monomers[1].concentration == 70  # 10 + (20 * 3)
            assert monomers[2].concentration == 10  # ((15 + 5) * 2) / 4

            Path(f.name).unlink()

    def test_floating_point_arithmetic(self):
        """Test floating point arithmetic."""
        content = """\\UNITS: nM
monomer1: a b, {{x * 1.5 + y * 0.5}}
monomer2: c d, {{(a + b) / 2.0}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 40, "y": 30, "a": 75, "b": 45}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 75  # 40 * 1.5 + 30 * 0.5 = 60 + 15
            assert monomers[1].concentration == 60  # (75 + 45) / 2.0 = 120 / 2

            Path(f.name).unlink()

    def test_mixed_templates_and_literals(self):
        """Test mixing template expressions with literal values."""
        content = """\\UNITS: nM
monomer1: a b, {{base * 2}}
monomer2: c d, 50.5
monomer3: e f, {{base + offset}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"base": 25, "offset": 10}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 50  # 25 * 2
            assert monomers[1].concentration == 50.5  # literal
            assert monomers[2].concentration == 35  # 25 + 10

            Path(f.name).unlink()

    def test_backward_compatibility_simple_var(self):
        """Test that simple variable syntax still works."""
        content = """\\UNITS: nM
monomer1: a b, {{conc1}}
monomer2: c d, {{conc2}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"conc1": 100, "conc2": 75.5}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 100
            assert monomers[1].concentration == 75.5
            assert "conc1" in used_vars
            assert "conc2" in used_vars

            Path(f.name).unlink()

    def test_spaces_in_expression(self):
        """Test expressions with various spacing."""
        content = """\\UNITS: nM
monomer1: a b, {{ x + y }}
monomer2: c d, {{a+b}}
monomer3: e f, {{  (m * n)  +  p  }}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 30, "y": 20, "a": 10, "b": 15, "m": 5, "n": 8, "p": 10}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 50  # 30 + 20
            assert monomers[1].concentration == 25  # 10 + 15
            assert monomers[2].concentration == 50  # (5 * 8) + 10

            Path(f.name).unlink()

    def test_negative_results_with_units(self):
        """Test that negative results are allowed when UNITS is specified."""
        content = """\\UNITS: nM
monomer1: a b, {{x - y}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 20, "y": 50}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == -30  # 20 - 50

            Path(f.name).unlink()

    def test_error_missing_variable_in_expression(self):
        """Test error when variable in expression is not provided."""
        content = """\\UNITS: nM
monomer1: a b, {{x + missing_var}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 50}
            with pytest.raises(ValueError) as excinfo:
                TBNParser.parse_file(f.name, variables)

            assert "missing_var" in str(excinfo.value)
            assert "not provided" in str(excinfo.value)

            Path(f.name).unlink()

    def test_error_invalid_expression(self):
        """Test error for invalid expressions."""
        content = """\\UNITS: nM
monomer1: a b, {{x +* y}}"""  # Invalid syntax

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 50, "y": 30}
            with pytest.raises(ValueError) as excinfo:
                TBNParser.parse_file(f.name, variables)

            assert "Invalid expression" in str(excinfo.value) or "Error evaluating" in str(excinfo.value)

            Path(f.name).unlink()

    def test_error_division_by_zero(self):
        """Test error for division by zero."""
        content = """\\UNITS: nM
monomer1: a b, {{x / y}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 100, "y": 0}
            with pytest.raises(ValueError) as excinfo:
                TBNParser.parse_file(f.name, variables)

            assert "Error evaluating" in str(excinfo.value) or "Invalid expression" in str(excinfo.value)

            Path(f.name).unlink()

    def test_complex_expression(self):
        """Test a more complex real-world expression."""
        content = """\\UNITS: nM
# Base concentration with scaling factor
monomer1: a b, {{base_conc * scale_factor}}
# Average of two concentrations
monomer2: c d, {{(conc1 + conc2) / 2}}
# Dilution calculation
monomer3: e f, {{stock_conc * (dilution_vol / total_vol)}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {
                "base_conc": 100,
                "scale_factor": 1.5,
                "conc1": 80,
                "conc2": 120,
                "stock_conc": 1000,
                "dilution_vol": 5,
                "total_vol": 50,
            }
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert monomers[0].concentration == 150  # 100 * 1.5
            assert monomers[1].concentration == 100  # (80 + 120) / 2
            assert monomers[2].concentration == 100  # 1000 * (5 / 50)

            Path(f.name).unlink()

    def test_unused_variables_not_tracked(self):
        """Test that variables not used in expressions are not tracked."""
        content = """\\UNITS: nM
monomer1: a b, {{x + y}}"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tbn", delete=False) as f:
            f.write(content)
            f.flush()

            variables = {"x": 50, "y": 30, "z": 100, "unused": 200}
            monomers, _, units, used_vars = TBNParser.parse_file(f.name, variables)

            assert "x" in used_vars
            assert "y" in used_vars
            assert "z" not in used_vars
            assert "unused" not in used_vars

            Path(f.name).unlink()

import pytest
import tempfile
import os
from tbnexplorer2.parser import TBNParser
from tbnexplorer2.model import BindingSite


class TestTBNParser:
    
    def test_parse_simple_monomer(self):
        """Test parsing a simple monomer without name or concentration."""
        content = "a b c"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 1
        assert monomers[0].name is None
        assert monomers[0].concentration is None
        assert len(monomers[0].binding_sites) == 3
        assert monomers[0].binding_sites[0] == BindingSite("a", False)
        assert monomers[0].binding_sites[1] == BindingSite("b", False)
        assert monomers[0].binding_sites[2] == BindingSite("c", False)
    
    def test_parse_star_binding_sites(self):
        """Test parsing star and unstar binding sites."""
        content = "a a* b b*"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 1
        assert monomers[0].binding_sites[0] == BindingSite("a", False)
        assert monomers[0].binding_sites[1] == BindingSite("a", True)
        assert monomers[0].binding_sites[2] == BindingSite("b", False)
        assert monomers[0].binding_sites[3] == BindingSite("b", True)
    
    def test_parse_named_monomer(self):
        """Test parsing a monomer with a name."""
        content = "monomer1: a b c"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 1
        assert monomers[0].name == "monomer1"
        assert len(monomers[0].binding_sites) == 3
    
    def test_parse_monomer_with_concentration(self):
        """Test parsing a monomer with concentration."""
        content = "UNITS: nM\na b c, 100.5"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 1
        assert monomers[0].concentration == 100.5
        assert units == "nM"
    
    def test_parse_named_monomer_with_concentration(self):
        """Test parsing a named monomer with concentration."""
        content = "UNITS: nM\nmol1: a b c, 50.7"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 1
        assert monomers[0].name == "mol1"
        assert monomers[0].concentration == 50.7
    
    def test_parse_multiple_monomers(self):
        """Test parsing multiple monomers."""
        content = """
        monomer1: a a* b2 b1
        b2* b2* b1 b1
        a b2* b1
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 3
        assert monomers[0].name == "monomer1"
        assert monomers[0].concentration is None
        assert monomers[1].name is None
        assert monomers[1].concentration is None
        assert monomers[2].concentration is None
    
    def test_parse_with_comments(self):
        """Test parsing file with comments."""
        content = """
        # This is a comment
        monomer1: a b # inline comment
        # Another comment
        c d e
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 2
        assert monomers[0].name == "monomer1"
        assert len(monomers[0].binding_sites) == 2
        assert monomers[1].name is None
        assert len(monomers[1].binding_sites) == 3
    
    def test_parse_empty_lines(self):
        """Test parsing file with empty lines."""
        content = """
        
        a b c
        
        
        d e f
        
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 2
    
    def test_inconsistent_concentration_error(self):
        """Test that mixing concentration specifications raises an error."""
        content = """
        a b c, 100
        d e f
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Monomer has concentration but no UNITS specified"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_invalid_binding_site_error(self):
        """Test that invalid binding site characters raise errors."""
        content = "a b|c d"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Invalid binding site"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_negative_concentration_error(self):
        """Test that negative concentrations raise an error."""
        content = "a b c, -10"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Negative concentration"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_binding_site_index(self):
        """Test that binding site index is correctly built."""
        content = """
        a b c
        b c d
        a* d* e
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        # Should have 5 unique binding sites: a, b, c, d, e
        assert len(binding_site_index) == 5
        assert "a" in binding_site_index
        assert "b" in binding_site_index
        assert "c" in binding_site_index
        assert "d" in binding_site_index
        assert "e" in binding_site_index
    
    def test_monomer_name_with_spaces_error(self):
        """Test that monomer names with spaces raise an error."""
        content = "my monomer: a b c"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="cannot contain spaces"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_monomer_name_with_prohibited_chars_error(self):
        """Test that monomer names with prohibited characters raise an error."""
        test_cases = [
            ("mon,omer: a b c", "cannot contain ,"),
            ("mon*omer: a b c", "cannot contain \\*"),
            ("mon|omer: a b c", "cannot contain \\|"),
            ("mon:omer: a b c", "cannot contain :"),
        ]
        
        for content, expected_msg in test_cases:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
                f.write(content)
                f.flush()
                
                with pytest.raises(ValueError, match="cannot contain"):
                    TBNParser.parse_file(f.name)
                
            os.unlink(f.name)
    
    def test_monomer_binding_site_name_conflict_error(self):
        """Test that conflicting monomer and binding site names raise an error."""
        # Monomer name conflicts with later binding site
        content = """
        mysite: a b c
        d e mysite
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Binding site 'mysite' conflicts with monomer name"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        # Binding site conflicts with later monomer name
        content = """
        a b myname
        myname: d e f
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Monomer name 'myname' conflicts with binding site"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_case_sensitive_names_allowed(self):
        """Test that case differences are allowed for monomer vs binding site names."""
        content = """
        C: a b d
        a b c d
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            # Should not raise error since 'C' (monomer) != 'c' (binding site)
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert len(monomers) == 2
        assert monomers[0].name == "C"
        assert "c" in binding_sites  # lowercase c as binding site
        # C should not be in binding sites since it's only used as monomer name
        assert "a" in binding_sites
        assert "b" in binding_sites
        assert "d" in binding_sites
    
    def test_units_parsing(self):
        """Test that UNITS keyword is parsed correctly."""
        content = "UNITS: uM\na b c, 100"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert units == "uM"
        assert len(monomers) == 1
        assert monomers[0].concentration == 100
    
    def test_units_with_comments(self):
        """Test that UNITS works with comments before it."""
        content = """# This is a comment
# Another comment
UNITS: mM
# Comment after units
a b c, 50"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert units == "mM"
        assert len(monomers) == 1
    
    def test_no_units_no_concentrations(self):
        """Test file without UNITS and without concentrations."""
        content = "a b c\nd e f"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_sites, units = TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
        
        assert units is None
        assert len(monomers) == 2
        assert all(m.concentration is None for m in monomers)
    
    def test_invalid_units_error(self):
        """Test that invalid units raise an error."""
        content = "UNITS: invalid\na b c, 100"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Invalid units 'invalid'"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_units_without_concentrations_error(self):
        """Test that UNITS specified but no concentrations raises error."""
        content = "UNITS: nM\na b c"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="UNITS specified but monomer lacks concentration"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_multiple_units_error(self):
        """Test that multiple UNITS specifications raise an error."""
        content = """UNITS: nM
UNITS: uM
a b c, 100"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Multiple UNITS specifications found"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)


class TestMonomerRepetition:
    """Test cases for monomer repetition handling."""
    
    def test_identical_monomers_with_units_aggregate(self):
        """Test that identical monomers with UNITS aggregate their concentrations."""
        content = """UNITS: nM
a b c*, 100
a b c*, 50
d e*, 25"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
            # Should have 2 monomers (one aggregated, one unique)
            assert len(monomers) == 2
            assert units == "nM"
            
            # Find the aggregated monomer
            aggregated = None
            unique = None
            for m in monomers:
                if len([site.name for site in m.binding_sites if site.name in ['a', 'b', 'c']]) == 3:
                    aggregated = m
                else:
                    unique = m
            
            assert aggregated is not None
            assert unique is not None
            assert aggregated.concentration == 150.0  # 100 + 50
            assert unique.concentration == 25.0
            
        os.unlink(f.name)
    
    def test_identical_monomers_different_order_aggregate(self):
        """Test that monomers with same binding sites in different order aggregate."""
        content = """UNITS: nM
a b c*, 100
c* b a, 50"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
            # Should have 1 aggregated monomer
            assert len(monomers) == 1
            assert monomers[0].concentration == 150.0  # 100 + 50
            
        os.unlink(f.name)
    
    def test_identical_monomers_with_negative_concentrations(self):
        """Test that negative concentrations are allowed and sum correctly with UNITS."""
        content = """UNITS: nM
a b c*, 100
a b c*, -30
a b c*, 50"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
            # Should have 1 aggregated monomer
            assert len(monomers) == 1
            assert monomers[0].concentration == 120.0  # 100 + (-30) + 50
            
        os.unlink(f.name)
    
    def test_negative_final_concentration_error(self):
        """Test that negative final concentration after aggregation raises error."""
        content = """UNITS: nM
a b c*, 50
a b c*, -100"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            with pytest.raises(ValueError, match="Negative final concentration.*after aggregating"):
                TBNParser.parse_file(f.name)
            
        os.unlink(f.name)
    
    def test_identical_named_monomers_aggregate(self):
        """Test that identical named monomers aggregate."""
        content = """UNITS: nM
monomer1: a b c*, 100
monomer1: a b c*, 50"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
            # Should have 1 aggregated monomer
            assert len(monomers) == 1
            assert monomers[0].name == "monomer1"
            assert monomers[0].concentration == 150.0  # 100 + 50
            
        os.unlink(f.name)
    
    def test_different_monomers_dont_aggregate(self):
        """Test that different monomers don't aggregate."""
        content = """UNITS: nM
a b c*, 100
a b d*, 50"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
            # Should have 2 separate monomers
            assert len(monomers) == 2
            concentrations = [m.concentration for m in monomers]
            assert 100.0 in concentrations
            assert 50.0 in concentrations
            
        os.unlink(f.name)
    
    def test_identical_monomers_without_units_stay_separate(self):
        """Test that identical monomers without UNITS stay separate."""
        content = """a b c*
a b c*
d e*"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
            # Should have 3 separate monomers (no aggregation)
            assert len(monomers) == 3
            assert units is None
            
        os.unlink(f.name)
    
    def test_complex_aggregation_scenario(self):
        """Test complex scenario with multiple groups of identical monomers."""
        content = """UNITS: nM
# Group 1: identical monomers
a b*, 100
b* a, 50
# Group 2: different monomer 
c d*, 25
# Group 3: another set of identical monomers
e f* g, 30
g e f*, 20"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tbn', delete=False) as f:
            f.write(content)
            f.flush()
            
            monomers, binding_site_index, units = TBNParser.parse_file(f.name)
            
            # Should have 3 monomers (2 aggregated groups + 1 unique)
            assert len(monomers) == 3
            
            # Sort by concentration to make testing easier
            monomers.sort(key=lambda m: m.concentration)
            
            assert monomers[0].concentration == 25.0  # unique monomer
            assert monomers[1].concentration == 50.0  # aggregated from 30 + 20
            assert monomers[2].concentration == 150.0  # aggregated from 100 + 50
            
        os.unlink(f.name)
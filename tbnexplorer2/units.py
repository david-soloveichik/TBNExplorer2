"""
Concentration unit conversion utilities for TBN Explorer 2.

Supports conversion between common concentration units:
- pM (picoMolar): 10^-12 M
- nM (nanoMolar): 10^-9 M  
- uM (microMolar): 10^-6 M
- mM (milliMolar): 10^-3 M
- M (Molar): 1 M
"""

from typing import Union
import numpy as np


# Conversion factors to Molar
UNIT_TO_MOLAR = {
    'pM': 1e-12,
    'nM': 1e-9,
    'uM': 1e-6,
    'mM': 1e-3,
    'M': 1.0
}

# Valid concentration units
VALID_UNITS = list(UNIT_TO_MOLAR.keys())


def validate_unit(unit: str) -> str:
    """
    Validate that a concentration unit is supported.
    
    Args:
        unit: Unit string to validate
        
    Returns:
        The validated unit string
        
    Raises:
        ValueError: If unit is not supported
    """
    if unit not in VALID_UNITS:
        valid_units_str = ', '.join(VALID_UNITS)
        raise ValueError(f"Invalid concentration unit '{unit}'. Supported units: {valid_units_str}")
    return unit


def to_molar(value: Union[float, np.ndarray], from_unit: str) -> Union[float, np.ndarray]:
    """
    Convert concentration value(s) from given unit to Molar.
    
    Args:
        value: Concentration value(s) to convert
        from_unit: Source unit (pM, nM, uM, mM, M)
        
    Returns:
        Concentration value(s) in Molar
        
    Raises:
        ValueError: If from_unit is not supported
    """
    validate_unit(from_unit)
    factor = UNIT_TO_MOLAR[from_unit]
    return value * factor


def from_molar(value: Union[float, np.ndarray], to_unit: str) -> Union[float, np.ndarray]:
    """
    Convert concentration value(s) from Molar to given unit.
    
    Args:
        value: Concentration value(s) in Molar to convert
        to_unit: Target unit (pM, nM, uM, mM, M)
        
    Returns:
        Concentration value(s) in target unit
        
    Raises:
        ValueError: If to_unit is not supported
    """
    validate_unit(to_unit)
    factor = UNIT_TO_MOLAR[to_unit]
    return value / factor


def convert_concentration(value: Union[float, np.ndarray], from_unit: str, to_unit: str) -> Union[float, np.ndarray]:
    """
    Convert concentration value(s) between units.
    
    Args:
        value: Concentration value(s) to convert
        from_unit: Source unit (pM, nM, uM, mM, M)  
        to_unit: Target unit (pM, nM, uM, mM, M)
        
    Returns:
        Concentration value(s) in target unit
        
    Raises:
        ValueError: If either unit is not supported
    """
    if from_unit == to_unit:
        return value
    
    # Convert through Molar
    molar_value = to_molar(value, from_unit)
    return from_molar(molar_value, to_unit)


def get_unit_display_name(unit: str) -> str:
    """
    Get the full display name for a concentration unit.
    
    Args:
        unit: Unit abbreviation (pM, nM, uM, mM, M)
        
    Returns:
        Full unit name with abbreviation
        
    Raises:
        ValueError: If unit is not supported
    """
    validate_unit(unit)
    
    unit_names = {
        'pM': 'picoMolar (pM)',
        'nM': 'nanoMolar (nM)',
        'uM': 'microMolar (uM)', 
        'mM': 'milliMolar (mM)',
        'M': 'Molar (M)'
    }
    
    return unit_names[unit]
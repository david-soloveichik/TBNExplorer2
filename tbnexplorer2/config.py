"""Configuration management for TBNExplorer2."""

import os
from pathlib import Path

# Try to load .env file if it exists
env_file = Path(__file__).parent.parent / ".env"
if env_file.exists():
    with open(env_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                key, value = line.split("=", 1)
                os.environ.setdefault(key.strip(), value.strip())

# Configuration values
NORMALIZ_PATH = os.environ.get("NORMALIZ_PATH", "normaliz")
COFFEE_CLI_PATH = os.environ.get("COFFEE_CLI_PATH", "coffee-cli")
FOURTI2_PATH = os.environ.get("FOURTI2_PATH", "4ti2")

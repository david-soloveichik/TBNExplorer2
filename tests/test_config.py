import os
from unittest.mock import patch

from tbnexplorer2 import config


class TestConfig:
    def test_normaliz_path_from_env(self):
        """Test NORMALIZ_PATH from environment variable."""
        test_path = "/custom/path/to/normaliz"
        with patch.dict(os.environ, {"NORMALIZ_PATH": test_path}):
            import importlib

            importlib.reload(config)
            assert test_path == config.NORMALIZ_PATH

    def test_coffee_cli_path_from_env(self):
        """Test COFFEE_CLI_PATH from environment variable."""
        test_path = "/custom/path/to/coffee-cli"
        with patch.dict(os.environ, {"COFFEE_CLI_PATH": test_path}):
            import importlib

            importlib.reload(config)
            assert test_path == config.COFFEE_CLI_PATH

    def test_fourti2_path_from_env(self):
        """Test FOURTI2_PATH from environment variable."""
        test_path = "/custom/path/to/4ti2"
        with patch.dict(os.environ, {"FOURTI2_PATH": test_path}):
            import importlib

            importlib.reload(config)
            assert test_path == config.FOURTI2_PATH

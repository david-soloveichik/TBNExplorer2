import os
import tempfile
from pathlib import Path
from unittest.mock import patch

from tbnexplorer2 import config


class TestConfig:
    def test_default_normaliz_path(self):
        """Test default NORMALIZ_PATH value."""
        with patch.dict(os.environ, {}, clear=True):
            # Reload the module to get fresh configuration
            import importlib

            importlib.reload(config)
            assert config.NORMALIZ_PATH == "normaliz"

    def test_normaliz_path_from_env(self):
        """Test NORMALIZ_PATH from environment variable."""
        test_path = "/custom/path/to/normaliz"
        with patch.dict(os.environ, {"NORMALIZ_PATH": test_path}):
            import importlib

            importlib.reload(config)
            assert test_path == config.NORMALIZ_PATH

    def test_default_coffee_cli_path(self):
        """Test default COFFEE_CLI_PATH value."""
        with patch.dict(os.environ, {}, clear=True):
            import importlib

            importlib.reload(config)
            assert config.COFFEE_CLI_PATH == "coffee-cli"

    def test_coffee_cli_path_from_env(self):
        """Test COFFEE_CLI_PATH from environment variable."""
        test_path = "/custom/path/to/coffee-cli"
        with patch.dict(os.environ, {"COFFEE_CLI_PATH": test_path}):
            import importlib

            importlib.reload(config)
            assert test_path == config.COFFEE_CLI_PATH

    def test_default_fourti2_path(self):
        """Test default FOURTI2_PATH value."""
        with patch.dict(os.environ, {}, clear=True):
            import importlib

            importlib.reload(config)
            assert config.FOURTI2_PATH == "/Users/dsolov/Documents/ResearchTools/4ti2"

    def test_fourti2_path_from_env(self):
        """Test FOURTI2_PATH from environment variable."""
        test_path = "/custom/path/to/4ti2"
        with patch.dict(os.environ, {"FOURTI2_PATH": test_path}):
            import importlib

            importlib.reload(config)
            assert test_path == config.FOURTI2_PATH

    def test_env_file_loading(self):
        """Test loading configuration from .env file."""
        # Create a temporary .env file
        with tempfile.TemporaryDirectory() as temp_dir:
            env_content = """
# Test environment file
NORMALIZ_PATH=/test/normaliz
COFFEE_CLI_PATH=/test/coffee
FOURTI2_PATH=/test/4ti2
"""
            env_file = Path(temp_dir) / ".env"
            env_file.write_text(env_content)

            # Mock the env file path
            with patch("tbnexplorer2.config.Path") as mock_path:
                mock_path.return_value.parent.parent = Path(temp_dir)
                mock_path_instance = mock_path.return_value
                mock_path_instance.parent.parent.__truediv__ = lambda self, x: env_file

                # Clear environment and reload
                with patch.dict(os.environ, {}, clear=True):
                    import importlib

                    importlib.reload(config)

                    # These will be set from the .env file
                    assert os.environ.get("NORMALIZ_PATH") == "/test/normaliz"
                    assert os.environ.get("COFFEE_CLI_PATH") == "/test/coffee"
                    assert os.environ.get("FOURTI2_PATH") == "/test/4ti2"

    def test_env_file_with_comments_and_whitespace(self):
        """Test .env file parsing with comments and whitespace."""
        with tempfile.TemporaryDirectory() as temp_dir:
            env_content = """
# Comment line

NORMALIZ_PATH = /path/with/spaces  
# Another comment
  COFFEE_CLI_PATH=/another/path
"""
            env_file = Path(temp_dir) / ".env"
            env_file.write_text(env_content)

            with patch("tbnexplorer2.config.Path") as mock_path:
                mock_path.return_value.parent.parent = Path(temp_dir)
                mock_path_instance = mock_path.return_value
                mock_path_instance.parent.parent.__truediv__ = lambda self, x: env_file

                with patch.dict(os.environ, {}, clear=True):
                    import importlib

                    importlib.reload(config)

                    assert os.environ.get("NORMALIZ_PATH") == "/path/with/spaces"
                    assert os.environ.get("COFFEE_CLI_PATH") == "/another/path"

    def test_env_file_not_exists(self):
        """Test behavior when .env file doesn't exist."""
        with patch("tbnexplorer2.config.Path") as mock_path:
            # Make the env file not exist
            mock_path.return_value.parent.parent.__truediv__ = lambda self, x: Path("/nonexistent/.env")

            with patch.dict(os.environ, {}, clear=True):
                import importlib

                importlib.reload(config)

                # Should use defaults
                assert config.NORMALIZ_PATH == "normaliz"
                assert config.COFFEE_CLI_PATH == "coffee-cli"

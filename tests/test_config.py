import pytest

from optimhc.core.config import DEFAULT_CONFIG, Config, _deep_merge


class TestConfig:
    """Test cases for Config class."""

    def test_default_config_creation(self):
        """Test creating config with default values."""
        config = Config()
        assert config["inputType"] == "pepxml"
        assert config["outputDir"] == "./results"
        assert config["rescore"]["testFDR"] == 0.01
        assert config["rescore"]["model"] == "Percolator"

    def test_config_from_dict(self):
        """Test creating config from dictionary."""
        test_dict = {
            "inputType": "pin",
            "inputFile": ["test.pin"],
            "outputDir": "./test_results",
        }
        config = Config(test_dict)
        assert config["inputType"] == "pin"
        assert config["inputFile"] == ["test.pin"]
        assert config["outputDir"] == "./test_results"
        # Should still have default values for unspecified keys
        assert config["rescore"]["testFDR"] == 0.01

    def test_config_dict_access(self):
        """Test dictionary-like access methods."""
        config = Config()

        # Test __getitem__
        assert config["inputType"] == "pepxml"

        # Test __setitem__
        config["inputType"] = "pin"
        assert config["inputType"] == "pin"

        # Test get method
        assert config.get("inputType") == "pin"
        assert config.get("nonexistent", "default") == "default"

        # Test __contains__
        assert "inputType" in config
        assert "nonexistent" not in config

    def test_config_validation_missing_required(self):
        """Test validation fails with missing required fields."""
        config = Config()
        # Remove required fields to test missing field validation
        del config._config["inputFile"]
        del config._config["outputDir"]
        del config._config["rescore"]

        with pytest.raises(ValueError, match="Missing required configuration"):
            config.validate()

    def test_config_validation_invalid_input_type(self):
        """Test validation fails with invalid inputType."""
        invalid_config = {
            "inputType": "invalid",
            "inputFile": ["test.xml"],
            "outputDir": "./output",
            "rescore": {"testFDR": 0.01},
        }
        config = Config(invalid_config)

        with pytest.raises(ValueError, match="inputType must be 'pepxml' or 'pin'"):
            config.validate()

    def test_config_validation_empty_input_file(self):
        """Test validation fails with empty inputFile."""
        invalid_config = {
            "inputType": "pepxml",
            "inputFile": [],
            "outputDir": "./output",
            "rescore": {"testFDR": 0.01},
        }
        config = Config(invalid_config)

        with pytest.raises(ValueError, match="inputFile list cannot be empty"):
            config.validate()


class TestConfigUtilities:
    """Test utility functions."""

    def test_deep_merge(self):
        """Test deep merge functionality."""
        default = {"a": 1, "b": {"c": 2, "d": 3}, "e": [1, 2, 3]}
        override = {"b": {"c": 20}, "f": 4}

        result = _deep_merge(default, override)

        assert result["a"] == 1  # unchanged
        assert result["b"]["c"] == 20  # overridden
        assert result["b"]["d"] == 3  # preserved
        assert result["e"] == [1, 2, 3]  # unchanged
        assert result["f"] == 4  # added

    def test_config_to_dict(self):
        """Test converting config to dictionary."""
        config = Config()
        config_dict = config.to_dict()

        assert isinstance(config_dict, dict)
        assert config_dict["inputType"] == "pepxml"
        assert config_dict["rescore"]["testFDR"] == 0.01


class TestDefaultConfig:
    """Test default configuration."""

    def test_default_config_structure(self):
        """Test that DEFAULT_CONFIG has expected structure."""
        assert "inputType" in DEFAULT_CONFIG
        assert "inputFile" in DEFAULT_CONFIG
        assert "outputDir" in DEFAULT_CONFIG
        assert "rescore" in DEFAULT_CONFIG

        assert isinstance(DEFAULT_CONFIG["rescore"], dict)
        assert "testFDR" in DEFAULT_CONFIG["rescore"]
        assert "model" in DEFAULT_CONFIG["rescore"]

    def test_default_config_values(self):
        """Test default configuration values."""
        assert DEFAULT_CONFIG["inputType"] == "pepxml"
        assert DEFAULT_CONFIG["outputDir"] == "./results"
        assert DEFAULT_CONFIG["decoyPrefix"] == "DECOY_"
        assert DEFAULT_CONFIG["numProcess"] == 4
        assert DEFAULT_CONFIG["rescore"]["testFDR"] == 0.01
        assert DEFAULT_CONFIG["rescore"]["model"] == "Percolator"

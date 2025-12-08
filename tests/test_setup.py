"""
Minimal test file to verify pytest setup is working correctly.
This file contains basic tests that should always pass to confirm the test framework is configured properly.
"""

import pytest
import sys
import os

# Add src directory to path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


#def test_pytest_is_working():
#    """Basic test to verify pytest is functioning."""
#    assert True


def test_python_version():
    """Test that we're running on a supported Python version."""
    assert sys.version_info >= (3, 9), "Python 3.9+ is required"


def test_src_imports():
    """Test that we can import from the src directory."""
    try:
        from src import classes
        from src import utils
        assert True
    except ImportError as e:
        pytest.fail(f"Failed to import from src: {e}")


def test_data_directory_exists():
    """Test that the test data directory exists."""
    test_data_dir = os.path.join(os.path.dirname(__file__), 'data')
    assert os.path.exists(test_data_dir), "Test data directory should exist"
    assert os.path.isdir(test_data_dir), "Test data path should be a directory"


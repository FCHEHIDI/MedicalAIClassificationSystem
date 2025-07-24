# Test Configuration
import pytest
from pathlib import Path
import sys

# Add project root to Python path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Test settings
API_BASE_URL = "http://localhost:8001"
MODELS_DIR = PROJECT_ROOT / "models"
DATA_DIR = PROJECT_ROOT / "data"

# Test constants
TEST_TIMEOUT = 30  # seconds
EXPECTED_ACCURACY_THRESHOLD = 0.90  # 90% minimum for tests to pass

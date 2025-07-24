# Tests Package
from pathlib import Path

# Add project root to path for test imports
import sys
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

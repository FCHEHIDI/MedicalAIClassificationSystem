#!/usr/bin/env python3
"""
Test Runner for Medical Classification Engine
===========================================
Simple script to run all tests and display results.
"""

import subprocess
import sys
from pathlib import Path

def run_tests():
    """Run all tests and display results"""
    print("🧪 Medical Classification Engine - Test Suite")
    print("=" * 60)
    
    # Change to project root
    project_root = Path(__file__).parent.parent
    
    # Run pytest with detailed output
    cmd = [
        sys.executable, "-m", "pytest",
        "tests/",
        "-v",
        "--tb=short",
        "--color=yes",
        "-s"  # Don't capture print statements
    ]
    
    try:
        result = subprocess.run(cmd, cwd=project_root, capture_output=False)
        return result.returncode == 0
    except Exception as e:
        print(f"❌ Error running tests: {e}")
        return False

def run_specific_test(test_name):
    """Run a specific test file"""
    project_root = Path(__file__).parent.parent
    
    cmd = [
        sys.executable, "-m", "pytest",
        f"tests/{test_name}",
        "-v",
        "--tb=short",
        "--color=yes",
        "-s"
    ]
    
    try:
        result = subprocess.run(cmd, cwd=project_root, capture_output=False)
        return result.returncode == 0
    except Exception as e:
        print(f"❌ Error running test {test_name}: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Run specific test
        test_file = sys.argv[1]
        if not test_file.startswith("test_"):
            test_file = f"test_{test_file}"
        if not test_file.endswith(".py"):
            test_file = f"{test_file}.py"
            
        print(f"🎯 Running specific test: {test_file}")
        success = run_specific_test(test_file)
    else:
        # Run all tests
        print("🚀 Running all tests...")
        success = run_tests()
    
    if success:
        print("\n✅ All tests passed!")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed!")
        sys.exit(1)

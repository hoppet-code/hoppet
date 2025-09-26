#!/usr/bin/env python3
"""
Test script for hoppet Python interface import and basic functionality.

This script tests:
1. Basic import of the hoppet module
2. Availability of key functions like Start()
3. Basic functionality test with Start/DeleteAll
"""

import sys
import traceback

def main():
    print("Python path:", sys.path)
    print("Testing hoppet import...")
    
    try:
        import hoppet as hp
        print("PASS: Import successful!")
        
        # List available functions (first 10)
        available_functions = sorted([x for x in dir(hp) if not x.startswith('_')])
        print("Available functions:", available_functions[:10])
        
        # Test for Start function specifically
        if hasattr(hp, 'Start'):
            print("PASS: Start function found!")
            
            # Test basic functionality
            hp.Start(0.1, 2)
            print("PASS: Start function works!")
            
            # Clean up
            hp.DeleteAll()
            print("PASS: DeleteAll works!")
            
        else:
            print("FAIL: Start function NOT found!")
            start_like = [x for x in available_functions if 'start' in x.lower() or 'Start' in x]
            print("Functions containing 'start':", start_like)
            return 1
            
    except ImportError as e:
        print("FAIL: Import failed:", e)
        print("This suggests the hoppet module was not built or is not in the Python path")
        traceback.print_exc()
        return 1
        
    except Exception as e:
        print("FAIL: Error during testing:", e)
        traceback.print_exc()
        return 1
    
    print("PASS: All tests passed!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
"""
Edge Case Testing Script
========================
Tests the specific edge cases that were failing in our previous test suite
to validate the improved model performance.
"""

import requests
import json
import time

# API endpoint
API_URL = "http://localhost:8001"

# Test cases that were failing
TEST_CASES = [
    {
        "name": "Cardiology Edge Case (was failing)",
        "text": "Patient with chest pain and shortness of breath, but also reports dizziness and palpitations. ECG shows minimal changes. Troponin slightly elevated. Stress test pending.",
        "expected": "Cardiology",
        "previous_result": "Gastroenterology (0.272 confidence)"
    },
    {
        "name": "Dermatology Edge Case (was failing)", 
        "text": "Rapidly spreading erythematous rash with fever and malaise following antibiotic use. Suspected drug reaction versus Stevens-Johnson syndrome. Dermatology emergency consultation.",
        "expected": "Dermatology",
        "previous_result": "Emergency (0.489 confidence)"
    },
    # Additional similar cases to test generalization
    {
        "name": "Similar Cardiology Edge Case",
        "text": "67-year-old male with palpitations and lightheadedness during exercise. Holter monitor shows episodes of supraventricular tachycardia. Electrophysiology study recommended.",
        "expected": "Cardiology",
        "previous_result": "N/A (new test)"
    },
    {
        "name": "Similar Dermatology Edge Case",
        "text": "Drug-induced hypersensitivity reaction with widespread skin rash, fever, and constitutional symptoms following antibiotic therapy. Dermatologic consultation required.",
        "expected": "Dermatology", 
        "previous_result": "N/A (new test)"
    }
]

def test_api_endpoint():
    """Test if API is running"""
    try:
        response = requests.get(f"{API_URL}/health", timeout=5)
        return response.status_code == 200
    except:
        return False

def make_prediction(text):
    """Make prediction via API"""
    try:
        response = requests.post(
            f"{API_URL}/predict",
            json={"text": text},
            timeout=10
        )
        if response.status_code == 200:
            return response.json(), None
        else:
            return None, f"API Error: {response.status_code}"
    except requests.exceptions.ConnectionError:
        return None, "Cannot connect to API"
    except Exception as e:
        return None, f"Error: {str(e)}"

def run_edge_case_tests():
    """Run edge case validation tests"""
    
    print("🧪 Edge Case Validation Tests")
    print("=" * 50)
    
    # Check API status
    if not test_api_endpoint():
        print("❌ API is not running. Please start the API first:")
        print("   python start_api.py")
        return
    
    print("✅ API is running")
    print()
    
    results = []
    
    for i, test_case in enumerate(TEST_CASES, 1):
        print(f"🔍 Test {i}: {test_case['name']}")
        print(f"📝 Text: {test_case['text'][:100]}...")
        print(f"🎯 Expected: {test_case['expected']}")
        print(f"📊 Previous: {test_case['previous_result']}")
        
        # Make prediction
        result, error = make_prediction(test_case['text'])
        
        if result:
            predicted = result['specialty']
            confidence = result['confidence']
            correct = predicted.lower() == test_case['expected'].lower()
            
            print(f"🤖 Predicted: {predicted} ({confidence:.1%} confidence)")
            
            if correct:
                if confidence >= 0.7:
                    print("✅ PASS - Correct with high confidence")
                    status = "PASS_HIGH"
                else:
                    print("⚠️  PASS - Correct but low confidence")
                    status = "PASS_LOW"
            else:
                print("❌ FAIL - Incorrect prediction")
                status = "FAIL"
            
            results.append({
                'test': test_case['name'],
                'expected': test_case['expected'],
                'predicted': predicted,
                'confidence': confidence,
                'correct': correct,
                'status': status
            })
        else:
            print(f"❌ ERROR: {error}")
            results.append({
                'test': test_case['name'],
                'expected': test_case['expected'],
                'predicted': 'ERROR',
                'confidence': 0.0,
                'correct': False,
                'status': 'ERROR'
            })
        
        print("-" * 50)
        time.sleep(0.5)  # Small delay between tests
    
    # Summary
    print("\n📊 EDGE CASE TEST SUMMARY")
    print("=" * 30)
    
    total_tests = len(results)
    correct_tests = sum(1 for r in results if r['correct'])
    high_conf_tests = sum(1 for r in results if r['status'] == 'PASS_HIGH')
    
    print(f"Total Tests: {total_tests}")
    print(f"Correct Predictions: {correct_tests}/{total_tests} ({correct_tests/total_tests*100:.1f}%)")
    print(f"High Confidence Correct: {high_conf_tests}/{total_tests} ({high_conf_tests/total_tests*100:.1f}%)")
    
    print("\nDetailed Results:")
    for result in results:
        emoji = "✅" if result['status'] == 'PASS_HIGH' else "⚠️" if result['status'] == 'PASS_LOW' else "❌"
        print(f"{emoji} {result['test']}: {result['predicted']} ({result['confidence']:.1%})")
    
    # Improvement analysis
    print("\n🎯 IMPROVEMENT ANALYSIS")
    print("-" * 25)
    
    previously_failing = [r for r in results if "was failing" in r['test']]
    if previously_failing:
        fixed_count = sum(1 for r in previously_failing if r['correct'])
        print(f"Previously failing cases: {len(previously_failing)}")
        print(f"Now fixed: {fixed_count}/{len(previously_failing)}")
        
        if fixed_count == len(previously_failing):
            print("🎉 ALL PREVIOUSLY FAILING CASES ARE NOW FIXED!")
        elif fixed_count > 0:
            print(f"✅ {fixed_count} cases improved, {len(previously_failing)-fixed_count} still need work")
        else:
            print("❌ No improvement in previously failing cases")
    
    # Recommendations
    print("\n💡 RECOMMENDATIONS")
    print("-" * 18)
    
    failing_cases = [r for r in results if not r['correct']]
    low_conf_cases = [r for r in results if r['status'] == 'PASS_LOW']
    
    if not failing_cases and not low_conf_cases:
        print("🎯 PERFECT! All edge cases passing with high confidence.")
        print("✅ Ready for production deployment!")
    elif failing_cases:
        print("🔄 Need more training data for:")
        for case in failing_cases:
            print(f"   - {case['test']}")
    elif low_conf_cases:
        print("⚡ Consider adding more similar cases to boost confidence:")
        for case in low_conf_cases:
            print(f"   - {case['test']}")
    
    return results

if __name__ == "__main__":
    results = run_edge_case_tests()

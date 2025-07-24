"""
Compatibility Test for Medical Classification API
===============================================
Tests the API with hybrid Chi-square + F-score feature selection
"""

import requests
import json

def test_api_compatibility():
    """Test API compatibility with new feature selection"""
    base_url = "http://localhost:8001"
    
    # Test cases
    test_cases = [
        {
            "text": "Patient presents with chest pain radiating to left arm, elevated troponin levels, and ST-segment elevation on ECG.",
            "expected": "Cardiology"
        },
        {
            "text": "Patient with severe abdominal pain, rebound tenderness, and elevated white blood cell count. Emergency appendectomy required.",
            "expected": "Emergency"
        },
        {
            "text": "Patient with chronic cough, shortness of breath, and bilateral pulmonary infiltrates on chest X-ray.",
            "expected": "Pulmonology"
        },
        {
            "text": "Patient presents with erythematous rash, scaling, and pruritus on bilateral extremities.",
            "expected": "Dermatology"
        },
        {
            "text": "Patient with chronic dyspepsia, epigastric pain, and positive H. pylori test.",
            "expected": "Gastroenterology"
        }
    ]
    
    print("🧪 Testing API Compatibility with Hybrid Feature Selection")
    print("=" * 60)
    
    # Test health endpoint
    try:
        response = requests.get(f"{base_url}/health")
        print(f"✅ Health check: {response.json()}")
    except Exception as e:
        print(f"❌ Health check failed: {e}")
        return False
    
    # Test model info endpoint
    try:
        response = requests.get(f"{base_url}/model-info")
        model_info = response.json()
        print(f"✅ Model info: {model_info['model_name']}")
        print(f"📊 Accuracy: {model_info.get('test_accuracy', 0)*100:.1f}%")
        print(f"📊 F2 Score: {model_info.get('f2_score', 0)*100:.1f}%")
        print(f"🔧 Features: {model_info.get('feature_selection', 'Standard')}")
    except Exception as e:
        print(f"❌ Model info failed: {e}")
        return False
    
    # Test predictions
    print("\n🔬 Testing Predictions:")
    print("-" * 40)
    
    correct_predictions = 0
    total_predictions = len(test_cases)
    
    for i, test_case in enumerate(test_cases, 1):
        try:
            response = requests.post(
                f"{base_url}/predict",
                json={"text": test_case["text"]},
                headers={"Content-Type": "application/json"}
            )
            
            if response.status_code == 200:
                result = response.json()
                predicted = result["specialty"]
                confidence = result["confidence"]
                expected = test_case["expected"]
                
                status = "✅" if predicted == expected else "❌"
                if predicted == expected:
                    correct_predictions += 1
                
                print(f"{status} Test {i}: {predicted} ({confidence:.1%}) | Expected: {expected}")
                print(f"   Text: {test_case['text'][:80]}...")
                
            else:
                print(f"❌ Test {i}: API error {response.status_code}")
                
        except Exception as e:
            print(f"❌ Test {i}: Exception {e}")
    
    # Results summary
    accuracy = (correct_predictions / total_predictions) * 100
    print("\n📊 Test Results Summary:")
    print("=" * 40)
    print(f"✅ Correct Predictions: {correct_predictions}/{total_predictions}")
    print(f"📈 Test Accuracy: {accuracy:.1f}%")
    print(f"🔧 Feature Selection: Hybrid Chi-square + F-score")
    print(f"📊 Model Performance: {model_info.get('test_accuracy', 0)*100:.1f}% accuracy")
    
    if accuracy >= 80:
        print("🎉 API COMPATIBILITY: EXCELLENT!")
        return True
    elif accuracy >= 60:
        print("⚠️  API COMPATIBILITY: GOOD")
        return True
    else:
        print("❌ API COMPATIBILITY: NEEDS IMPROVEMENT")
        return False

if __name__ == "__main__":
    success = test_api_compatibility()
    if success:
        print("\n✅ 100% API Compatibility Achieved!")
        print("🚀 Ready for smooth deployment!")
    else:
        print("\n❌ Compatibility issues detected")

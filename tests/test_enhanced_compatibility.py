"""
Enhanced Medical Classification Compatibility Test
=================================================
Tests with realistic medical cases matching PubMed training data style
"""

import requests
import json

def test_enhanced_compatibility():
    """Test API with realistic medical cases"""
    base_url = "http://localhost:8001"
    
    # Realistic test cases matching PubMed medical literature style
    test_cases = [
        {
            "text": "A 58-year-old male presented with acute onset chest pain with ST-segment elevation on electrocardiogram. Coronary angiography revealed complete occlusion of the left anterior descending artery. Percutaneous coronary intervention with stent placement was performed.",
            "expected": "Cardiology"
        },
        {
            "text": "A 45-year-old female presented to the emergency department following a motor vehicle collision with multiple trauma including rib fractures, pneumothorax, and hemodynamic instability requiring immediate surgical intervention.",
            "expected": "Emergency"
        },
        {
            "text": "A 62-year-old male with a history of chronic obstructive pulmonary disease presented with acute exacerbation. Chest X-ray showed bilateral infiltrates and arterial blood gas analysis revealed respiratory acidosis requiring mechanical ventilation.",
            "expected": "Pulmonology"
        },
        {
            "text": "A 55-year-old female underwent upper endoscopy for evaluation of chronic dyspepsia. Biopsy specimens revealed Helicobacter pylori-associated gastritis with intestinal metaplasia.",
            "expected": "Gastroenterology"
        },
        {
            "text": "A 40-year-old patient presented with a pigmented lesion on the back showing asymmetry and irregular borders. Dermatoscopic examination and subsequent biopsy confirmed melanoma with Breslow thickness of 1.2 mm.",
            "expected": "Dermatology"
        },
        {
            "text": "Echocardiographic assessment demonstrated severe aortic stenosis with mean gradient of 45 mmHg. Left ventricular ejection fraction was preserved at 60%. Aortic valve replacement was recommended.",
            "expected": "Cardiology"
        },
        {
            "text": "Patient with acute respiratory failure secondary to pneumonia was admitted to the intensive care unit. Bronchoscopy with bronchoalveolar lavage was performed for microbiologic diagnosis.",
            "expected": "Pulmonology"
        },
        {
            "text": "Colonoscopy revealed multiple adenomatous polyps throughout the colon. Histopathologic examination showed high-grade dysplasia requiring total colectomy.",
            "expected": "Gastroenterology"
        }
    ]
    
    print("🧪 Enhanced Medical Classification Compatibility Test")
    print("=" * 60)
    print("🎯 Using realistic PubMed-style medical cases")
    
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
        print(f"✅ Model: {model_info['model_name']}")
        print(f"📊 Accuracy: {model_info.get('test_accuracy', 0)*100:.1f}%")
        print(f"📊 F2 Score: {model_info.get('f2_score', 0)*100:.1f}%")
        print(f"📊 OOB Score: {model_info.get('oob_score', 0)*100:.1f}%")
        print(f"🔧 Features: {model_info.get('feature_selection', 'Standard')}")
        print(f"📚 Training Data: {model_info.get('training_size', 0)} samples")
    except Exception as e:
        print(f"❌ Model info failed: {e}")
        return False
    
    # Test predictions
    print("\n🔬 Testing Enhanced Medical Cases:")
    print("-" * 50)
    
    correct_predictions = 0
    total_predictions = len(test_cases)
    results = []
    
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
                
                results.append({
                    "test": i,
                    "predicted": predicted,
                    "expected": expected,
                    "confidence": confidence,
                    "correct": predicted == expected
                })
                
                print(f"{status} Test {i}: {predicted} ({confidence:.1%}) | Expected: {expected}")
                print(f"   📝 Medical text: {test_case['text'][:100]}...")
                
            else:
                print(f"❌ Test {i}: API error {response.status_code}")
                
        except Exception as e:
            print(f"❌ Test {i}: Exception {e}")
    
    # Detailed Results Analysis
    accuracy = (correct_predictions / total_predictions) * 100
    print("\n📊 Enhanced Test Results Summary:")
    print("=" * 50)
    print(f"✅ Correct Predictions: {correct_predictions}/{total_predictions}")
    print(f"📈 Test Accuracy: {accuracy:.1f}%")
    print(f"🔧 Feature Selection: Fine-tuned Hybrid")
    print(f"📊 Model Performance: {model_info.get('test_accuracy', 0)*100:.1f}% accuracy")
    print(f"📚 Training Scale: {model_info.get('training_size', 0)} samples")
    
    # Specialty breakdown
    specialty_results = {}
    for result in results:
        specialty = result["expected"]
        if specialty not in specialty_results:
            specialty_results[specialty] = {"correct": 0, "total": 0}
        specialty_results[specialty]["total"] += 1
        if result["correct"]:
            specialty_results[specialty]["correct"] += 1
    
    print("\n📋 Performance by Specialty:")
    print("-" * 30)
    for specialty, stats in specialty_results.items():
        specialty_acc = (stats["correct"] / stats["total"]) * 100
        print(f"   {specialty}: {stats['correct']}/{stats['total']} ({specialty_acc:.0f}%)")
    
    if accuracy >= 80:
        print("\n🎉 EXCELLENT API COMPATIBILITY!")
        return True
    elif accuracy >= 60:
        print("\n✅ GOOD API COMPATIBILITY")
        return True
    else:
        print("\n⚠️  API COMPATIBILITY NEEDS IMPROVEMENT")
        return False

if __name__ == "__main__":
    success = test_enhanced_compatibility()
    if success:
        print("\n🚀 Ready for production deployment!")
    else:
        print("\n🔧 Consider additional model tuning")

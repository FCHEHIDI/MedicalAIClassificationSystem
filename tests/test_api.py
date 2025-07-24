"""
Medical API Test Client
======================

Simple test client to demonstrate the medical classification API.
"""

import json
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

print("🧪 Medical Classification API Test Client")
print("=" * 50)

# Test medical texts
test_cases = [
    {
        "text": "Patient presents with chest pain radiating to left arm, elevated troponin levels, and ST-segment elevation on ECG. Urgent cardiac catheterization recommended.",
        "expected": "Cardiology"
    },
    {
        "text": "Trauma patient with multiple rib fractures, pneumothorax, and hemodynamic instability following motor vehicle accident. Emergency surgical intervention required.",
        "expected": "Emergency"
    },
    {
        "text": "Patient with chronic cough, shortness of breath, and bilateral pulmonary infiltrates on chest X-ray. Bronchoscopy reveals inflammatory changes.",
        "expected": "Pulmonology"
    },
    {
        "text": "Patient presents with abdominal pain, bloody diarrhea, and weight loss. Colonoscopy shows ulcerative changes consistent with inflammatory bowel disease.",
        "expected": "Gastroenterology"
    },
    {
        "text": "Suspicious pigmented lesion on patient's back with irregular borders, asymmetry, and recent changes. Dermoscopy indicates possible malignant melanoma.",
        "expected": "Dermatology"
    }
]

try:
    # Test local model prediction (without API server)
    print("🔬 Testing local model predictions...")
    
    import joblib
    from pathlib import Path
    
    # Load models
    models_dir = Path("models")
    model = joblib.load(models_dir / "regularized_medical_classifier.joblib")
    vectorizer = joblib.load(models_dir / "regularized_tfidf_vectorizer.joblib")
    feature_selector = joblib.load(models_dir / "medical_feature_selector.joblib")
    label_encoder = joblib.load(models_dir / "regularized_label_encoder.joblib")
    
    print("✅ Models loaded successfully")
    
    correct_predictions = 0
    
    for i, test_case in enumerate(test_cases):
        print(f"\n📋 Test Case {i+1}: {test_case['expected']}")
        print(f"   Text: {test_case['text'][:80]}...")
        
        # Process text
        text_tfidf = vectorizer.transform([test_case['text']])
        text_features = feature_selector.transform(text_tfidf)
        
        # Predict
        prediction = model.predict(text_features)[0]
        probabilities = model.predict_proba(text_features)[0]
        
        predicted_specialty = label_encoder.inverse_transform([prediction])[0]
        confidence = max(probabilities)
        
        # Check if correct
        is_correct = predicted_specialty == test_case['expected']
        if is_correct:
            correct_predictions += 1
        
        status = "✅" if is_correct else "❌"
        print(f"   {status} Predicted: {predicted_specialty} (confidence: {confidence:.3f})")
        
        # Show top predictions
        specialty_scores = {}
        for j, prob in enumerate(probabilities):
            specialty = label_encoder.inverse_transform([j])[0]
            specialty_scores[specialty] = prob
        
        sorted_scores = sorted(specialty_scores.items(), key=lambda x: x[1], reverse=True)
        print(f"   Top predictions:")
        for k, (specialty, score) in enumerate(sorted_scores[:3]):
            print(f"     {k+1}. {specialty}: {score:.3f}")
    
    accuracy = correct_predictions / len(test_cases)
    print(f"\n🎯 Test Results:")
    print(f"   Correct predictions: {correct_predictions}/{len(test_cases)}")
    print(f"   Accuracy: {accuracy:.1%}")
    
    if accuracy >= 0.8:
        print(f"   ✅ EXCELLENT - Model performing well on test cases")
    elif accuracy >= 0.6:
        print(f"   ⚠️  GOOD - Model shows reasonable performance")
    else:
        print(f"   🚨 NEEDS IMPROVEMENT - Consider retraining")
    
    # Generate API request examples
    print(f"\n📝 Example API Requests (for when server is running):")
    
    api_examples = []
    for test_case in test_cases[:2]:  # Show first 2 examples
        request_data = {
            "text": test_case['text'],
            "include_confidence_scores": True,
            "minimum_confidence": 0.1
        }
        api_examples.append(request_data)
    
    print(f"\nCURL example:")
    print(f'curl -X POST "http://localhost:8000/classify" \\')
    print(f'  -H "Content-Type: application/json" \\')
    print(f'  -d \'{json.dumps(api_examples[0], indent=2)}\'')
    
    print(f"\n💡 To start the API server:")
    print(f"   python start_api.py")
    print(f"\n📖 Then visit: http://localhost:8000/docs")
    
except ImportError as e:
    print(f"❌ Missing dependencies: {e}")
    print("💡 Install with: pip install joblib scikit-learn")
    
except FileNotFoundError as e:
    print(f"❌ Model files not found: {e}")
    print("💡 Run professional_model_training.py first")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print(f"\n" + "=" * 50)

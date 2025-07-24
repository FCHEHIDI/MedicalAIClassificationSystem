#!/usr/bin/env python3
"""
Test script to check if models are properly fitted
"""
import joblib
import json
from pathlib import Path

def test_models():
    """Test if models are properly loaded and fitted"""
    try:
        # Use relative path from tests directory
        models_dir = Path(__file__).parent.parent / "models"
        
        print("🔍 Testing regularized models...")
        
        # Load regularized models
        model = joblib.load(models_dir / "regularized_medical_classifier.joblib")
        vectorizer = joblib.load(models_dir / "regularized_tfidf_vectorizer.joblib")
        feature_selector = joblib.load(models_dir / "medical_feature_selector.joblib")
        label_encoder = joblib.load(models_dir / "regularized_label_encoder.joblib")
        
        print("✅ Models loaded successfully")
        
        # Check if vectorizer is fitted
        print(f"📊 Vectorizer vocabulary size: {len(vectorizer.vocabulary_) if hasattr(vectorizer, 'vocabulary_') else 'No vocabulary'}")
        print(f"📊 Vectorizer features: {vectorizer.max_features}")
        print(f"📊 Vectorizer fitted: {hasattr(vectorizer, 'vocabulary_') and len(vectorizer.vocabulary_) > 0}")
        
        # Test prediction
        test_text = "Patient presents with chest pain and shortness of breath. ECG shows irregular rhythm."
        print(f"\n🧪 Testing with: '{test_text[:50]}...'")
        
        # Transform text
        text_tfidf = vectorizer.transform([test_text])
        print(f"✅ Text transformation successful. Shape: {text_tfidf.shape}")
        
        # Feature selection
        text_features = feature_selector.transform(text_tfidf)
        print(f"✅ Feature selection successful. Shape: {text_features.shape}")
        
        # Prediction
        prediction = model.predict(text_features)[0]
        probabilities = model.predict_proba(text_features)[0]
        predicted_specialty = label_encoder.inverse_transform([prediction])[0]
        
        print(f"✅ Prediction successful: {predicted_specialty}")
        print(f"📊 Confidence: {probabilities.max():.3f}")
        
        # Load model info
        with open(models_dir / "regularized_model_info.json", 'r') as f:
            model_info = json.load(f)
        
        print(f"\n📋 Model Info:")
        print(f"   - Accuracy: {model_info.get('accuracy', 'N/A')}")
        print(f"   - Features: {model_info.get('n_features', 'N/A')}")
        print(f"   - Training samples: {model_info.get('n_samples', 'N/A')}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        print(f"📋 Full traceback:\n{traceback.format_exc()}")
        return False

if __name__ == "__main__":
    success = test_models()
    print(f"\n{'✅ All tests passed!' if success else '❌ Tests failed!'}")

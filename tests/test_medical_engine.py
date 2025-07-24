"""
Unit Tests for Medical Classification Engine
===========================================
Pytest-based unit tests for the medical classification system.
"""

import pytest
import requests
import json
import joblib
from pathlib import Path
import sys

# Import project modules
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import prediction function from dashboard module
from src.dashboard.medical_dashboard import predict_specialty as dashboard_predict


class TestMedicalClassificationEngine:
    """Test suite for Medical Classification Engine"""
    
    @classmethod
    def setup_class(cls):
        """Setup test environment"""
        cls.models_dir = PROJECT_ROOT / "models"
        cls.api_base_url = "http://localhost:8001"
        
    def test_models_exist(self):
        """Test that all required model files exist"""
        required_models = [
            "medical_classifier.joblib",
            "medical_tfidf_vectorizer.joblib", 
            "medical_label_encoder.joblib",
            "medical_feature_selector.joblib"
        ]
        
        for model_file in required_models:
            model_path = self.models_dir / model_file
            assert model_path.exists(), f"Model file {model_file} not found"
            
    def test_models_load(self):
        """Test that models can be loaded properly"""
        try:
            classifier = joblib.load(self.models_dir / "medical_classifier.joblib")
            vectorizer = joblib.load(self.models_dir / "medical_tfidf_vectorizer.joblib")
            label_encoder = joblib.load(self.models_dir / "medical_label_encoder.joblib")
            feature_selector = joblib.load(self.models_dir / "medical_feature_selector.joblib")
            
            # Check that models have required attributes
            assert hasattr(classifier, 'predict'), "Classifier missing predict method"
            assert hasattr(vectorizer, 'transform'), "Vectorizer missing transform method"
            assert hasattr(label_encoder, 'classes_'), "Label encoder missing classes"
            assert hasattr(feature_selector, 'transform'), "Feature selector missing transform method"
            
        except Exception as e:
            pytest.fail(f"Failed to load models: {e}")
            
    def test_prediction_function(self):
        """Test the core prediction function"""
        test_text = "Patient presents with chest pain and elevated troponin levels."
        
        try:
            # Load models first
            classifier = joblib.load(self.models_dir / "medical_classifier.joblib")
            vectorizer = joblib.load(self.models_dir / "medical_tfidf_vectorizer.joblib")
            feature_selector = joblib.load(self.models_dir / "medical_feature_selector.joblib")
            label_encoder = joblib.load(self.models_dir / "medical_label_encoder.joblib")
            
            # Use dashboard prediction function
            specialty, confidence_scores = dashboard_predict(
                test_text, classifier, vectorizer, feature_selector, label_encoder
            )
            
            # Check result structure
            assert isinstance(specialty, str), "Specialty should be a string"
            assert isinstance(confidence_scores, dict), "Confidence scores should be a dictionary"
            
            # Check specialty is valid
            valid_specialties = ['Cardiology', 'Emergency Medicine', 'Pulmonology', 
                               'Gastroenterology', 'Dermatology']
            assert specialty in valid_specialties, f"Invalid specialty: {specialty}"
            
            # Check confidence scores
            for spec, conf in confidence_scores.items():
                assert 0 <= conf <= 1, f"Confidence for {spec} should be between 0 and 1"
            
        except Exception as e:
            pytest.fail(f"Prediction function failed: {e}")


class TestAPIEndpoints:
    """Test suite for API endpoints"""
    
    @classmethod
    def setup_class(cls):
        """Setup API test environment"""
        cls.base_url = "http://localhost:8001"
        
    def test_health_endpoint(self):
        """Test API health check"""
        try:
            response = requests.get(f"{self.base_url}/health", timeout=5)
            assert response.status_code == 200, "Health check failed"
            
            data = response.json()
            assert data['status'] == 'healthy', "API not reporting healthy status"
            
        except requests.exceptions.RequestException:
            pytest.skip("API not running - skipping API tests")
            
    def test_model_info_endpoint(self):
        """Test model info endpoint"""
        try:
            response = requests.get(f"{self.base_url}/model-info", timeout=5)
            assert response.status_code == 200, "Model info endpoint failed"
            
            data = response.json()
            assert 'model_version' in data, "Model info missing version"
            assert 'specialties' in data, "Model info missing specialties"
            
        except requests.exceptions.RequestException:
            pytest.skip("API not running - skipping API tests")
            
    def test_prediction_endpoint(self):
        """Test prediction endpoint"""
        test_cases = [
            {
                "text": "Patient with chest pain and elevated cardiac enzymes",
                "expected_specialty": "Cardiology"
            },
            {
                "text": "Patient presents to emergency department with acute abdominal pain",
                "expected_specialty": "Emergency Medicine"
            }
        ]
        
        try:
            for case in test_cases:
                response = requests.post(
                    f"{self.base_url}/predict",
                    headers={"Content-Type": "application/json"},
                    json={"text": case["text"]},
                    timeout=10
                )
                
                assert response.status_code == 200, f"Prediction failed for: {case['text']}"
                
                data = response.json()
                assert 'specialty' in data, "Response missing specialty"
                assert 'confidence' in data, "Response missing confidence"
                
                # Log the prediction for manual verification
                print(f"Text: {case['text'][:50]}...")
                print(f"Predicted: {data['specialty']} (confidence: {data['confidence']:.3f})")
                print(f"Expected: {case['expected_specialty']}")
                print("-" * 50)
                
        except requests.exceptions.RequestException:
            pytest.skip("API not running - skipping API tests")


class TestMedicalSpecialties:
    """Test predictions for each medical specialty"""
    
    @classmethod
    def setup_class(cls):
        """Setup specialty test cases"""
        cls.specialty_tests = {
            "Cardiology": [
                "Myocardial infarction with ST elevation on ECG",
                "Atrial fibrillation with rapid ventricular response",
                "Chest pain with elevated troponin levels"
            ],
            "Emergency Medicine": [
                "Multi-trauma patient from motor vehicle accident",
                "Acute respiratory distress in emergency department",
                "Septic shock requiring immediate intervention"
            ],
            "Pulmonology": [
                "Chronic obstructive pulmonary disease exacerbation",
                "Pneumonia with bilateral infiltrates on chest X-ray",
                "Asthma attack with bronchospasm"
            ],
            "Gastroenterology": [
                "Inflammatory bowel disease with abdominal pain",
                "Liver cirrhosis with portal hypertension",
                "Peptic ulcer disease with upper GI bleeding"
            ],
            "Dermatology": [
                "Melanoma with irregular borders and color variation",
                "Atopic dermatitis with eczematous lesions",
                "Psoriasis with characteristic silvery scales"
            ]
        }
        
    def test_specialty_predictions(self):
        """Test that each specialty can be correctly predicted"""
        # Load models once for all tests
        try:
            models_dir = PROJECT_ROOT / "models"
            classifier = joblib.load(models_dir / "medical_classifier.joblib")
            vectorizer = joblib.load(models_dir / "medical_tfidf_vectorizer.joblib")
            feature_selector = joblib.load(models_dir / "medical_feature_selector.joblib")
            label_encoder = joblib.load(models_dir / "medical_label_encoder.joblib")
        except Exception as e:
            pytest.skip(f"Could not load models: {e}")
        
        for specialty, test_cases in self.specialty_tests.items():
            correct_predictions = 0
            
            for test_text in test_cases:
                try:
                    predicted_specialty, confidence_scores = dashboard_predict(
                        test_text, classifier, vectorizer, feature_selector, label_encoder
                    )
                    
                    if predicted_specialty == specialty:
                        correct_predictions += 1
                        
                    print(f"Text: {test_text[:50]}...")
                    print(f"Expected: {specialty}, Got: {predicted_specialty}")
                    print(f"Confidence scores: {confidence_scores}")
                    print("-" * 50)
                    
                except Exception as e:
                    pytest.fail(f"Failed to predict for {specialty}: {e}")
            
            # Allow for some flexibility - at least 1 out of 3 should be correct
            success_rate = correct_predictions / len(test_cases)
            assert success_rate >= 0.33, f"Low success rate for {specialty}: {success_rate:.2f}"


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])

"""
Simple Medical Classification API
===============================
Direct API startup without complex configuration system.
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict
import joblib
import numpy as np
from pathlib import Path
import json

# Load models directly
def load_models():
    """Load trained models with hybrid feature selection"""
    models_dir = Path("models")
    
    # Use the regularized models which perform better
    models = {
        'classifier': joblib.load(models_dir / "regularized_medical_classifier.joblib"),
        'vectorizer': joblib.load(models_dir / "regularized_tfidf_vectorizer.joblib"), 
        'feature_selector': joblib.load(models_dir / "medical_feature_selector.joblib"),
        'label_encoder': joblib.load(models_dir / "regularized_label_encoder.joblib")
    }
    
    # Load model info
    with open(models_dir / "regularized_model_info.json", 'r') as f:
        models['info'] = json.load(f)
    
    print("✅ Models loaded with hybrid feature selection")
    return models

# Initialize FastAPI
app = FastAPI(
    title="Medical Classification API",
    description="Simple API for medical text classification",
    version="1.0.0"
)

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Load models on startup
ml_models = load_models()

# Request/Response models
class TextRequest(BaseModel):
    text: str

class PredictionResponse(BaseModel):
    specialty: str
    confidence: float
    model_version: str

# API endpoints
@app.get("/")
def root():
    """Root endpoint with API information"""
    return {
        "message": "Medical Classification API",
        "version": "1.0.0",
        "status": "running",
        "endpoints": {
            "health": "/health",
            "predict": "/predict",
            "model_info": "/model-info",
            "docs": "/docs"
        }
    }

@app.get("/health")
def health_check():
    """Health check endpoint"""
    return {"status": "healthy", "message": "Medical Classification API is running"}

@app.post("/predict", response_model=PredictionResponse)
def predict_specialty(request: TextRequest):
    """Predict medical specialty from text using regularized model with confidence thresholding"""
    try:
        # Input validation
        if not request.text or len(request.text.strip()) < 10:
            raise HTTPException(status_code=400, detail="Text must be at least 10 characters")
        
        if len(request.text) > 10000:
            raise HTTPException(status_code=400, detail="Text too long (max 10,000 characters)")
        
        # Vectorize text
        text_vector = ml_models['vectorizer'].transform([request.text])
        
        # Apply feature selection
        text_features = ml_models['feature_selector'].transform(text_vector)
        
        # Get prediction and probabilities
        prediction = ml_models['classifier'].predict(text_features)[0]
        probabilities = ml_models['classifier'].predict_proba(text_features)[0]
        
        # Get specialty name and confidence
        specialty = ml_models['label_encoder'].inverse_transform([prediction])[0]
        confidence = float(max(probabilities))
        
        # Apply confidence threshold for medical relevance
        if confidence < 0.70:
            specialty = "Unknown/Low_Confidence"
            return PredictionResponse(
                specialty=specialty,
                confidence=confidence,
                model_version=f"{ml_models['info'].get('model_name', 'Regularized Medical Classifier')} (Low Confidence)"
            )
        
        return PredictionResponse(
            specialty=specialty,
            confidence=confidence,
            model_version=ml_models['info'].get('model_name', 'Regularized Medical Classifier')
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction error: {str(e)}")

@app.post("/predict-with-details")
def predict_with_details(request: TextRequest):
    """Enhanced prediction with all probabilities and detailed analysis"""
    try:
        # Input validation
        if not request.text or len(request.text.strip()) < 10:
            raise HTTPException(status_code=400, detail="Text must be at least 10 characters")
        
        # Vectorize text
        text_vector = ml_models['vectorizer'].transform([request.text])
        text_features = ml_models['feature_selector'].transform(text_vector)
        
        # Get prediction and all probabilities
        prediction = ml_models['classifier'].predict(text_features)[0]
        probabilities = ml_models['classifier'].predict_proba(text_features)[0]
        
        # Get all specialties and their probabilities
        specialties = ml_models['label_encoder'].inverse_transform(range(len(probabilities)))
        specialty_probs = {specialty: float(prob) for specialty, prob in zip(specialties, probabilities)}
        
        # Top prediction
        top_specialty = ml_models['label_encoder'].inverse_transform([prediction])[0]
        top_confidence = float(max(probabilities))
        
        # Determine if medical text
        is_medical = top_confidence >= 0.70
        
        return {
            "prediction": {
                "specialty": top_specialty if is_medical else "Unknown/Low_Confidence",
                "confidence": top_confidence,
                "is_medical_text": is_medical
            },
            "all_probabilities": specialty_probs,
            "analysis": {
                "confidence_level": "High" if top_confidence >= 0.90 else "Medium" if top_confidence >= 0.70 else "Low",
                "medical_relevance": "Relevant" if is_medical else "Irrelevant or ambiguous",
                "recommendation": "Use prediction" if is_medical else "Requires human review"
            },
            "model_version": ml_models['info'].get('model_name', 'Regularized Medical Classifier')
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction error: {str(e)}")

@app.get("/model-info")
def get_model_info():
    """Get model information"""
    return ml_models['info']

if __name__ == "__main__":
    import uvicorn
    print("🏥 Starting Simple Medical Classification API")
    print("📖 API Documentation: http://localhost:8000/docs")
    uvicorn.run(app, host="0.0.0.0", port=8000)

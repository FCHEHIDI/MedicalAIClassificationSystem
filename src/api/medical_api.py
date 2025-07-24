"""
Medical Classification FastAPI Service
=====================================

Professional production-ready API for medical text classification.
Built with healthcare industry best practices.
"""

from fastapi import FastAPI, HTTPException, Depends, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.security import HTTPBearer
from pydantic import BaseModel, Field, validator
from typing import List, Dict, Optional, Any
import joblib
import numpy as np
from pathlib import Path
import logging
from datetime import datetime
import uvicorn
import json
import os

# Configure logging
os.makedirs('logs', exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/medical_api.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Initialize FastAPI with metadata
app = FastAPI(
    title="Medical Text Classification API",
    description="Professional API for classifying medical texts into specialties using ML",
    version="1.0.0",
    contact={
        "name": "Medical AI Team",
        "email": "fareschehidi7@gmail.com"
    },
    license_info={
        "name": "MIT License",
        "url": "https://opensource.org/licenses/MIT"
    },
    docs_url="/docs",
    redoc_url="/redoc"
)

# CORS middleware for web frontend integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Security
security = HTTPBearer()

# Global model variables
model = None
vectorizer = None
feature_selector = None
label_encoder = None
model_info = None

class MedicalTextRequest(BaseModel):
    """Request model for medical text classification"""
    text: str = Field(
        ...,
        min_length=10,
        max_length=5000,
        description="Medical text to classify (10-5000 characters)"
    )
    include_confidence_scores: bool = Field(
        default=True,
        description="Include confidence scores for all specialties"
    )
    minimum_confidence: float = Field(
        default=0.1,
        ge=0.0,
        le=1.0,
        description="Minimum confidence threshold (0.0-1.0)"
    )
    
    @validator('text')
    def validate_text(cls, v):
        """Validate medical text input"""
        if not v.strip():
            raise ValueError("Text cannot be empty")
        
        # Check for potentially sensitive information (basic HIPAA compliance)
        sensitive_patterns = ['ssn', 'social security', 'patient id', 'mrn']
        text_lower = v.lower()
        for pattern in sensitive_patterns:
            if pattern in text_lower:
                raise ValueError(f"Text may contain sensitive information: {pattern}")
        
        return v.strip()

class SpecialtyPrediction(BaseModel):
    """Individual specialty prediction"""
    specialty: str
    confidence: float = Field(..., ge=0.0, le=1.0)
    risk_level: str

class MedicalTextResponse(BaseModel):
    """Response model for medical text classification"""
    predicted_specialty: str
    confidence: float = Field(..., ge=0.0, le=1.0)
    risk_assessment: str
    all_predictions: Optional[List[SpecialtyPrediction]] = None
    metadata: Dict[str, Any] = {}
    timestamp: datetime
    model_version: str
    processing_time_ms: float

class HealthResponse(BaseModel):
    """Health check response"""
    status: str
    timestamp: datetime
    model_loaded: bool
    uptime_seconds: float
    version: str

# Startup event
@app.on_event("startup")
async def load_models():
    """Load trained models on startup"""
    global model, vectorizer, feature_selector, label_encoder, model_info
    
    try:
        logger.info("Loading medical classification models...")
        
        # Create logs directory
        os.makedirs("logs", exist_ok=True)
        
        models_dir = Path("models")
        
        # Load models (5 medical specialties - our trained models)
        model = joblib.load(models_dir / "best_medical_classifier.joblib")
        vectorizer = joblib.load(models_dir / "medical_tfidf_vectorizer.joblib") 
        feature_selector = joblib.load(models_dir / "medical_feature_selector.joblib")
        label_encoder = joblib.load(models_dir / "medical_label_encoder.joblib")
        
        # Load model metadata
        with open(models_dir / "model_info.json", 'r') as f:
            model_info = json.load(f)
        
        logger.info(f"Models loaded successfully:")
        logger.info(f"  - Model: {model_info['model_name']}")
        logger.info(f"  - Test Accuracy: {model_info['accuracy']:.3f}")
        logger.info(f"  - Features: {model_info['feature_count']}")
        logger.info(f"  - Classes: {len(model_info['classes'])}")
        
    except Exception as e:
        logger.error(f"Failed to load models: {e}")
        raise RuntimeError(f"Model loading failed: {e}")

def get_risk_level(confidence: float) -> str:
    """Determine risk level based on confidence score"""
    if confidence >= 0.8:
        return "LOW"
    elif confidence >= 0.6:
        return "MODERATE"
    else:
        return "HIGH"

def validate_models_loaded():
    """Dependency to ensure models are loaded"""
    if model is None or vectorizer is None:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Models not loaded. Please check server status."
        )

@app.get("/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    return HealthResponse(
        status="healthy" if model is not None else "unhealthy",
        timestamp=datetime.now(),
        model_loaded=model is not None,
        uptime_seconds=0,  # Can implement actual uptime tracking
        version="1.0.0"
    )

@app.get("/model-info")
async def get_model_info(models_loaded: None = Depends(validate_models_loaded)):
    """Get information about the loaded model"""
    return {
        "model_info": model_info,
        "specialties": model_info["classes"],
        "training_data_source": model_info["dataset_source"],
        "regularization_applied": model_info["regularization"]
    }

@app.post("/classify", response_model=MedicalTextResponse)
async def classify_medical_text(
    request: MedicalTextRequest,
    models_loaded: None = Depends(validate_models_loaded)
):
    """
    Classify medical text into medical specialties.
    
    This endpoint uses a trained machine learning model to classify
    medical text into one of 5 specialties: Cardiology, Dermatology,
    Emergency, Gastroenterology, or Pulmonology.
    """
    start_time = datetime.now()
    
    try:
        logger.info(f"Processing classification request: {len(request.text)} characters")
        
        # Preprocess text
        processed_text = request.text.strip()
        
        # Transform text through pipeline
        text_tfidf = vectorizer.transform([processed_text])
        text_features = feature_selector.transform(text_tfidf)
        
        # Make prediction
        prediction = model.predict(text_features)[0]
        probabilities = model.predict_proba(text_features)[0]
        
        # Get specialty name
        predicted_specialty = label_encoder.inverse_transform([prediction])[0]
        max_confidence = float(np.max(probabilities))
        
        # Create all predictions with confidence scores
        all_predictions = []
        if request.include_confidence_scores:
            for i, prob in enumerate(probabilities):
                if prob >= request.minimum_confidence:
                    specialty_name = label_encoder.inverse_transform([i])[0]
                    all_predictions.append(SpecialtyPrediction(
                        specialty=specialty_name,
                        confidence=float(prob),
                        risk_level=get_risk_level(float(prob))
                    ))
            
            # Sort by confidence
            all_predictions.sort(key=lambda x: x.confidence, reverse=True)
        
        # Calculate processing time
        processing_time = (datetime.now() - start_time).total_seconds() * 1000
        
        # Risk assessment
        risk_assessment = get_risk_level(max_confidence)
        if max_confidence < 0.5:
            risk_assessment += " - Consider manual review"
        
        # Metadata
        metadata = {
            "text_length": len(processed_text),
            "feature_count": text_features.shape[1],
            "model_features": model_info["features_selected"],
            "confidence_threshold": request.minimum_confidence
        }
        
        response = MedicalTextResponse(
            predicted_specialty=predicted_specialty,
            confidence=max_confidence,
            risk_assessment=risk_assessment,
            all_predictions=all_predictions if request.include_confidence_scores else None,
            metadata=metadata,
            timestamp=datetime.now(),
            model_version=model_info["model_name"],
            processing_time_ms=processing_time
        )
        
        logger.info(f"Classification completed: {predicted_specialty} ({max_confidence:.3f})")
        return response
        
    except Exception as e:
        logger.error(f"Classification error: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Classification failed: {str(e)}"
        )

class BatchClassifyRequest(BaseModel):
    """Request model for batch classification"""
    texts: List[str] = Field(..., max_items=50, description="List of texts to classify (max 50)")
    include_confidence_scores: bool = Field(
        default=True,
        description="Include confidence scores for all specialties"
    )
    minimum_confidence: float = Field(
        default=0.1,
        ge=0.0,
        le=1.0,
        description="Minimum confidence threshold (0.0-1.0)"
    )

@app.post("/batch-classify")
async def batch_classify(
    request: BatchClassifyRequest,
    models_loaded: None = Depends(validate_models_loaded)
):
    """
    Batch classify multiple medical texts.
    Maximum 50 texts per request for performance.
    """
    start_time = datetime.now()
    
    try:
        logger.info(f"Processing batch classification: {len(request.texts)} texts")
        
        if len(request.texts) > 50:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Maximum 50 texts allowed per batch request"
            )
        
        results = []
        
        for i, text in enumerate(request.texts):
            if not text.strip():
                continue
                
            # Process each text
            text_tfidf = vectorizer.transform([text.strip()])
            text_features = feature_selector.transform(text_tfidf)
            
            prediction = model.predict(text_features)[0]
            probabilities = model.predict_proba(text_features)[0]
            
            predicted_specialty = label_encoder.inverse_transform([prediction])[0]
            confidence = float(np.max(probabilities))
            
            results.append({
                "index": i,
                "text_preview": text[:100] + "..." if len(text) > 100 else text,
                "predicted_specialty": predicted_specialty,
                "confidence": confidence,
                "risk_level": get_risk_level(confidence)
            })
        
        processing_time = (datetime.now() - start_time).total_seconds() * 1000
        
        return {
            "results": results,
            "total_processed": len(results),
            "processing_time_ms": processing_time,
            "timestamp": datetime.now(),
            "model_version": model_info["model_name"]
        }
        
    except Exception as e:
        logger.error(f"Batch classification error: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Batch classification failed: {str(e)}"
        )

@app.get("/specialties")
async def get_specialties(models_loaded: None = Depends(validate_models_loaded)):
    """Get list of available medical specialties"""
    return {
        "specialties": model_info["classes"],
        "total_count": len(model_info["classes"]),
        "descriptions": {
            "Cardiology": "Heart and cardiovascular diseases",
            "Dermatology": "Skin diseases and conditions",
            "Emergency": "Emergency medicine and trauma care",
            "Gastroenterology": "Digestive system diseases",
            "Pulmonology": "Lung and respiratory diseases"
        }
    }

@app.exception_handler(HTTPException)
async def http_exception_handler(request, exc):
    """Custom HTTP exception handler"""
    logger.warning(f"HTTP exception: {exc.status_code} - {exc.detail}")
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": True,
            "message": exc.detail,
            "status_code": exc.status_code,
            "timestamp": datetime.now().isoformat()
        }
    )

@app.exception_handler(Exception)
async def general_exception_handler(request, exc):
    """General exception handler"""
    logger.error(f"Unexpected error: {exc}")
    return JSONResponse(
        status_code=500,
        content={
            "error": True,
            "message": "Internal server error",
            "timestamp": datetime.now().isoformat()
        }
    )

if __name__ == "__main__":
    # Production configuration
    uvicorn.run(
        "medical_api:app",
        host="0.0.0.0",
        port=8000,
        reload=False,  # Set to False in production
        workers=1,     # Adjust based on your server capacity
        log_level="info",
        access_log=True
    )

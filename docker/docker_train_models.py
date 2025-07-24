#!/usr/bin/env python3
"""
Docker-compatible model training script
Ensures models are trained with the same scikit-learn version as deployment
"""
import json
import joblib
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score
import numpy as np

def train_models_docker_compatible():
    """Train models with Docker-compatible versions"""
    print("🐳 Docker-Compatible Model Training")
    print("=" * 50)
    
    # Load data - adjust path for docker subdirectory
    data_path = Path("../data/processed/complete_medical_dataset.json")
    with open(data_path, 'r') as f:
        dataset = json.load(f)
    
    # Extract documents from the dataset structure
    data = dataset.get('documents', [])
    if not data and 'data' in dataset:
        data = dataset['data']
    
    print(f"📚 Loaded {len(data)} medical documents")
    
    # Prepare data
    texts = [item['text'] for item in data]
    labels = [item['specialty'] for item in data]
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        texts, labels, test_size=0.2, random_state=42, stratify=labels
    )
    
    print(f"✅ Training: {len(X_train)} documents")
    print(f"✅ Testing: {len(X_test)} documents")
    
    # Create vectorizer with same parameters
    print("🔧 Creating TF-IDF Vectorizer...")
    vectorizer = TfidfVectorizer(
        max_features=1000,
        min_df=2,
        max_df=0.8,
        stop_words='english',
        ngram_range=(1, 2),
        lowercase=True,
        strip_accents='unicode'
    )
    
    # Fit and transform training data
    print("📊 Fitting TF-IDF on training data...")
    X_train_tfidf = vectorizer.fit_transform(X_train)
    X_test_tfidf = vectorizer.transform(X_test)
    
    print(f"📊 TF-IDF Matrix Shape: {X_train_tfidf.shape}")
    print(f"📊 Vocabulary Size: {len(vectorizer.vocabulary_)}")
    
    # Feature selection
    print("🎯 Selecting top 500 features...")
    feature_selector = SelectKBest(chi2, k=500)
    X_train_selected = feature_selector.fit_transform(X_train_tfidf, y_train)
    X_test_selected = feature_selector.transform(X_test_tfidf)
    
    print(f"📊 Selected Features Shape: {X_train_selected.shape}")
    
    # Label encoding
    print("🏷️ Encoding labels...")
    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)
    y_test_encoded = label_encoder.transform(y_test)
    
    print(f"🏷️ Classes: {label_encoder.classes_}")
    
    # Train Random Forest
    print("🌲 Training Random Forest...")
    model = RandomForestClassifier(
        n_estimators=100,
        max_depth=10,
        min_samples_split=5,
        min_samples_leaf=2,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    
    model.fit(X_train_selected, y_train_encoded)
    
    # Test the model
    print("🧪 Testing model performance...")
    y_pred = model.predict(X_test_selected)
    accuracy = accuracy_score(y_test_encoded, y_pred)
    
    print(f"🎯 Test Accuracy: {accuracy:.1%}")
    print("\n📊 Classification Report:")
    print(classification_report(y_test_encoded, y_pred, target_names=label_encoder.classes_))
    
    # Test TF-IDF vectorizer is working
    print("✅ Testing vectorizer compatibility...")
    test_text = "Patient presents with chest pain and shortness of breath"
    test_tfidf = vectorizer.transform([test_text])
    test_features = feature_selector.transform(test_tfidf)
    test_pred = model.predict(test_features)
    predicted_specialty = label_encoder.inverse_transform(test_pred)[0]
    print(f"✅ Test prediction: {predicted_specialty}")
    
    # Save models - adjust path for docker subdirectory
    models_dir = Path("../models")
    models_dir.mkdir(exist_ok=True)
    
    print("💾 Saving Docker-compatible models...")
    joblib.dump(model, models_dir / "docker_medical_classifier.joblib")
    joblib.dump(vectorizer, models_dir / "docker_tfidf_vectorizer.joblib")
    joblib.dump(feature_selector, models_dir / "docker_feature_selector.joblib")
    joblib.dump(label_encoder, models_dir / "docker_label_encoder.joblib")
    
    # Save model info
    model_info = {
        "model_name": "Docker-Compatible Random Forest",
        "test_accuracy": float(accuracy),
        "features_selected": 500,
        "training_size": len(X_train),
        "test_size": len(X_test),
        "specialties": label_encoder.classes_.tolist(),
        "sklearn_version": "Compatible with Docker deployment"
    }
    
    with open(models_dir / "docker_model_info.json", 'w') as f:
        json.dump(model_info, f, indent=2)
    
    print("✅ Docker-compatible models saved successfully!")
    print("=" * 50)
    print(f"🎯 Final Accuracy: {accuracy:.1%}")
    print("🐳 Ready for Docker deployment!")
    
    return True

if __name__ == "__main__":
    train_models_docker_compatible()

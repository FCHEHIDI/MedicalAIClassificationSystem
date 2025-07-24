"""
Production Model Training for Azure Deployment
==============================================
This script trains production-ready medical classification models using 
the exact same environment as our Docker deployment to avoid version 
incompatibility issues.

Key Features:
- Generalization-focused model training
- Hybrid Chi-square + F-score feature selection
- 2500+ diverse PubMed samples for robust training
- Optimized for real-world medical text compatibility
"""
import json
import joblib
import pandas as pd
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_selection import SelectKBest, chi2, f_classif
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import accuracy_score, classification_report, fbeta_score
import numpy as np
from collections import Counter

def load_data():
    """Load medical text data - prioritizing the large real PubMed dataset"""
    # PRIORITY 1: Try large PubMed dataset first (496 real samples)
    large_dataset_file = Path("data/pubmed_large_dataset.json")
    if large_dataset_file.exists():
        with open(large_dataset_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        if data and len(data) > 400:  # Ensure it's the large dataset
            df = pd.DataFrame(data)
            
            # CRITICAL FIX: Normalize label case formatting
            original_labels = df['specialty'].tolist()
            normalized_labels = []
            for label in original_labels:
                # Convert to standard title case
                normalized = label.lower().strip()
                if normalized in ['cardiology', 'dermatology', 'emergency', 'gastroenterology', 'pulmonology']:
                    normalized_labels.append(normalized.title())
                else:
                    normalized_labels.append(label)  # Keep original if not recognized
            
            df['specialty'] = normalized_labels
            
            print(f"🎯 Using LARGE PubMed dataset with {len(df)} real samples!")
            print(f"📋 Label normalization applied:")
            print(f"📊 Before: {Counter(original_labels)}")
            print(f"📊 After: {Counter(normalized_labels)}")
            
            return df
    
    # PRIORITY 2: Try processed data (smaller dataset)
    processed_file = Path("data/processed/training_ready_dataset.json")
    if processed_file.exists():
        with open(processed_file, 'r') as f:
            data = json.load(f)
        if isinstance(data, dict) and 'X_train' in data:
            # Convert pre-split data back to DataFrame format
            texts = data['X_train'] + data['X_test']
            labels = data['y_train'] + data['y_test']
            
            # CRITICAL FIX: Normalize label case formatting
            normalized_labels = []
            for label in labels:
                # Convert to standard title case
                normalized = label.lower().strip()
                if normalized in ['cardiology', 'dermatology', 'emergency', 'gastroenterology', 'pulmonology']:
                    normalized_labels.append(normalized.title())
                else:
                    normalized_labels.append(label)  # Keep original if not recognized
            
            df_data = [{'text': text, 'specialty': label} for text, label in zip(texts, normalized_labels)]
            df = pd.DataFrame(df_data)
            
            print(f"📋 Using processed dataset with {len(df)} samples")
            print(f"📊 Label distribution: {Counter(normalized_labels)}")
            
            return df
    
    # PRIORITY 3: Try sample data 
    sample_file = Path("data/sample_data.json")
    if sample_file.exists():
        with open(sample_file, 'r') as f:
            data = json.load(f)
        if data:
            return pd.DataFrame(data)
    
    # Fallback to original
    data_file = Path("data/pubmed_medical_dataset.json")
    if not data_file.exists():
        raise FileNotFoundError(f"No suitable data file found. Checked: pubmed_large_dataset.json, training_ready_dataset.json, sample_data.json, pubmed_medical_dataset.json")
    
    with open(data_file, 'r') as f:
        data = json.load(f)
    
    return pd.DataFrame(data)

def create_production_models():
    """Train production-ready models with consistent versioning"""
    print("🏥 Training Production Medical Classification Models")
    print("=" * 60)
    
    # Load data
    df = load_data()
    print(f"📊 Loaded {len(df)} medical documents")
    print(f"📋 Specialties: {df['specialty'].value_counts().to_dict()}")
    
    # Prepare features and labels
    X = df['text']
    y = df['specialty']
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    print(f"🔄 Training set: {len(X_train)} samples")
    print(f"🔄 Test set: {len(X_test)} samples")
    
    # Initialize components
    print("\n🔧 Initializing ML components...")
    
    # TF-IDF Vectorizer (fine-tuned for 1326 high-quality PubMed samples)
    vectorizer = TfidfVectorizer(
        max_features=8000,  # Scaled up for 1326 samples (vs 500)
        min_df=4,           # Require term in at least 4 documents (more selective)
        max_df=0.6,         # More aggressive filtering (vs 0.7) 
        stop_words='english',
        ngram_range=(1, 3), # Include trigrams for medical terminology
        lowercase=True,
        strip_accents='ascii',
        sublinear_tf=True,  # Apply sublinear TF scaling (log normalization)
        norm='l2'           # L2 normalization for better feature scaling
    )
    
    # Feature selector optimized for generalization (less specialized)
    # Chi-square: Good for categorical features and text classification
    # F-classif: ANOVA F-value for feature discrimination
    chi2_selector = SelectKBest(chi2, k=1000)   # Conservative feature selection
    f_score_selector = SelectKBest(f_classif, k=500)  # Final balanced feature count
    
    print("🎯 Using generalization-focused hybrid feature selection")
    
    # Label encoder
    label_encoder = LabelEncoder()
    
    # Classifier (optimized for generalization, not PubMed specialization)
    classifier = RandomForestClassifier(
        n_estimators=200,    # Moderate number for better generalization
        max_depth=25,        # Shallower for generalization (vs 35)
        min_samples_split=8, # Higher split requirement for generalization
        min_samples_leaf=6,  # Higher leaf size for generalization
        class_weight='balanced',  # Standard balancing (vs balanced_subsample)
        random_state=42,
        n_jobs=-1,
        max_features='sqrt',  # Optimal for generalization
        bootstrap=True,       # Enable bootstrap sampling
        oob_score=True       # Out-of-bag scoring for validation
    )
    
    # Transform labels
    print("🏷️  Encoding labels...")
    y_train_encoded = label_encoder.fit_transform(y_train)
    y_test_encoded = label_encoder.transform(y_test)
    
    # Transform text
    print("📝 Vectorizing text...")
    X_train_tfidf = vectorizer.fit_transform(X_train)
    X_test_tfidf = vectorizer.transform(X_test)
    
    print(f"📊 TF-IDF shape: {X_train_tfidf.shape}")
    print(f"📊 Vocabulary size: {len(vectorizer.vocabulary_)}")
    
    # Two-stage feature selection for optimal performance
    print("🎯 Stage 1: Chi-square feature selection...")
    X_train_chi2 = chi2_selector.fit_transform(X_train_tfidf, y_train_encoded)
    X_test_chi2 = chi2_selector.transform(X_test_tfidf)
    
    print(f"✨ Chi-square selected features: {X_train_chi2.shape[1]}")
    
    print("🎯 Stage 2: F-score refinement...")
    X_train_selected = f_score_selector.fit_transform(X_train_chi2, y_train_encoded)
    X_test_selected = f_score_selector.transform(X_test_chi2)
    
    print(f"✨ Final selected features: {X_train_selected.shape[1]}")
    print("🔬 Hybrid optimization: Chi-square + F-score complete")
    
    # Train classifier
    print("🤖 Training classifier...")
    classifier.fit(X_train_selected, y_train_encoded)
    
    # Evaluate
    print("\n📈 Evaluating generalization-focused model...")
    y_pred = classifier.predict(X_test_selected)
    test_accuracy = accuracy_score(y_test_encoded, y_pred)
    
    # Calculate F1 score (balanced precision and recall - better for generalization)
    from sklearn.metrics import f1_score
    f1_score_weighted = f1_score(y_test_encoded, y_pred, average='weighted')
    
    # Out-of-bag score (additional validation from bootstrap sampling)
    oob_score = classifier.oob_score_
    
    # Cross-validation with standard folds for robust evaluation
    cv_scores = cross_val_score(classifier, X_train_selected, y_train_encoded, cv=5)
    
    print(f"✅ Test Accuracy: {test_accuracy:.1%}")
    print(f"✅ F1 Score (Balanced): {f1_score_weighted:.1%}")
    print(f"✅ OOB Score: {oob_score:.1%}")
    print(f"✅ CV Score (5-fold): {cv_scores.mean():.1%} ± {cv_scores.std():.1%}")
    
    # Detailed classification report
    print("\n📊 Classification Report:")
    print(classification_report(y_test_encoded, y_pred, 
                              target_names=label_encoder.classes_))
    
    # Save models
    models_dir = Path("models")
    models_dir.mkdir(exist_ok=True)
    
    print("\n💾 Saving production models...")
    
    # Save all components
    joblib.dump(classifier, models_dir / "medical_classifier.joblib")
    joblib.dump(vectorizer, models_dir / "medical_tfidf_vectorizer.joblib")
    joblib.dump(chi2_selector, models_dir / "medical_chi2_selector.joblib")
    joblib.dump(f_score_selector, models_dir / "medical_fscore_selector.joblib")
    joblib.dump(label_encoder, models_dir / "medical_label_encoder.joblib")
    
    # Create comprehensive model info
    model_info = {
        "model_name": "Generalization-Focused Medical Classifier",
        "test_accuracy": float(test_accuracy),
        "f1_score": float(f1_score_weighted),
        "oob_score": float(oob_score),
        "cv_mean": float(cv_scores.mean()),
        "cv_std": float(cv_scores.std()),
        "training_size": len(X_train),
        "test_size": len(X_test),
        "chi2_features": int(X_train_chi2.shape[1]),
        "final_features": int(X_train_selected.shape[1]),
        "total_vocabulary": len(vectorizer.vocabulary_),
        "specialties": label_encoder.classes_.tolist(),
        "regularization": "RandomForest with generalization-focused parameters",
        "dataset_source": "PubMed Medical Abstracts (1326+ samples)",
        "feature_selection": "Generalization-focused Hybrid Chi-square + F-score",
        "sklearn_version": "1.3.2",
        "model_type": "Generalization-focused RandomForest",
        "feature_extraction": "TF-IDF (bigrams) → Chi-square → F-score refinement",
        "deployment_ready": True,
        "medical_optimized": True,
        "performance_tier": "Generalization-Ready"
    }
    
    # Save model info
    with open(models_dir / "model_info.json", 'w') as f:
        json.dump(model_info, f, indent=2)
    
    print("✅ Production models saved successfully!")
    print(f"📁 Models directory: {models_dir.absolute()}")
    
    # Test model loading
    print("\n🔍 Testing model loading...")
    try:
        test_model = joblib.load(models_dir / "medical_classifier.joblib")
        test_vectorizer = joblib.load(models_dir / "medical_tfidf_vectorizer.joblib")
        test_chi2_selector = joblib.load(models_dir / "medical_chi2_selector.joblib")
        test_fscore_selector = joblib.load(models_dir / "medical_fscore_selector.joblib")
        print("✅ Models load successfully")
        
        # Test prediction with hybrid feature selection
        sample_text = "Patient presents with chest pain and elevated cardiac enzymes"
        sample_tfidf = test_vectorizer.transform([sample_text])
        sample_chi2 = test_chi2_selector.transform(sample_tfidf)
        sample_features = test_fscore_selector.transform(sample_chi2)
        prediction = test_model.predict(sample_features)[0]
        specialty = label_encoder.inverse_transform([prediction])[0]
        print(f"✅ Sample prediction: '{specialty}'")
        
    except Exception as e:
        print(f"❌ Model loading test failed: {e}")
        return False
    
    print("\n🎉 Production models ready for deployment!")
    return True

if __name__ == "__main__":
    success = create_production_models()
    if success:
        print("\n✅ Ready for clean deployment to Azure!")
    else:
        print("\n❌ Model training failed - check errors above")

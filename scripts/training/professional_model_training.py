"""
Professional Medical Model Training with Regularization
======================================================

Train robust medical classification models with proper regularization 
to prevent overfitting on the current 496-document dataset.
"""

import json
import numpy as np
from pathlib import Path
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

print("🏥 Professional Medical Model Training (Regularized)")
print("=" * 60)

# Load dataset
dataset_file = Path("data/pubmed_large_dataset.json")
with open(dataset_file, 'r', encoding='utf-8') as f:
    documents = json.load(f)

texts = [doc['text'] for doc in documents]
labels = [doc['specialty'] for doc in documents]

print(f"📚 Dataset: {len(documents)} documents")

# Professional train/validation/test split
print(f"🎯 Professional Data Splitting...")
from sklearn.model_selection import train_test_split

# First split: 70% train, 30% temp
X_train, X_temp, y_train, y_temp = train_test_split(
    texts, labels, test_size=0.3, random_state=42, stratify=labels
)

# Second split: 15% validation, 15% test
X_val, X_test, y_val, y_test = train_test_split(
    X_temp, y_temp, test_size=0.5, random_state=42, stratify=y_temp
)

print(f"✅ Training: {len(X_train)} documents")
print(f"✅ Validation: {len(X_val)} documents") 
print(f"✅ Test: {len(X_test)} documents")

# REGULARIZED Feature Extraction
print(f"\n🔧 Regularized Feature Engineering...")
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import SelectKBest, chi2

# Reduced feature set with better parameters
vectorizer = TfidfVectorizer(
    max_features=1000,      # REDUCED from 5000 to 1000
    ngram_range=(1, 2),     # Only unigrams and bigrams
    min_df=3,              # Must appear in at least 3 docs
    max_df=0.7,            # Ignore if in >70% of docs
    stop_words='english',
    lowercase=True,
    strip_accents='ascii',
    sublinear_tf=True,     # Log scaling
    norm='l2'              # L2 normalization
)

# Extract features
X_train_tfidf = vectorizer.fit_transform(X_train)
X_val_tfidf = vectorizer.transform(X_val)
X_test_tfidf = vectorizer.transform(X_test)

# Encode labels
label_encoder = LabelEncoder()
y_train_encoded = label_encoder.fit_transform(y_train)
y_val_encoded = label_encoder.transform(y_val)
y_test_encoded = label_encoder.transform(y_test)

print(f"✅ Feature matrix: {X_train_tfidf.shape}")
print(f"✅ Feature/Sample ratio: {X_train_tfidf.shape[1]/X_train_tfidf.shape[0]:.2f}")

# Further feature selection
print(f"🎯 Feature Selection with Chi-Square...")
feature_selector = SelectKBest(chi2, k=500)  # Top 500 most relevant features
X_train_selected = feature_selector.fit_transform(X_train_tfidf, y_train_encoded)
X_val_selected = feature_selector.transform(X_val_tfidf)
X_test_selected = feature_selector.transform(X_test_tfidf)

print(f"✅ Selected features: {X_train_selected.shape[1]}")
print(f"✅ New ratio: {X_train_selected.shape[1]/X_train_selected.shape[0]:.2f}")

# REGULARIZED Models
print(f"\n🤖 Training Regularized Models...")
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import classification_report, accuracy_score

models = {
    'Regularized Random Forest': RandomForestClassifier(
        n_estimators=50,           # Reduced from 100
        max_depth=10,             # Limited depth
        min_samples_split=10,     # Higher split threshold
        min_samples_leaf=5,       # Higher leaf threshold
        max_features='sqrt',      # Feature subsampling
        random_state=42,
        class_weight='balanced'   # Handle class imbalance
    ),
    'L2 Logistic Regression': LogisticRegression(
        C=0.1,                   # Strong regularization
        max_iter=1000,
        random_state=42,
        class_weight='balanced',
        solver='liblinear'
    ),
    'RBF SVM (Regularized)': SVC(
        C=1.0,                   # Moderate regularization
        gamma='scale',           # Auto gamma
        kernel='rbf',
        random_state=42,
        probability=True,
        class_weight='balanced'
    )
}

results = {}

for name, model in models.items():
    print(f"\n🔬 Training {name}...")
    
    # Train
    model.fit(X_train_selected, y_train_encoded)
    
    # Predictions
    train_pred = model.predict(X_train_selected)
    val_pred = model.predict(X_val_selected)
    test_pred = model.predict(X_test_selected)
    
    # Accuracies
    train_acc = accuracy_score(y_train_encoded, train_pred)
    val_acc = accuracy_score(y_val_encoded, val_pred)
    test_acc = accuracy_score(y_test_encoded, test_pred)
    
    # Overfitting gap
    gap = train_acc - val_acc
    
    print(f"   📊 Train Accuracy: {train_acc:.3f}")
    print(f"   📊 Val Accuracy: {val_acc:.3f}")
    print(f"   📊 Test Accuracy: {test_acc:.3f}")
    print(f"   📈 Overfitting Gap: {gap:.3f}")
    
    if gap < 0.05:
        print(f"   ✅ LOW OVERFITTING")
    elif gap < 0.10:
        print(f"   ⚠️ MODERATE OVERFITTING")
    else:
        print(f"   🚨 HIGH OVERFITTING")
    
    results[name] = {
        'model': model,
        'train_acc': train_acc,
        'val_acc': val_acc,
        'test_acc': test_acc,
        'gap': gap
    }

# Best model selection based on validation accuracy
best_model_name = max(results.keys(), key=lambda x: results[x]['val_acc'])
best_model_info = results[best_model_name]

print(f"\n🏆 Best Model: {best_model_name}")
print(f"🎯 Validation Accuracy: {best_model_info['val_acc']:.3f}")
print(f"🎯 Test Accuracy: {best_model_info['test_acc']:.3f}")
print(f"📈 Overfitting Gap: {best_model_info['gap']:.3f}")

# Detailed evaluation
print(f"\n📊 Detailed Classification Report:")
test_pred = best_model_info['model'].predict(X_test_selected)
report = classification_report(
    y_test_encoded, test_pred,
    target_names=label_encoder.classes_,
    digits=3
)
print(report)

# Cross-validation on best model
print(f"\n🔄 Cross-Validation on Best Model...")
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.pipeline import Pipeline

# Create full pipeline for CV
cv_pipeline = Pipeline([
    ('tfidf', vectorizer),
    ('selector', feature_selector),
    ('classifier', best_model_info['model'])
])

cv_scores = cross_val_score(
    cv_pipeline, texts, labels, 
    cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
    scoring='accuracy'
)

print(f"   Mean CV Accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
print(f"   CV Standard Deviation: {cv_scores.std():.3f}")

# Save regularized models
print(f"\n💾 Saving Regularized Models...")
import joblib

model_dir = Path("models")
joblib.dump(best_model_info['model'], model_dir / 'regularized_medical_classifier.joblib')
joblib.dump(vectorizer, model_dir / 'regularized_tfidf_vectorizer.joblib')
joblib.dump(feature_selector, model_dir / 'medical_feature_selector.joblib')
joblib.dump(label_encoder, model_dir / 'regularized_label_encoder.joblib')

# Model metadata
model_info = {
    'model_name': best_model_name,
    'train_accuracy': float(best_model_info['train_acc']),
    'validation_accuracy': float(best_model_info['val_acc']),
    'test_accuracy': float(best_model_info['test_acc']),
    'overfitting_gap': float(best_model_info['gap']),
    'cv_mean': float(cv_scores.mean()),
    'cv_std': float(cv_scores.std()),
    'features_selected': int(X_train_selected.shape[1]),
    'feature_sample_ratio': float(X_train_selected.shape[1]/X_train_selected.shape[0]),
    'classes': list(label_encoder.classes_),
    'regularization': 'Strong L2 regularization + feature selection',
    'training_size': len(X_train),
    'dataset_source': 'PubMed NCBI Medical Literature'
}

with open(model_dir / 'regularized_model_info.json', 'w') as f:
    json.dump(model_info, f, indent=2)

print(f"✅ Regularized models saved")

# Final assessment
print(f"\n" + "="*60)
print(f"🎯 PROFESSIONAL MEDICAL MODEL ASSESSMENT")
print(f"="*60)

print(f"📊 Model Performance:")
print(f"   Best Model: {best_model_name}")
print(f"   Test Accuracy: {best_model_info['test_acc']:.1%}")
print(f"   Overfitting Gap: {best_model_info['gap']:.3f}")
print(f"   CV Stability: ±{cv_scores.std():.3f}")

print(f"\n🔧 Regularization Applied:")
print(f"   ✅ Reduced features: 5000 → 500")
print(f"   ✅ Feature/Sample ratio: 10.08 → {X_train_selected.shape[1]/X_train_selected.shape[0]:.2f}")
print(f"   ✅ Strong L2 regularization")
print(f"   ✅ Feature selection (Chi-Square)")
print(f"   ✅ Balanced class weights")
print(f"   ✅ Conservative model hyperparameters")

# Honest assessment
if best_model_info['gap'] < 0.05 and cv_scores.std() < 0.05:
    assessment = "✅ PRODUCTION READY - Well regularized"
elif best_model_info['gap'] < 0.10:
    assessment = "⚠️ DEVELOPMENT READY - Monitor performance"  
else:
    assessment = "🚨 PROOF-OF-CONCEPT - Expand dataset"

print(f"\n🏥 FINAL ASSESSMENT: {assessment}")

print(f"\n💡 PROFESSIONAL RECOMMENDATION:")
if best_model_info['gap'] < 0.10:
    print(f"   ✅ Current regularized model is suitable for development")
    print(f"   📊 Can be used for FastAPI deployment with monitoring")
    print(f"   🔍 Collect real-world data to validate performance")
    print(f"   📈 Expand dataset when possible for even better robustness")
else:
    print(f"   📈 Expand dataset to 1000+ samples per class")
    print(f"   🌍 Add diverse data sources beyond PubMed")
    print(f"   🔄 Continue with current model for proof-of-concept")

print(f"\n🚀 NEXT STEPS:")
print(f"   1. Deploy regularized model to FastAPI")
print(f"   2. Build Streamlit dashboard with confidence thresholds")
print(f"   3. Monitor real-world performance")
print(f"   4. Collect additional data iteratively")
print(f"   5. Retrain with expanded datasets")

print(f"\n" + "="*60)

"""
Medical Model Validation & Overfitting Detection
===============================================

Comprehensive validation to detect overfitting and assess model robustness.
"""

import json
import numpy as np
from pathlib import Path
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

print("🔍 Medical Model Validation & Overfitting Detection")
print("=" * 60)

# Load dataset and model
print("📚 Loading dataset and trained models...")

dataset_file = Path("data/pubmed_large_dataset.json")
with open(dataset_file, 'r', encoding='utf-8') as f:
    documents = json.load(f)

# Load trained models
import joblib
model = joblib.load("models/best_medical_classifier.joblib")
vectorizer = joblib.load("models/medical_tfidf_vectorizer.joblib")
label_encoder = joblib.load("models/medical_label_encoder.joblib")

texts = [doc['text'] for doc in documents]
labels = [doc['specialty'] for doc in documents]

print(f"✅ Dataset: {len(documents)} documents")
print(f"✅ Models loaded successfully")

# 1. K-FOLD CROSS VALIDATION
print(f"\n🔄 K-Fold Cross Validation (k=5)...")
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier

# Create pipeline
pipeline = Pipeline([
    ('tfidf', vectorizer),
    ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
])

# 5-fold stratified cross-validation
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = cross_val_score(pipeline, texts, labels, cv=skf, scoring='accuracy')

print(f"📊 Cross-Validation Results:")
print(f"   Mean Accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
print(f"   Individual Folds: {[f'{score:.3f}' for score in cv_scores]}")
print(f"   Std Deviation: {cv_scores.std():.3f}")

# High std deviation indicates overfitting
if cv_scores.std() > 0.05:
    print(f"   ⚠️  HIGH VARIANCE - Possible overfitting detected!")
else:
    print(f"   ✅ LOW VARIANCE - Model appears stable")

# 2. LEARNING CURVES
print(f"\n📈 Learning Curves Analysis...")
from sklearn.model_selection import learning_curve

# Generate learning curves
train_sizes = np.linspace(0.1, 1.0, 10)
train_sizes_abs, train_scores, val_scores = learning_curve(
    pipeline, texts, labels, 
    train_sizes=train_sizes,
    cv=5,
    scoring='accuracy',
    random_state=42
)

train_mean = train_scores.mean(axis=1)
train_std = train_scores.std(axis=1)
val_mean = val_scores.mean(axis=1)
val_std = val_scores.std(axis=1)

print(f"📊 Learning Curve Analysis:")
print(f"   Final Training Accuracy: {train_mean[-1]:.3f} ± {train_std[-1]:.3f}")
print(f"   Final Validation Accuracy: {val_mean[-1]:.3f} ± {val_std[-1]:.3f}")
print(f"   Gap (Train - Val): {train_mean[-1] - val_mean[-1]:.3f}")

# Large gap indicates overfitting
gap = train_mean[-1] - val_mean[-1]
if gap > 0.10:
    print(f"   🚨 LARGE GAP - Strong overfitting detected!")
elif gap > 0.05:
    print(f"   ⚠️  MODERATE GAP - Some overfitting possible")
else:
    print(f"   ✅ SMALL GAP - Good generalization")

# 3. DATASET SIZE ANALYSIS
print(f"\n📊 Dataset Size Analysis...")

# Current per-class distribution
label_counts = Counter(labels)
min_samples = min(label_counts.values())
max_samples = max(label_counts.values())

print(f"   Samples per class:")
for specialty, count in sorted(label_counts.items()):
    print(f"     {specialty}: {count}")

print(f"   Minimum samples: {min_samples}")
print(f"   Maximum samples: {max_samples}")
print(f"   Balance ratio: {min_samples/max_samples:.2f}")

# Industry benchmarks
print(f"\n🎯 Industry Benchmarks:")
print(f"   Minimum viable: 200-500 samples/class")
print(f"   Production ready: 1000+ samples/class") 
print(f"   Your dataset: ~{min_samples} samples/class")

if min_samples < 200:
    print(f"   📊 STATUS: SMALL - Suitable for proof-of-concept")
    recommendation = "EXPAND"
elif min_samples < 500:
    print(f"   📊 STATUS: MODERATE - Good for development")
    recommendation = "EXPAND"
else:
    print(f"   📊 STATUS: LARGE - Production ready")
    recommendation = "SUFFICIENT"

# 4. FEATURE/SAMPLE RATIO
print(f"\n🔢 Feature-to-Sample Ratio Analysis...")
n_features = vectorizer.get_feature_names_out().shape[0] if hasattr(vectorizer, 'get_feature_names_out') else 5000
n_samples = len(texts)

feature_sample_ratio = n_features / n_samples
print(f"   Features: {n_features}")
print(f"   Samples: {n_samples}")
print(f"   Ratio: {feature_sample_ratio:.2f}")

if feature_sample_ratio > 1.0:
    print(f"   🚨 HIGH-DIMENSIONAL - Overfitting risk very high!")
elif feature_sample_ratio > 0.5:
    print(f"   ⚠️  MODERATE-DIMENSIONAL - Some overfitting risk")
else:
    print(f"   ✅ LOW-DIMENSIONAL - Good ratio")

# 5. VALIDATION ON UNSEEN DATA
print(f"\n🧪 Generating Synthetic Test Cases...")

# Create challenging test cases
synthetic_tests = [
    # Cardiology
    ("Acute myocardial infarction with ST-segment elevation requiring immediate percutaneous coronary intervention", "Cardiology"),
    
    # Emergency  
    ("Multiple trauma victim with hemothorax and pelvic fracture requiring emergent surgical intervention", "Emergency"),
    
    # Pulmonology
    ("Severe acute exacerbation of chronic obstructive pulmonary disease with respiratory failure", "Pulmonology"),
    
    # Gastroenterology
    ("Severe ulcerative colitis with toxic megacolon requiring emergency colectomy", "Gastroenterology"),
    
    # Dermatology
    ("Malignant melanoma with irregular borders and asymmetric pigmentation requiring wide excision", "Dermatology"),
    
    # Edge cases - medical terms but different context
    ("Patient discusses heart problems during dermatology consultation for skin rash", "Dermatology"),
    ("Emergency room visit for cardiac symptoms but primary issue was anxiety attack", "Emergency")
]

print(f"🔮 Testing on synthetic medical cases...")
correct_predictions = 0
total_predictions = len(synthetic_tests)

for i, (text, true_label) in enumerate(synthetic_tests):
    # Predict
    text_vector = vectorizer.transform([text])
    prediction = model.predict(text_vector)[0]
    probabilities = model.predict_proba(text_vector)[0]
    predicted_label = label_encoder.inverse_transform([prediction])[0]
    confidence = max(probabilities)
    
    is_correct = predicted_label == true_label
    if is_correct:
        correct_predictions += 1
    
    status = "✅" if is_correct else "❌"
    print(f"   {i+1}. {status} Predicted: {predicted_label} (confidence: {confidence:.2f})")
    print(f"      Expected: {true_label}")
    if not is_correct:
        print(f"      Text: {text[:80]}...")
    print()

synthetic_accuracy = correct_predictions / total_predictions
print(f"🎯 Synthetic Test Accuracy: {synthetic_accuracy:.2%}")

# 6. FINAL ASSESSMENT
print(f"\n" + "="*60)
print(f"🏥 MEDICAL MODEL VALIDATION SUMMARY")
print(f"="*60)

print(f"📊 Dataset Metrics:")
print(f"   Total Documents: {len(documents)}")
print(f"   Samples per Class: {min_samples}")
print(f"   Feature/Sample Ratio: {feature_sample_ratio:.2f}")

print(f"\n🔍 Overfitting Assessment:")
print(f"   Cross-Val Accuracy: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
print(f"   Train-Val Gap: {gap:.3f}")
print(f"   Synthetic Test Accuracy: {synthetic_accuracy:.2%}")

# Overall assessment
print(f"\n🎯 OVERALL ASSESSMENT:")

risk_factors = []
if cv_scores.std() > 0.05:
    risk_factors.append("High CV variance")
if gap > 0.05:
    risk_factors.append("Large train-val gap")
if feature_sample_ratio > 0.5:
    risk_factors.append("High feature/sample ratio")
if min_samples < 200:
    risk_factors.append("Small dataset size")
if synthetic_accuracy < 0.8:
    risk_factors.append("Poor synthetic test performance")

if len(risk_factors) == 0:
    assessment = "✅ LOW OVERFITTING RISK - Model appears robust"
elif len(risk_factors) <= 2:
    assessment = "⚠️  MODERATE OVERFITTING RISK - Monitor closely"
else:
    assessment = "🚨 HIGH OVERFITTING RISK - Expand dataset recommended"

print(f"   {assessment}")

if risk_factors:
    print(f"\n   Risk Factors:")
    for factor in risk_factors:
        print(f"     • {factor}")

print(f"\n🚀 RECOMMENDATIONS:")
if recommendation == "EXPAND":
    print(f"   1. 📈 EXPAND dataset to 1000+ samples per class")
    print(f"   2. 🌍 Add diverse data sources beyond PubMed")
    print(f"   3. 🔄 Implement regularization techniques")
    print(f"   4. 📊 Monitor performance on real clinical data")
    print(f"   5. 🎯 Consider ensemble methods for robustness")
else:
    print(f"   ✅ Current dataset size appears adequate")
    print(f"   📊 Continue monitoring with real-world data")

print(f"\n💡 HONEST PROFESSIONAL ANSWER:")
print(f"   Your 496 documents are GOOD for proof-of-concept")
print(f"   but INSUFFICIENT for high-confidence production use.")
print(f"   Expand to 1000+ samples per class for robustness.")

print(f"\n" + "="*60)

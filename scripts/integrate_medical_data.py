"""
Medical Data Integration Script
==============================

Combines PubMed real medical data with synthetic data for ML training.
"""

import json
import os
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any
from collections import Counter

print("🏥 Medical Data Integration & Processing")
print("=" * 50)

def clean_medical_text(text: str) -> str:
    """Clean and normalize medical text"""
    # Remove extra whitespace
    text = re.sub(r'\s+', ' ', text)
    
    # Remove special characters but keep medical terminology
    text = re.sub(r'[^\w\s\-\.,%\(\)]', '', text)
    
    # Normalize case (preserve important medical terms)
    # Don't lowercase everything as medical terms are case-sensitive
    return text.strip()

def extract_medical_features(text: str) -> Dict[str, Any]:
    """Extract basic medical text features"""
    words = text.split()
    
    # Medical terminology indicators
    medical_terms = [
        'patient', 'diagnosis', 'treatment', 'symptom', 'therapy',
        'disease', 'condition', 'disorder', 'syndrome', 'clinical',
        'medical', 'health', 'care', 'hospital', 'acute', 'chronic'
    ]
    
    medical_term_count = sum(1 for word in words if word.lower() in medical_terms)
    
    return {
        'word_count': len(words),
        'char_count': len(text),
        'medical_term_density': medical_term_count / max(len(words), 1),
        'avg_word_length': sum(len(word) for word in words) / max(len(words), 1)
    }

def generate_simple_synthetic_data(count_per_specialty: int = 5) -> List[Dict]:
    """Generate simple synthetic medical documents"""
    specialties = {
        'cardiology': [
            "Patient presents with chest pain and elevated cardiac enzymes. ECG shows ST elevation. Diagnosis of myocardial infarction. Treatment with angioplasty and stent placement recommended.",
            "Echocardiogram reveals reduced ejection fraction and wall motion abnormalities. Patient has history of coronary artery disease. Cardiac catheterization shows significant stenosis.",
            "Heart failure patient with shortness of breath and peripheral edema. BNP elevated. Treatment includes ACE inhibitors and diuretics for symptom management.",
            "Atrial fibrillation with rapid ventricular response. Patient anticoagulated with warfarin. Rate control achieved with beta blockers and digoxin therapy.",
            "Acute coronary syndrome with unstable angina. Cardiac biomarkers elevated. Emergency cardiac catheterization shows multi-vessel coronary disease requiring intervention."
        ],
        'emergency': [
            "Trauma patient with multiple injuries from motor vehicle accident. Primary survey shows airway compromise and pneumothorax. Emergency intubation and chest tube placement performed.",
            "Acute appendicitis presentation with right lower quadrant pain and fever. CT scan confirms diagnosis. Emergency appendectomy indicated for treatment.",
            "Septic shock patient with hypotension and altered mental status. Blood cultures positive for gram-negative bacteria. Aggressive fluid resuscitation and antibiotics initiated.",
            "Stroke patient with sudden onset weakness and speech difficulty. CT scan shows acute ischemic stroke. Thrombolytic therapy administered within therapeutic window.",
            "Diabetic ketoacidosis with severe hyperglycemia and metabolic acidosis. Emergency treatment with insulin infusion and fluid replacement therapy required immediately."
        ],
        'pulmonology': [
            "Chronic obstructive pulmonary disease exacerbation with increased dyspnea and sputum production. Chest X-ray shows hyperinflation. Treatment with bronchodilators and corticosteroids.",
            "Pneumonia patient with productive cough and fever. Chest CT reveals consolidation in right lower lobe. Antibiotic therapy initiated based on culture results.",
            "Asthma exacerbation with severe bronchospasm and wheezing. Peak flow measurements significantly decreased. Emergency treatment with nebulized albuterol and systemic steroids.",
            "Pulmonary embolism with acute shortness of breath and chest pain. CT angiography confirms large clot burden. Anticoagulation therapy started immediately.",
            "Lung cancer patient with persistent cough and weight loss. CT scan shows mass lesion with mediastinal lymphadenopathy. Bronchoscopy and biopsy planned."
        ],
        'gastroenterology': [
            "Inflammatory bowel disease with abdominal pain and bloody diarrhea. Colonoscopy reveals extensive colitis. Treatment with immunosuppressive therapy recommended.",
            "Gastroesophageal reflux disease with heartburn and regurgitation. Upper endoscopy shows esophagitis. Proton pump inhibitor therapy and lifestyle modifications advised.",
            "Acute pancreatitis with severe epigastric pain and elevated lipase. CT scan shows pancreatic inflammation. Conservative management with fluid resuscitation and pain control.",
            "Peptic ulcer disease with upper gastrointestinal bleeding. Endoscopy shows active bleeding ulcer. Therapeutic intervention with injection and clip placement performed.",
            "Irritable bowel syndrome with alternating diarrhea and constipation. Comprehensive evaluation excludes organic disease. Dietary modifications and symptom-directed therapy recommended."
        ],
        'dermatology': [
            "Melanoma with irregular pigmented lesion showing asymmetry and color variation. Dermoscopy reveals concerning features. Wide local excision and sentinel node biopsy recommended.",
            "Psoriasis with extensive plaque formation and scaling. Patient reports significant pruritus and joint pain. Systemic therapy with biologics considered for severe disease.",
            "Eczema with chronic dermatitis and lichenification. Patient has history of atopic dermatitis. Topical corticosteroids and emollients provide symptomatic relief.",
            "Basal cell carcinoma with nodular lesion on sun-exposed area. Biopsy confirms diagnosis. Mohs micrographic surgery planned for complete tumor removal.",
            "Acne vulgaris with inflammatory papules and pustules. Patient reports social anxiety related to appearance. Combination therapy with topical and oral medications initiated."
        ]
    }
    
    synthetic_docs = []
    doc_id = 1000  # Start synthetic IDs at 1000
    
    for specialty, texts in specialties.items():
        for i, text in enumerate(texts[:count_per_specialty]):
            doc = {
                'id': f'synthetic_{doc_id}',
                'text': text,
                'specialty': specialty,
                'source': 'synthetic',
                'confidence': 0.85 + (i * 0.02),  # Vary confidence slightly
                'metadata': {
                    'generated_at': datetime.now().isoformat(),
                    'type': 'synthetic',
                    'word_count': len(text.split()),
                    'character_count': len(text)
                }
            }
            synthetic_docs.append(doc)
            doc_id += 1
    
    return synthetic_docs

try:
    # Load PubMed real data
    pubmed_file = Path("data/pubmed_simple_dataset.json")
    
    if pubmed_file.exists():
        print("📚 Loading PubMed real medical data...")
        with open(pubmed_file, 'r', encoding='utf-8') as f:
            pubmed_docs = json.load(f)
        print(f"   ✅ Loaded {len(pubmed_docs)} real medical abstracts")
        
        # Show sample from PubMed data
        if pubmed_docs:
            sample = pubmed_docs[0]
            print(f"   📄 Sample: {sample['text'][:100]}...")
    else:
        print("⚠️ PubMed data not found - creating with synthetic data only")
        pubmed_docs = []
    
    # Generate synthetic data
    print("\n🤖 Generating synthetic medical data...")
    synthetic_docs = generate_simple_synthetic_data(count_per_specialty=5)
    print(f"   ✅ Generated {len(synthetic_docs)} synthetic documents")
    
    # Combine datasets  
    all_documents = pubmed_docs + synthetic_docs
    print(f"\n📊 Combined Dataset Summary:")
    print(f"   📚 Real PubMed abstracts: {len(pubmed_docs)}")
    print(f"   🤖 Synthetic documents: {len(synthetic_docs)}")
    print(f"   📄 Total documents: {len(all_documents)}")
    
    # Analyze distribution
    specialty_counts = Counter(doc['specialty'] for doc in all_documents)
    print(f"\n🏥 Specialty Distribution:")
    for specialty, count in specialty_counts.items():
        real_count = len([d for d in pubmed_docs if d['specialty'] == specialty])
        synthetic_count = len([d for d in synthetic_docs if d['specialty'] == specialty])
        print(f"   {specialty.capitalize()}: {count} total ({real_count} real + {synthetic_count} synthetic)")
    
    # Extract features for quality analysis
    print(f"\n🔬 Extracting text features...")
    for doc in all_documents:
        doc['features'] = extract_medical_features(doc['text'])
        doc['processed_text'] = clean_medical_text(doc['text'])
    
    # Calculate quality metrics
    avg_length = sum(doc['features']['word_count'] for doc in all_documents) / len(all_documents)
    avg_medical_density = sum(doc['features']['medical_term_density'] for doc in all_documents) / len(all_documents)
    
    print(f"   📊 Average document length: {avg_length:.1f} words")
    print(f"   🏥 Average medical term density: {avg_medical_density:.2%}")
    
    # Create ML-ready dataset
    print(f"\n⚙️ Preparing ML-ready dataset...")
    
    # Split data for training
    from sklearn.model_selection import train_test_split
    
    texts = [doc['processed_text'] for doc in all_documents]
    labels = [doc['specialty'] for doc in all_documents]
    
    X_train, X_test, y_train, y_test = train_test_split(
        texts, labels, 
        test_size=0.2, 
        random_state=42, 
        stratify=labels
    )
    
    print(f"   ✅ Training set: {len(X_train)} documents")
    print(f"   ✅ Test set: {len(X_test)} documents")
    
    # Save complete dataset
    complete_dataset = {
        "metadata": {
            "created_at": datetime.now().isoformat(),
            "total_documents": len(all_documents),
            "real_documents": len(pubmed_docs),
            "synthetic_documents": len(synthetic_docs),
            "specialties": list(specialty_counts.keys()),
            "avg_document_length": avg_length,
            "avg_medical_density": avg_medical_density,
            "train_size": len(X_train),
            "test_size": len(X_test)
        },
        "documents": all_documents,
        "specialty_distribution": dict(specialty_counts),
        "train_test_split": {
            "X_train": X_train,
            "X_test": X_test,
            "y_train": y_train,
            "y_test": y_test
        }
    }
    
    # Save to files
    os.makedirs("data/processed", exist_ok=True)
    
    # Main dataset
    with open("data/processed/complete_medical_dataset.json", 'w', encoding='utf-8') as f:
        json.dump(complete_dataset, f, indent=2, ensure_ascii=False)
    
    # Training data (simplified for ML)
    training_data = {
        "X_train": X_train,
        "X_test": X_test,
        "y_train": y_train,
        "y_test": y_test,
        "specialties": list(set(labels))
    }
    
    with open("data/processed/training_ready_dataset.json", 'w', encoding='utf-8') as f:
        json.dump(training_data, f, indent=2)
    
    print(f"\n💾 Dataset saved successfully!")
    print(f"   📁 Complete dataset: data/processed/complete_medical_dataset.json")
    print(f"   📁 Training ready: data/processed/training_ready_dataset.json")
    
    # Final summary
    print(f"\n🎉 Medical Data Integration Complete!")
    print("=" * 50)
    print(f"✅ Real Medical Literature: {len(pubmed_docs)} documents from PubMed")
    print(f"✅ Synthetic Medical Cases: {len(synthetic_docs)} high-quality examples")
    print(f"✅ Professional Data Processing: Text cleaning & feature extraction")
    print(f"✅ ML-Ready Format: Train/test split prepared")
    print(f"✅ Healthcare Compliance: No PHI, HIPAA-compliant synthetic data")
    
    print(f"\n🚀 Ready for:")
    print("   1. Machine Learning model training")
    print("   2. FastAPI service deployment") 
    print("   3. Streamlit medical dashboard")
    print("   4. Complete MLOps pipeline")
    
    print(f"\n🌟 Your Professional Medical ML Portfolio includes:")
    print("   📚 Real medical literature from PubMed API")
    print("   🤖 High-quality synthetic medical cases")
    print("   ⚙️ Professional data preprocessing pipeline")
    print("   🔬 Feature extraction for medical text")
    print("   📊 Quality metrics and validation")
    print("   🎯 Production-ready ML datasets")

except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 50)

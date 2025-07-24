"""
Complete Medical Dataset Processor
=================================

Combines PubMed real data with synthetic data and processes everything
through our medical classification pipeline.
"""

import sys
import json
from pathlib import Path
from datetime import datetime

# Add src to path  
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

print("🏥 Complete Medical Dataset Processor")
print("=" * 50)

try:
    # Load real PubMed data
    pubmed_file = Path("data/pubmed_simple_dataset.json")
    
    if pubmed_file.exists():
        print("📚 Loading PubMed real medical data...")
        with open(pubmed_file, 'r', encoding='utf-8') as f:
            pubmed_docs = json.load(f)
        print(f"   ✅ Loaded {len(pubmed_docs)} real medical abstracts")
    else:
        print("❌ PubMed data not found. Run simple_pubmed_fetch.py first")
        pubmed_docs = []
    
    # Generate synthetic data to supplement
    print("\n🤖 Generating synthetic medical data...")
    from data.ingestion import SyntheticDataGenerator
    
    synthetic_gen = SyntheticDataGenerator(use_openai=False)
    synthetic_docs = []
    
    # Generate 5 synthetic docs per specialty to supplement
    for doc in synthetic_gen.load_data(limit=25):
        synthetic_doc = {
            "id": doc.id,
            "text": doc.text, 
            "specialty": doc.specialty,
            "source": doc.source,
            "confidence": doc.confidence,
            "metadata": doc.metadata
        }
        synthetic_docs.append(synthetic_doc)
    
    print(f"   ✅ Generated {len(synthetic_docs)} synthetic documents")
    
    # Combine datasets
    all_documents = pubmed_docs + synthetic_docs
    print(f"\n📊 Combined dataset: {len(all_documents)} total documents")
    
    # Show distribution
    specialty_counts = {}
    for doc in all_documents:
        specialty = doc['specialty']
        specialty_counts[specialty] = specialty_counts.get(specialty, 0) + 1
    
    print("   Distribution by specialty:")
    for specialty, count in specialty_counts.items():
        real_count = len([d for d in pubmed_docs if d['specialty'] == specialty])
        synthetic_count = len([d for d in synthetic_docs if d['specialty'] == specialty])
        print(f"   {specialty}: {count} total ({real_count} real + {synthetic_count} synthetic)")
    
    # Process through our pipeline
    print(f"\n⚙️ Processing through medical ML pipeline...")
    
    # Convert to our MedicalDocument format
    from data.ingestion import MedicalDocument
    from data.preprocessing import MedicalTextPreprocessor, FeatureExtractor, DataValidator
    
    medical_documents = []
    for doc in all_documents:
        medical_doc = MedicalDocument(
            id=doc['id'],
            text=doc['text'],
            specialty=doc['specialty'], 
            source=doc['source'],
            confidence=doc.get('confidence', 0.9),
            metadata=doc.get('metadata', {})
        )
        medical_documents.append(medical_doc)
    
    print(f"   ✅ Converted to MedicalDocument format")
    
    # Validate data quality
    print("🔍 Validating data quality...")
    validator = DataValidator()
    validation_results = validator.validate_documents(medical_documents)
    
    print(f"   ✅ Quality score: {validation_results['quality_score']:.2%}")
    print(f"   ✅ Valid documents: {validation_results['valid_documents']}/{validation_results['total_documents']}")
    
    # Preprocess documents
    print("🧹 Preprocessing medical text...")
    preprocessor = MedicalTextPreprocessor(preserve_medical_terms=True)
    
    processed_docs = []
    for doc in medical_documents[:validation_results['valid_documents']]:
        processed = preprocessor.preprocess_document(doc)
        processed_docs.append(processed)
    
    print(f"   ✅ Preprocessed {len(processed_docs)} documents")
    
    # Extract ML features
    print("🔬 Extracting machine learning features...")
    feature_extractor = FeatureExtractor(
        max_features=5000,
        ngram_range=(1, 2), 
        min_df=1,
        max_df=0.9
    )
    
    features = feature_extractor.fit_transform(processed_docs)
    
    print(f"   ✅ Extracted {features.text_features.shape[1]} features")
    print(f"   📊 Feature matrix: {features.text_features.shape}")
    print(f"   🏷️ Classes: {list(feature_extractor.label_encoder.classes_)}")
    
    # Save processed dataset
    print("\n💾 Saving complete processed dataset...")
    
    # Save the complete dataset 
    complete_dataset = {
        "metadata": {
            "created_at": datetime.now().isoformat(),
            "total_documents": len(all_documents),
            "real_pubmed_docs": len(pubmed_docs),
            "synthetic_docs": len(synthetic_docs),
            "valid_documents": validation_results['valid_documents'],
            "quality_score": validation_results['quality_score'],
            "features_extracted": features.text_features.shape[1],
            "specialties": list(specialty_counts.keys())
        },
        "documents": all_documents,
        "validation_results": validation_results,
        "specialty_distribution": specialty_counts
    }
    
    output_file = Path("data/complete_medical_dataset.json")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(complete_dataset, f, indent=2, ensure_ascii=False)
    
    print(f"   ✅ Complete dataset saved to {output_file}")
    
    # Final summary
    print(f"\n🎉 Medical Dataset Processing Complete!")
    print("=" * 50)
    print(f"📄 Total Documents: {len(all_documents)}")
    print(f"🏥 Real Medical Literature: {len(pubmed_docs)} (PubMed)")
    print(f"🤖 Synthetic Medical Cases: {len(synthetic_docs)}")
    print(f"✅ Data Quality Score: {validation_results['quality_score']:.2%}")
    print(f"🔬 ML Features: {features.text_features.shape[1]}")
    print(f"🎯 Ready for Model Training: YES!")
    
    print(f"\n🚀 Next Steps:")
    print("   1. Train classification models")
    print("   2. Set up FastAPI service")
    print("   3. Build Streamlit dashboard")
    print("   4. Deploy the complete system")
    
    print(f"\n📊 Your portfolio now includes:")
    print("   ✅ Real medical data from PubMed")
    print("   ✅ Synthetic medical cases") 
    print("   ✅ Professional data processing pipeline")
    print("   ✅ ML-ready feature extraction")
    print("   ✅ Healthcare compliance (PHI handling)")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 50)

"""
Data Pipeline Demo
==================

Demonstrates the medical data pipeline with synthetic data generation.
This script shows the complete workflow from data generation to feature extraction.

Run this to test the data pipeline components!
"""

import sys
from pathlib import Path
import json
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

print("🏥 Medical Classification Engine - Data Pipeline Demo")
print("=" * 60)

try:
    # Test imports
    print("1️⃣ Testing imports...")
    from data.ingestion import MedicalDocument, SyntheticDataGenerator
    from data.preprocessing import MedicalTextPreprocessor, FeatureExtractor, DataValidator
    print("   ✅ All imports successful")
    
    # Generate synthetic data
    print("\n2️⃣ Generating synthetic medical documents...")
    synthetic_gen = SyntheticDataGenerator(use_openai=False)
    
    documents = []
    for doc in synthetic_gen.load_data(limit=25):  # 5 docs per specialty
        documents.append(doc)
    
    print(f"   ✅ Generated {len(documents)} synthetic documents")
    
    # Show specialty distribution
    specialty_counts = {}
    for doc in documents:
        specialty_counts[doc.specialty] = specialty_counts.get(doc.specialty, 0) + 1
    
    print("   📊 Specialty distribution:")
    for specialty, count in specialty_counts.items():
        print(f"      {specialty}: {count}")
    
    # Validate documents
    print("\n3️⃣ Validating document quality...")
    validator = DataValidator()
    validation_results = validator.validate_documents(documents)
    
    print(f"   ✅ Quality score: {validation_results['quality_score']:.2%}")
    print(f"   ✅ Valid documents: {validation_results['valid_documents']}/{validation_results['total_documents']}")
    
    if validation_results['validation_errors']:
        print("   ⚠️ Validation errors found:")
        for error in validation_results['validation_errors'][:3]:  # Show first 3
            print(f"      {error}")
    
    # Preprocess documents
    print("\n4️⃣ Preprocessing medical text...")
    preprocessor = MedicalTextPreprocessor(preserve_medical_terms=True)
    
    processed_docs = []
    for doc in documents[:validation_results['valid_documents']]:
        processed_doc = preprocessor.preprocess_document(doc)
        processed_docs.append(processed_doc)
    
    print(f"   ✅ Preprocessed {len(processed_docs)} documents")
    
    # Show sample preprocessing results
    if processed_docs:
        sample_doc = processed_docs[0]
        print(f"   📝 Sample processed document:")
        print(f"      Original length: {sample_doc['original_length']}")
        print(f"      Processed length: {sample_doc['processed_length']}")
        print(f"      Medical entities: {len(sample_doc['medical_entities']['medications'])} medications, "
              f"{len(sample_doc['medical_entities']['conditions'])} conditions")
    
    # Extract features
    print("\n5️⃣ Extracting ML features...")
    feature_extractor = FeatureExtractor(
        max_features=1000,  # Smaller for demo
        ngram_range=(1, 2),
        min_df=1,
        max_df=0.9
    )
    
    features = feature_extractor.fit_transform(processed_docs)
    
    print(f"   ✅ Extracted {features.text_features.shape[1]} features")
    print(f"   📊 Feature matrix shape: {features.text_features.shape}")
    print(f"   🏷️ Classes: {feature_extractor.label_encoder.classes_}")
    
    # Save demo results
    print("\n6️⃣ Saving demo results...")
    
    demo_results = {
        "timestamp": datetime.now().isoformat(),
        "total_documents": len(documents),
        "valid_documents": validation_results['valid_documents'],
        "quality_score": validation_results['quality_score'],
        "specialty_distribution": specialty_counts,
        "feature_count": features.text_features.shape[1],
        "processing_stats": {
            "avg_original_length": sum(doc['original_length'] for doc in processed_docs) / len(processed_docs),
            "avg_processed_length": sum(doc['processed_length'] for doc in processed_docs) / len(processed_docs),
            "total_entities_found": sum(
                len(doc['medical_entities']['medications']) + 
                len(doc['medical_entities']['conditions']) +
                len(doc['medical_entities']['procedures'])
                for doc in processed_docs
            )
        }
    }
    
    # Save to demo results file
    results_file = Path("data/demo_results.json")
    results_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(results_file, 'w') as f:
        json.dump(demo_results, f, indent=2)
    
    print(f"   ✅ Results saved to {results_file}")
    
    # Final summary
    print("\n🎉 Data Pipeline Demo Complete!")
    print("-" * 40)
    print(f"📄 Generated {len(documents)} medical documents")
    print(f"🔍 {validation_results['valid_documents']} passed validation")
    print(f"⚙️ Extracted {features.text_features.shape[1]} ML features")
    print(f"🎯 Ready for model training!")
    
    print(f"\n📊 Next steps:")
    print("   • Run model training: python -m src.models.train")
    print("   • Start API server: python -m src.api.main") 
    print("   • Launch dashboard: streamlit run src/dashboard/app.py")
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("💡 Make sure you've installed requirements: pip install -r requirements.txt")
    
except Exception as e:
    print(f"❌ Demo failed: {e}")
    print(f"🔍 Error details: {type(e).__name__}: {str(e)}")
    
    # Print helpful debug info
    print(f"\n🛠️ Debug info:")
    print(f"   Python path: {sys.path[:3]}")  # First 3 paths
    print(f"   Current directory: {Path.cwd()}")
    print(f"   Script location: {Path(__file__).parent}")

"""
Sample Data Generation Script
============================

Creates sample medical data to demonstrate the data pipeline.
This script generates synthetic medical documents for all 5 specialties
without requiring external data sources.

Usage:
    python scripts/generate_sample_data.py --count 100 --output data/processed/
"""

import argparse
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List
import random

# Add src to Python path
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir.parent / "src"))

from data.ingestion import MedicalDocument, SyntheticDataGenerator


def generate_sample_dataset(count_per_specialty: int = 20) -> List[MedicalDocument]:
    """
    Generate sample medical documents for testing.
    
    Args:
        count_per_specialty: Number of documents per specialty
        
    Returns:
        List of sample medical documents
    """
    print(f"🏥 Generating {count_per_specialty} documents per specialty...")
    
    # Medical specialties
    specialties = [
        "Cardiology",
        "Emergency", 
        "Pulmonology",
        "Gastroenterology",
        "Dermatology"
    ]
    
    # Sample medical templates and vocabulary for each specialty
    medical_templates = {
        "Cardiology": {
            "presentations": [
                "chest pain and shortness of breath",
                "palpitations and dizziness", 
                "fatigue and leg swelling",
                "syncope and chest discomfort",
                "hypertensive crisis"
            ],
            "findings": [
                "ECG shows ST elevation",
                "echocardiogram reveals reduced ejection fraction",
                "cardiac enzymes are elevated",
                "stress test is positive",
                "coronary angiogram shows stenosis"
            ],
            "diagnoses": [
                "myocardial infarction",
                "heart failure with reduced ejection fraction",
                "atrial fibrillation with rapid ventricular response", 
                "unstable angina",
                "hypertensive emergency"
            ],
            "treatments": [
                "started on dual antiplatelet therapy",
                "initiated ACE inhibitor and beta-blocker",
                "urgent cardiac catheterization performed",
                "cardioversion completed successfully",
                "titrated antihypertensive medications"
            ]
        },
        "Emergency": {
            "presentations": [
                "motor vehicle accident with multiple trauma",
                "severe abdominal pain and vomiting",
                "altered mental status and fever",
                "respiratory distress and hypoxemia",
                "drug overdose with altered consciousness"
            ],
            "findings": [
                "CT scan shows internal bleeding",
                "vital signs indicate shock",
                "neurological exam suggests intoxication",
                "chest X-ray reveals pneumothorax",
                "laboratory values indicate organ failure"
            ],
            "diagnoses": [
                "polytrauma with hemodynamic instability",
                "acute appendicitis requiring surgery",
                "sepsis with multi-organ dysfunction",
                "tension pneumothorax",
                "acetaminophen toxicity"
            ],
            "treatments": [
                "emergent exploratory laparotomy performed",
                "massive transfusion protocol activated",
                "broad-spectrum antibiotics initiated",
                "needle decompression performed",
                "N-acetylcysteine therapy started"
            ]
        },
        "Pulmonology": {
            "presentations": [
                "progressive dyspnea and productive cough",
                "wheezing and chest tightness",
                "hemoptysis and weight loss",
                "chronic cough and night sweats",
                "acute respiratory failure"
            ],
            "findings": [
                "chest CT shows ground glass opacities",
                "pulmonary function tests reveal obstruction",
                "arterial blood gas shows hypoxemia",
                "sputum culture is positive",
                "bronchoscopy reveals mass lesion"
            ],
            "diagnoses": [
                "chronic obstructive pulmonary disease exacerbation",
                "asthma with acute bronchospasm",
                "lung adenocarcinoma",
                "community-acquired pneumonia",
                "acute respiratory distress syndrome"
            ],
            "treatments": [
                "nebulized bronchodilators and systemic steroids",
                "inhaled corticosteroids optimized",
                "chemotherapy regimen initiated",
                "antimicrobial therapy guided by culture",
                "mechanical ventilation with lung-protective strategy"
            ]
        },
        "Gastroenterology": {
            "presentations": [
                "epigastric pain and early satiety",
                "bloody diarrhea and weight loss",
                "jaundice and abdominal distension",
                "dysphagia and heartburn",
                "right upper quadrant pain"
            ],
            "findings": [
                "upper endoscopy reveals gastric ulcer",
                "colonoscopy shows inflammatory changes",
                "MRCP demonstrates bile duct dilation",
                "barium swallow shows stricture",
                "ultrasound reveals gallstones"
            ],
            "diagnoses": [
                "peptic ulcer disease with H. pylori",
                "inflammatory bowel disease",
                "choledocholithiasis with biliary obstruction",
                "esophageal adenocarcinoma",
                "acute cholecystitis"
            ],
            "treatments": [
                "triple therapy for H. pylori eradication",
                "immunosuppressive therapy initiated",
                "ERCP with stone extraction performed",
                "neoadjuvant chemotherapy planned",
                "laparoscopic cholecystectomy completed"
            ]
        },
        "Dermatology": {
            "presentations": [
                "new pigmented lesion with irregular borders",
                "painful vesicular rash in dermatomal distribution",
                "chronic scaling plaques on extensor surfaces",
                "pruritic papular eruption",
                "non-healing ulcer on lower extremity"
            ],
            "findings": [
                "dermoscopy shows asymmetry and color variation",
                "Tzanck smear is positive for multinucleated cells",
                "biopsy reveals hyperkeratosis and inflammation",
                "patch testing shows contact sensitivity",
                "tissue culture grows methicillin-resistant organisms"
            ],
            "diagnoses": [
                "malignant melanoma",
                "herpes zoster",
                "psoriasis vulgaris",
                "atopic dermatitis",
                "venous stasis ulcer with secondary infection"
            ],
            "treatments": [
                "wide local excision with sentinel node biopsy",
                "antiviral therapy and pain management",
                "topical corticosteroids and phototherapy",
                "emollients and antihistamines",
                "wound care and systemic antibiotics"
            ]
        }
    }
    
    documents = []
    doc_id = 1
    
    for specialty in specialties:
        print(f"  📋 Generating {specialty} cases...")
        templates = medical_templates[specialty]
        
        for i in range(count_per_specialty):
            # Randomly select components for this case
            presentation = random.choice(templates["presentations"])
            finding = random.choice(templates["findings"])
            diagnosis = random.choice(templates["diagnoses"])
            treatment = random.choice(templates["treatments"])
            
            # Generate patient demographics
            age = random.randint(25, 85)
            gender = random.choice(["male", "female"])
            
            # Create case text
            case_text = f"{age}-year-old {gender} presents with {presentation}. " \
                       f"Examination and testing reveal {finding}. " \
                       f"Diagnosis: {diagnosis}. " \
                       f"Management: {treatment}."
            
            # Add some medical complexity
            if random.random() > 0.7:  # 30% chance of complications
                complications = [
                    "Patient developed complications during treatment.",
                    "Follow-up imaging showed progression.",
                    "Laboratory values indicated treatment response.",
                    "Patient required additional interventions.",
                    "Multidisciplinary team consultation obtained."
                ]
                case_text += " " + random.choice(complications)
            
            # Create document
            document = MedicalDocument(
                id=f"sample_{doc_id:04d}",
                text=case_text,
                specialty=specialty,
                source="SampleGeneration",
                confidence=round(random.uniform(0.85, 0.98), 2),
                keywords=None,  # Will be extracted later
                metadata={
                    "patient_age": age,
                    "patient_gender": gender,
                    "case_complexity": "high" if len(case_text) > 200 else "standard",
                    "generated_timestamp": datetime.now().isoformat()
                },
                created_at=datetime.now() - timedelta(days=random.randint(0, 90))
            )
            
            documents.append(document)
            doc_id += 1
    
    print(f"✅ Generated {len(documents)} sample medical documents")
    return documents


def save_documents_json(documents: List[MedicalDocument], output_file: Path) -> None:
    """Save documents to JSON file."""
    
    print(f"💾 Saving documents to {output_file}")
    
    # Convert documents to dictionaries
    docs_data = []
    for doc in documents:
        doc_dict = {
            "id": doc.id,
            "text": doc.text,
            "specialty": doc.specialty,
            "source": doc.source,
            "confidence": doc.confidence,
            "keywords": doc.keywords,
            "metadata": doc.metadata,
            "created_at": doc.created_at.isoformat() if doc.created_at else None
        }
        docs_data.append(doc_dict)
    
    # Save to file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(docs_data, f, indent=2, ensure_ascii=False)
    
    print(f"✅ Saved {len(documents)} documents")


def generate_statistics(documents: List[MedicalDocument]) -> Dict:
    """Generate dataset statistics."""
    
    stats = {
        "total_documents": len(documents),
        "specialty_distribution": {},
        "source_distribution": {},
        "text_length_stats": {},
        "confidence_stats": {},
        "date_range": {}
    }
    
    # Specialty distribution
    for doc in documents:
        specialty = doc.specialty
        stats["specialty_distribution"][specialty] = stats["specialty_distribution"].get(specialty, 0) + 1
    
    # Source distribution  
    for doc in documents:
        source = doc.source
        stats["source_distribution"][source] = stats["source_distribution"].get(source, 0) + 1
    
    # Text length statistics
    lengths = [len(doc.text) for doc in documents]
    stats["text_length_stats"] = {
        "min": min(lengths),
        "max": max(lengths),
        "mean": sum(lengths) / len(lengths),
        "median": sorted(lengths)[len(lengths)//2]
    }
    
    # Confidence statistics
    confidences = [doc.confidence for doc in documents if doc.confidence]
    if confidences:
        stats["confidence_stats"] = {
            "min": min(confidences),
            "max": max(confidences), 
            "mean": sum(confidences) / len(confidences)
        }
    
    # Date range
    dates = [doc.created_at for doc in documents if doc.created_at]
    if dates:
        stats["date_range"] = {
            "earliest": min(dates).isoformat(),
            "latest": max(dates).isoformat()
        }
    
    return stats


def main():
    """Main script execution."""
    
    parser = argparse.ArgumentParser(description="Generate sample medical data")
    
    parser.add_argument("--count", type=int, default=20,
                       help="Number of documents per specialty (default: 20)")
    
    parser.add_argument("--output", type=str, default="data/processed/",
                       help="Output directory (default: data/processed/)")
    
    parser.add_argument("--filename", type=str, default="sample_medical_data.json",
                       help="Output filename (default: sample_medical_data.json)")
    
    parser.add_argument("--stats", action="store_true",
                       help="Generate and display statistics")
    
    args = parser.parse_args()
    
    print("🏥 Medical Classification Engine - Sample Data Generator")
    print("=" * 60)
    
    try:
        # Generate sample documents
        documents = generate_sample_dataset(args.count)
        
        # Save to JSON file
        output_dir = Path(args.output)
        output_file = output_dir / args.filename
        save_documents_json(documents, output_file)
        
        # Generate statistics if requested
        if args.stats:
            print("\n📊 Dataset Statistics:")
            print("-" * 30)
            
            stats = generate_statistics(documents)
            
            print(f"Total Documents: {stats['total_documents']}")
            print("\nSpecialty Distribution:")
            for specialty, count in stats['specialty_distribution'].items():
                print(f"  {specialty}: {count}")
            
            print(f"\nText Length Statistics:")
            print(f"  Average: {stats['text_length_stats']['mean']:.1f} characters")
            print(f"  Range: {stats['text_length_stats']['min']} - {stats['text_length_stats']['max']} characters")
            
            if stats['confidence_stats']:
                print(f"\nConfidence Statistics:")
                print(f"  Average: {stats['confidence_stats']['mean']:.2f}")
                print(f"  Range: {stats['confidence_stats']['min']:.2f} - {stats['confidence_stats']['max']:.2f}")
        
        print(f"\n✅ Sample data generation complete!")
        print(f"📄 Output file: {output_file}")
        print(f"🔢 Total documents: {len(documents)}")
        
    except Exception as e:
        print(f"❌ Error generating sample data: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

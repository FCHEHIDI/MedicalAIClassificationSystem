"""
Edge Case Enhancement for Medical Classification Engine
=====================================================
Adds targeted edge case training data to address specific classification issues
identified in testing, without creating a new dataset.

Target Issues:
1. Cardiology edge cases with neurological symptoms
2. Dermatology edge cases with systemic symptoms
"""

import json
import random
from pathlib import Path
from datetime import datetime

# Enhanced edge case data targeting specific classification issues
CARDIOLOGY_EDGE_CASES = [
    {
        "text": "Patient presents with chest pain and dizziness during exertion. ECG shows sinus tachycardia. Troponin negative but stress test reveals ischemic changes. Coronary angiography recommended for suspected stable angina with vasovagal response.",
        "specialty": "Cardiology",
        "edge_case_type": "cardiac_neurological"
    },
    {
        "text": "67-year-old male with palpitations and lightheadedness. Holter monitor shows episodes of supraventricular tachycardia with brief periods of altered mental status. Electrophysiology study planned for arrhythmia evaluation.",
        "specialty": "Cardiology",
        "edge_case_type": "cardiac_neurological"
    },
    {
        "text": "Patient with heart failure presenting with confusion and altered mental status. BNP elevated at 1500 pg/mL. Echocardiogram shows severe left ventricular dysfunction. Cognitive symptoms likely secondary to decreased cardiac output.",
        "specialty": "Cardiology",
        "edge_case_type": "cardiac_neurological"
    },
    {
        "text": "Acute chest pain with syncope in 45-year-old female. ECG shows prolonged QT interval. Torsades de pointes documented on telemetry. Electrolyte correction and cardiac monitoring initiated for ventricular arrhythmia.",
        "specialty": "Cardiology",
        "edge_case_type": "cardiac_neurological"
    },
    {
        "text": "Patient reports chest discomfort with episodes of near-fainting during physical activity. Exercise stress test positive for ischemia. Cardiac catheterization reveals significant LAD stenosis requiring intervention.",
        "specialty": "Cardiology",
        "edge_case_type": "cardiac_neurological"
    },
    {
        "text": "Hypertensive crisis with headache and visual disturbances. Blood pressure 220/120 mmHg. Cardiac evaluation shows left ventricular hypertrophy. Antihypertensive therapy initiated for hypertensive emergency with cardiac involvement.",
        "specialty": "Cardiology",
        "edge_case_type": "cardiac_neurological"
    }
]

DERMATOLOGY_EDGE_CASES = [
    {
        "text": "Drug-induced hypersensitivity reaction with widespread erythematous rash, fever, and malaise following antibiotic therapy. Skin biopsy shows interface dermatitis. Systemic corticosteroids initiated for severe cutaneous adverse drug reaction.",
        "specialty": "Dermatology",
        "edge_case_type": "dermatologic_systemic"
    },
    {
        "text": "Stevens-Johnson syndrome with extensive mucocutaneous involvement, fever, and constitutional symptoms. Skin detachment less than 10% body surface area. Dermatology consultation for severe cutaneous adverse reaction management.",
        "specialty": "Dermatology",
        "edge_case_type": "dermatologic_systemic"
    },
    {
        "text": "Severe atopic dermatitis exacerbation with secondary bacterial infection, fever, and lymphadenopathy. Skin cultures positive for Staphylococcus aureus. Combined topical and systemic therapy for infected eczema.",
        "specialty": "Dermatology",
        "edge_case_type": "dermatologic_systemic"
    },
    {
        "text": "Pustular psoriasis with high fever, malaise, and widespread pustular eruption. Laboratory shows elevated neutrophils and inflammatory markers. Hospitalization required for severe inflammatory skin condition.",
        "specialty": "Dermatology",
        "edge_case_type": "dermatologic_systemic"
    },
    {
        "text": "Cellulitis with systemic inflammatory response, fever, and rapid spreading erythema. Blood cultures drawn for suspected bacteremia. IV antibiotics initiated for severe soft tissue infection with systemic involvement.",
        "specialty": "Dermatology",
        "edge_case_type": "dermatologic_systemic"
    },
    {
        "text": "Erythema multiforme major with fever, malaise, and mucosal involvement following HSV infection. Target lesions on extremities with oral ulcerations. Dermatologic emergency requiring immediate intervention.",
        "specialty": "Dermatology",
        "edge_case_type": "dermatologic_systemic"
    },
    {
        "text": "Necrotizing fasciitis with rapidly spreading skin necrosis, fever, and septic shock. Surgical debridement and broad-spectrum antibiotics initiated. Dermatologic emergency with life-threatening complications.",
        "specialty": "Dermatology",
        "edge_case_type": "dermatologic_systemic"
    }
]

def enhance_dataset():
    """Add edge case data to existing dataset"""
    
    # Load existing dataset
    dataset_file = Path("data/pubmed_large_dataset.json")
    if not dataset_file.exists():
        print("❌ Large dataset not found. Using simple dataset.")
        dataset_file = Path("data/pubmed_simple_dataset.json")
    
    if not dataset_file.exists():
        print("❌ No dataset found to enhance!")
        return
    
    with open(dataset_file, 'r', encoding='utf-8') as f:
        existing_data = json.load(f)
    
    print(f"📊 Current dataset size: {len(existing_data)} samples")
    
    # Count existing specialties
    specialty_counts = {}
    for item in existing_data:
        specialty = item['specialty']
        specialty_counts[specialty] = specialty_counts.get(specialty, 0) + 1
    
    print("📈 Current specialty distribution:")
    for specialty, count in specialty_counts.items():
        print(f"  {specialty}: {count} samples")
    
    # Prepare edge case additions
    new_cases = []
    
    # Add cardiology edge cases
    for i, case in enumerate(CARDIOLOGY_EDGE_CASES):
        new_case = {
            "id": f"edge_cardiology_{i+1}",
            "text": case["text"],
            "specialty": "Cardiology",
            "source": "Edge_Case_Enhancement",
            "confidence": 0.95,
            "metadata": {
                "edge_case_type": case["edge_case_type"],
                "added_at": datetime.now().isoformat(),
                "purpose": "address_cardiology_neurological_confusion"
            }
        }
        new_cases.append(new_case)
    
    # Add dermatology edge cases
    for i, case in enumerate(DERMATOLOGY_EDGE_CASES):
        new_case = {
            "id": f"edge_dermatology_{i+1}",
            "text": case["text"],
            "specialty": "Dermatology",
            "source": "Edge_Case_Enhancement",
            "confidence": 0.95,
            "metadata": {
                "edge_case_type": case["edge_case_type"],
                "added_at": datetime.now().isoformat(),
                "purpose": "address_dermatology_emergency_confusion"
            }
        }
        new_cases.append(new_case)
    
    # Add new cases to existing data
    enhanced_data = existing_data + new_cases
    
    print(f"➕ Adding {len(new_cases)} targeted edge cases:")
    print(f"  - {len(CARDIOLOGY_EDGE_CASES)} Cardiology edge cases (cardiac + neurological symptoms)")
    print(f"  - {len(DERMATOLOGY_EDGE_CASES)} Dermatology edge cases (dermatologic + systemic symptoms)")
    
    # Save enhanced dataset
    backup_file = dataset_file.with_suffix('.backup.json')
    print(f"💾 Creating backup: {backup_file}")
    
    # Create backup
    with open(backup_file, 'w', encoding='utf-8') as f:
        json.dump(existing_data, f, indent=2, ensure_ascii=False)
    
    # Save enhanced dataset
    with open(dataset_file, 'w', encoding='utf-8') as f:
        json.dump(enhanced_data, f, indent=2, ensure_ascii=False)
    
    print(f"✅ Enhanced dataset saved: {len(enhanced_data)} total samples")
    
    # Show new distribution
    new_specialty_counts = {}
    for item in enhanced_data:
        specialty = item['specialty']
        new_specialty_counts[specialty] = new_specialty_counts.get(specialty, 0) + 1
    
    print("📈 Enhanced specialty distribution:")
    for specialty, count in new_specialty_counts.items():
        old_count = specialty_counts.get(specialty, 0)
        added = count - old_count
        if added > 0:
            print(f"  {specialty}: {count} samples (+{added} edge cases)")
        else:
            print(f"  {specialty}: {count} samples")
    
    return enhanced_data

if __name__ == "__main__":
    print("🚀 Enhancing Medical Classification Dataset with Edge Cases")
    print("=" * 60)
    
    enhanced_data = enhance_dataset()
    
    if enhanced_data:
        print("\n✅ Dataset enhancement complete!")
        print("📋 Next steps:")
        print("  1. Run training script to create improved models")
        print("  2. Test edge cases with new models")
        print("  3. Validate that main case performance is maintained")
        print("\n🎯 Ready for improved edge case handling!")
    else:
        print("❌ Enhancement failed!")

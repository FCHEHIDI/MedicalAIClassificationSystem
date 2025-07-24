"""
Intensive Edge Case Training
============================
Adds significantly more edge case data with specific feature engineering
to address the classification issues we're seeing.
"""

import json
import random
from pathlib import Path
from datetime import datetime

# Much more comprehensive edge case data
CARDIOLOGY_INTENSIVE_CASES = [
    # Chest pain + neurological symptoms - very specific cardiac terminology
    {
        "text": "Acute chest pain with ST elevation in leads V2-V4 accompanied by dizziness and diaphoresis. Cardiac catheterization reveals 90% LAD stenosis. Emergent PCI performed with stent placement. Post-procedure echocardiogram shows anterior wall hypokinesis.",
        "specialty": "Cardiology"
    },
    {
        "text": "Patient with typical angina and presyncope during exertion. Stress echocardiography shows reversible perfusion defect in RCA territory. Coronary angiography planned for suspected ischemic heart disease with vasovagal response component.",
        "specialty": "Cardiology"
    },
    {
        "text": "Acute myocardial infarction with cardiogenic shock and altered mental status. Troponin I elevated to 25 ng/mL. Emergency cardiac catheterization shows total occlusion of left main coronary artery. IABP placement and urgent revascularization.",
        "specialty": "Cardiology"
    },
    {
        "text": "Ventricular tachycardia with hemodynamic compromise and confusion. Wide complex tachycardia at 200 bpm. Cardioversion performed with return to sinus rhythm. Electrophysiology study recommended for ICD evaluation.",
        "specialty": "Cardiology"
    },
    {
        "text": "Heart failure exacerbation with orthopnea, PND, and cognitive impairment. BNP 2800 pg/mL. Echocardiogram shows severe LV systolic dysfunction with EF 20%. Diuresis and ACE inhibitor optimization initiated.",
        "specialty": "Cardiology"
    },
    {
        "text": "Atrial fibrillation with rapid ventricular response and lightheadedness. Heart rate 160 bpm irregularly irregular. Rate control achieved with metoprolol. Anticoagulation with warfarin for stroke prevention.",
        "specialty": "Cardiology"
    },
    {
        "text": "Cardiac arrest with successful ROSC and post-arrest confusion. Initial rhythm was ventricular fibrillation. Therapeutic hypothermia protocol initiated. Coronary angiography shows patent vessels.",
        "specialty": "Cardiology"
    },
    {
        "text": "Aortic stenosis with effort syncope and chest pain. Peak gradient 80 mmHg on echocardiography. Valve area 0.8 cm². TAVR evaluation for severe symptomatic aortic stenosis.",
        "specialty": "Cardiology"
    },
    {
        "text": "Pericarditis with chest pain and pericardial friction rub. ECG shows diffuse ST elevation and PR depression. Echocardiogram reveals small pericardial effusion. NSAIDs and colchicine therapy initiated.",
        "specialty": "Cardiology"
    },
    {
        "text": "Hypertensive emergency with chest discomfort and headache. Blood pressure 220/120 mmHg. Cardiac enzymes elevated. Emergency department management with nicardipine drip for controlled reduction.",
        "specialty": "Cardiology"
    }
]

DERMATOLOGY_INTENSIVE_CASES = [
    # Skin conditions with systemic symptoms - very specific dermatologic terminology
    {
        "text": "Severe allergic contact dermatitis with widespread vesicular eruption, fever, and malaise. Patch testing reveals nickel sensitivity. Topical corticosteroids and oral antihistamines prescribed for contact allergic reaction.",
        "specialty": "Dermatology"
    },
    {
        "text": "Drug reaction with eosinophilia and systemic symptoms (DRESS) presenting with morbilliform eruption, fever, and hepatitis. Skin biopsy shows interface dermatitis. Offending medication discontinued, systemic corticosteroids initiated.",
        "specialty": "Dermatology"
    },
    {
        "text": "Erythroderma with generalized scaling, fever, and lymphadenopathy. Skin biopsy shows psoriasiform dermatitis. Temperature regulation impaired due to extensive cutaneous involvement. Hospitalization for supportive care.",
        "specialty": "Dermatology"
    },
    {
        "text": "Bullous pemphigoid with tense bullae, pruritus, and constitutional symptoms. Direct immunofluorescence shows linear IgG and C3 deposition. High-dose systemic corticosteroids for autoimmune blistering disease.",
        "specialty": "Dermatology"
    },
    {
        "text": "Pemphigus vulgaris with flaccid bullae, mucosal erosions, and fever. Positive Nikolsky sign. Indirect immunofluorescence shows intercellular IgG. Immunosuppressive therapy with methotrexate and prednisone.",
        "specialty": "Dermatology"
    },
    {
        "text": "Severe hidradenitis suppurativa with multiple abscesses, fever, and systemic inflammation. Hurley stage III disease. Surgical drainage and adalimumab therapy for chronic inflammatory skin condition.",
        "specialty": "Dermatology"
    },
    {
        "text": "Toxic epidermal necrolysis with skin detachment and fever. SCORTEN score calculated for mortality risk. Intensive care management for severe cutaneous adverse drug reaction with systemic involvement.",
        "specialty": "Dermatology"
    },
    {
        "text": "Acute generalized exanthematous pustulosis with sterile pustules and fever. Drug-related pustular eruption following antibiotic therapy. Supportive care and topical corticosteroids for acute pustular drug reaction.",
        "specialty": "Dermatology"
    },
    {
        "text": "Cutaneous T-cell lymphoma with tumor stage lesions and B symptoms. Skin biopsy shows atypical lymphoid infiltrate. Staging workup and radiation therapy planned for mycosis fungoides.",
        "specialty": "Dermatology"
    },
    {
        "text": "Pyoderma gangrenosum with rapidly expanding ulceration and systemic symptoms. Pathergy test positive. Associated with inflammatory bowel disease. Systemic immunosuppression with cyclosporine.",
        "specialty": "Dermatology"
    }
]

def intensive_dataset_enhancement():
    """Add intensive edge case training data"""
    
    # Load existing dataset
    dataset_file = Path("data/pubmed_large_dataset.json")
    if not dataset_file.exists():
        print("❌ Dataset not found!")
        return
    
    with open(dataset_file, 'r', encoding='utf-8') as f:
        existing_data = json.load(f)
    
    print(f"📊 Current dataset size: {len(existing_data)} samples")
    
    # Prepare intensive edge cases
    new_cases = []
    
    # Add cardiology intensive cases
    for i, case in enumerate(CARDIOLOGY_INTENSIVE_CASES):
        new_case = {
            "id": f"intensive_cardiology_{i+1}",
            "text": case["text"],
            "specialty": "Cardiology",
            "source": "Intensive_Edge_Case_Training",
            "confidence": 0.98,
            "metadata": {
                "edge_case_type": "cardiac_neurological_intensive",
                "added_at": datetime.now().isoformat(),
                "purpose": "intensive_cardiology_edge_case_training"
            }
        }
        new_cases.append(new_case)
    
    # Add dermatology intensive cases
    for i, case in enumerate(DERMATOLOGY_INTENSIVE_CASES):
        new_case = {
            "id": f"intensive_dermatology_{i+1}",
            "text": case["text"],
            "specialty": "Dermatology",
            "source": "Intensive_Edge_Case_Training",
            "confidence": 0.98,
            "metadata": {
                "edge_case_type": "dermatologic_systemic_intensive",
                "added_at": datetime.now().isoformat(),
                "purpose": "intensive_dermatology_edge_case_training"
            }
        }
        new_cases.append(new_case)
    
    # Add new cases to existing data
    enhanced_data = existing_data + new_cases
    
    print(f"➕ Adding {len(new_cases)} intensive edge cases:")
    print(f"  - {len(CARDIOLOGY_INTENSIVE_CASES)} Cardiology intensive cases")
    print(f"  - {len(DERMATOLOGY_INTENSIVE_CASES)} Dermatology intensive cases")
    
    # Save enhanced dataset
    backup_file = dataset_file.with_suffix('.intensive_backup.json')
    print(f"💾 Creating backup: {backup_file}")
    
    # Create backup
    with open(backup_file, 'w', encoding='utf-8') as f:
        json.dump(existing_data, f, indent=2, ensure_ascii=False)
    
    # Save enhanced dataset
    with open(dataset_file, 'w', encoding='utf-8') as f:
        json.dump(enhanced_data, f, indent=2, ensure_ascii=False)
    
    print(f"✅ Enhanced dataset saved: {len(enhanced_data)} total samples")
    return enhanced_data

if __name__ == "__main__":
    print("🎯 Intensive Edge Case Training Enhancement")
    print("=" * 50)
    
    enhanced_data = intensive_dataset_enhancement()
    
    if enhanced_data:
        print("\n✅ Intensive enhancement complete!")
        print("📋 Next step: Retrain models with intensive edge case data")
        print("🎯 Ready for intensive edge case training!")
    else:
        print("❌ Intensive enhancement failed!")

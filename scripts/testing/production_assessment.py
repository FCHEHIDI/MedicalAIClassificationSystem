"""
Production Readiness Assessment
==============================
Final assessment of our Medical Classification Engine for deployment.
"""

import json
import requests
import time
from pathlib import Path

API_URL = "http://localhost:8001"

# Representative test cases for each specialty (main cases)
MAIN_CASES = [
    {
        "specialty": "Cardiology",
        "text": "Patient presents with acute chest pain radiating to left arm. ECG shows ST elevation in leads II, III, aVF. Troponin levels elevated at 2.5 ng/mL. Immediate cardiac catheterization recommended for suspected STEMI."
    },
    {
        "specialty": "Emergency", 
        "text": "Trauma patient from motor vehicle accident. Primary survey shows airway patent, breathing adequate, circulation stable. Secondary survey reveals possible splenic laceration on CT scan."
    },
    {
        "specialty": "Pulmonology",
        "text": "Chronic obstructive pulmonary disease exacerbation with increased dyspnea, purulent sputum production. Chest X-ray shows hyperinflation. Treatment with bronchodilators and systemic corticosteroids."
    },
    {
        "specialty": "Gastroenterology",
        "text": "Upper endoscopy reveals multiple gastric ulcers with active bleeding. H. pylori positive. Treatment includes PPI therapy, clarithromycin, and amoxicillin for eradication."
    },
    {
        "specialty": "Dermatology",
        "text": "Asymmetric pigmented lesion on back measuring 8mm with irregular borders and color variation. Dermoscopy shows atypical network. Urgent excisional biopsy recommended to rule out melanoma."
    }
]

# Edge cases within scope
MANAGEABLE_EDGE_CASES = [
    {
        "specialty": "Cardiology",
        "text": "Heart failure patient presenting with acute exacerbation, bilateral pedal edema, and orthopnea. Chest X-ray shows pulmonary congestion. BNP markedly elevated."
    },
    {
        "specialty": "Emergency",
        "text": "Multi-trauma patient with head injury, chest trauma, and abdominal pain. Glasgow Coma Scale 12, possible internal bleeding, requiring immediate surgical evaluation."
    },
    {
        "specialty": "Pulmonology", 
        "text": "Patient with chronic cough, weight loss, and night sweats. Chest CT shows multiple pulmonary nodules. Sputum cytology and bronchoscopy with biopsy planned."
    },
    {
        "specialty": "Gastroenterology",
        "text": "Patient with severe abdominal pain, nausea, vomiting, and elevated lipase. CT shows pancreatic inflammation consistent with acute pancreatitis. NPO and pain management initiated."
    },
    {
        "specialty": "Dermatology",
        "text": "Drug-induced hypersensitivity reaction with widespread skin rash, fever, and constitutional symptoms following antibiotic therapy. Dermatologic consultation required."
    }
]

def test_api():
    """Test if API is running"""
    try:
        response = requests.get(f"{API_URL}/health", timeout=5)
        return response.status_code == 200
    except:
        return False

def make_prediction(text):
    """Make prediction via API"""
    try:
        response = requests.post(
            f"{API_URL}/predict",
            json={"text": text},
            timeout=10
        )
        if response.status_code == 200:
            return response.json(), None
        else:
            return None, f"API Error: {response.status_code}"
    except Exception as e:
        return None, f"Error: {str(e)}"

def run_production_assessment():
    """Run production readiness assessment"""
    
    print("🏥 Medical Classification Engine - Production Assessment")
    print("=" * 60)
    
    # Check API
    if not test_api():
        print("❌ API is not running!")
        return
    
    print("✅ API is running")
    print()
    
    # Test main cases
    print("🎯 MAIN CASE TESTING")
    print("-" * 25)
    
    main_results = []
    for case in MAIN_CASES:
        result, error = make_prediction(case["text"])
        if result:
            correct = result["specialty"].lower() == case["specialty"].lower()
            main_results.append({
                "specialty": case["specialty"],
                "predicted": result["specialty"],
                "confidence": result["confidence"],
                "correct": correct
            })
            
            emoji = "✅" if correct else "❌"
            print(f"{emoji} {case['specialty']}: {result['specialty']} ({result['confidence']:.1%})")
        else:
            print(f"❌ {case['specialty']}: Error - {error}")
        
        time.sleep(0.1)
    
    # Test manageable edge cases
    print("\n🔧 MANAGEABLE EDGE CASE TESTING")
    print("-" * 35)
    
    edge_results = []
    for case in MANAGEABLE_EDGE_CASES:
        result, error = make_prediction(case["text"])
        if result:
            correct = result["specialty"].lower() == case["specialty"].lower()
            edge_results.append({
                "specialty": case["specialty"],
                "predicted": result["specialty"],
                "confidence": result["confidence"],
                "correct": correct
            })
            
            emoji = "✅" if correct else "❌"
            print(f"{emoji} {case['specialty']} Edge: {result['specialty']} ({result['confidence']:.1%})")
        else:
            print(f"❌ {case['specialty']} Edge: Error - {error}")
        
        time.sleep(0.1)
    
    # Analysis
    print("\n📊 PRODUCTION READINESS ANALYSIS")
    print("=" * 35)
    
    # Main case performance
    main_correct = sum(1 for r in main_results if r["correct"])
    main_total = len(main_results)
    main_accuracy = (main_correct / main_total) * 100 if main_total > 0 else 0
    
    print(f"Main Cases: {main_correct}/{main_total} ({main_accuracy:.1f}%)")
    
    # Edge case performance  
    edge_correct = sum(1 for r in edge_results if r["correct"])
    edge_total = len(edge_results)
    edge_accuracy = (edge_correct / edge_total) * 100 if edge_total > 0 else 0
    
    print(f"Manageable Edge Cases: {edge_correct}/{edge_total} ({edge_accuracy:.1f}%)")
    
    # Overall assessment
    total_correct = main_correct + edge_correct
    total_tests = main_total + edge_total
    overall_accuracy = (total_correct / total_tests) * 100 if total_tests > 0 else 0
    
    print(f"Overall: {total_correct}/{total_tests} ({overall_accuracy:.1f}%)")
    
    # Confidence analysis
    all_results = main_results + edge_results
    if all_results:
        avg_confidence = sum(r["confidence"] for r in all_results) / len(all_results)
        high_conf_count = sum(1 for r in all_results if r["confidence"] >= 0.7)
        
        print(f"Average Confidence: {avg_confidence:.1%}")
        print(f"High Confidence (≥70%): {high_conf_count}/{total_tests} ({high_conf_count/total_tests*100:.1f}%)")
    
    # Production decision
    print("\n🚀 PRODUCTION DEPLOYMENT DECISION")
    print("=" * 35)
    
    if main_accuracy >= 90:
        print("✅ MAIN CASES: Ready for production")
    else:
        print("❌ MAIN CASES: Needs improvement")
    
    if edge_accuracy >= 70:
        print("✅ EDGE CASES: Acceptable performance")
    elif edge_accuracy >= 50:
        print("⚠️ EDGE CASES: Deploy with caution")
    else:
        print("❌ EDGE CASES: Significant issues")
    
    # Final recommendation
    print("\n💡 DEPLOYMENT RECOMMENDATION")
    print("-" * 30)
    
    if main_accuracy >= 90 and overall_accuracy >= 80:
        print("🎯 RECOMMENDED FOR DEPLOYMENT")
        print("   - Excellent main case performance")
        print("   - Acceptable overall performance")
        print("   - Document edge case limitations")
    elif main_accuracy >= 90:
        print("⚡ CONDITIONAL DEPLOYMENT")
        print("   - Strong main case performance")
        print("   - Deploy with edge case warnings")
        print("   - Continue training for edge cases")
    else:
        print("🔄 NEEDS MORE TRAINING")
        print("   - Main case performance insufficient")
        print("   - Delay deployment for improvements")
    
    # Specific recommendations
    print("\n📋 SPECIFIC RECOMMENDATIONS")
    print("-" * 30)
    
    failed_main = [r for r in main_results if not r["correct"]]
    failed_edge = [r for r in edge_results if not r["correct"]]
    
    if failed_main:
        print("🔥 HIGH PRIORITY - Main case failures:")
        for fail in failed_main:
            print(f"   - {fail['specialty']} → {fail['predicted']}")
    
    if failed_edge:
        print("📝 MEDIUM PRIORITY - Edge case failures:")
        for fail in failed_edge:
            print(f"   - {fail['specialty']} → {fail['predicted']}")
    
    if not failed_main and not failed_edge:
        print("🎉 ALL TESTS PASSING - Ready for production!")
    
    print("\n🎯 NEXT STEPS")
    print("-" * 15)
    if overall_accuracy >= 80:
        print("1. Deploy to staging environment")
        print("2. Run comprehensive testing")
        print("3. Document known limitations")
        print("4. Deploy to production")
        print("5. Monitor and collect feedback")
    else:
        print("1. Address main case failures")
        print("2. Collect more edge case data")
        print("3. Retrain models")
        print("4. Re-run assessment")
    
    return {
        "main_accuracy": main_accuracy,
        "edge_accuracy": edge_accuracy,
        "overall_accuracy": overall_accuracy,
        "deployment_ready": overall_accuracy >= 80
    }

if __name__ == "__main__":
    assessment = run_production_assessment()

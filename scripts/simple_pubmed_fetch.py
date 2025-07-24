"""
Simple PubMed Data Fetcher
=========================

Fetches medical abstracts from PubMed for each specialty.
Simplified version that focuses on getting the data working.
"""

import sys
import json
from pathlib import Path
from datetime import datetime

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

print("🏥 Simple PubMed Data Fetcher")
print("=" * 40)

try:
    from Bio import Entrez
    import time
    
    # Configure with your email
    Entrez.email = "fareschehidi7@gmail.com"
    print("✅ Email configured")
    
    # Simple search terms for each specialty - expanded for more variety
    searches = {
        "Cardiology": [
            "cardiology AND heart",
            "coronary artery disease",
            "myocardial infarction",
            "heart failure",
            "arrhythmia",
            "echocardiography",
            "cardiac catheterization",
            "atrial fibrillation"
        ],
        "Emergency": [
            "emergency medicine AND trauma",
            "emergency department",
            "acute care",
            "trauma surgery",
            "sepsis emergency",
            "stroke emergency",
            "cardiac arrest",
            "emergency procedures"
        ],
        "Pulmonology": [
            "pulmonology AND lung",
            "respiratory disease",
            "pneumonia",
            "asthma",
            "chronic obstructive pulmonary",
            "lung cancer",
            "pulmonary embolism",
            "bronchoscopy"
        ],
        "Gastroenterology": [
            "gastroenterology AND digestive",
            "inflammatory bowel disease",
            "hepatology",
            "endoscopy",
            "peptic ulcer",
            "pancreatitis",
            "liver disease",
            "colonoscopy"
        ],
        "Dermatology": [
            "dermatology AND skin",
            "melanoma",
            "psoriasis",
            "skin cancer",
            "atopic dermatitis",
            "dermatopathology",
            "skin biopsy",
            "acne treatment"
        ]
    }
    
    documents = []
    articles_per_specialty = 1000  # PROFESSIONAL SCALE: 1000 per specialty = 5000 total
    
    total_fetched = 0
    target_per_search = articles_per_specialty // len(searches["Cardiology"])  # ~125 per search term
    
    for specialty, search_terms in searches.items():
        print(f"\n📚 Fetching {specialty} abstracts...")
        specialty_docs = 0
        
        for search_term in search_terms:
            if specialty_docs >= articles_per_specialty:
                break
                
            print(f"   Search: {search_term}")
            
            try:
                # Search PubMed
                handle = Entrez.esearch(
                    db="pubmed",
                    term=search_term,
                    retmax=target_per_search,
                    sort="relevance"
                )
                search_results = Entrez.read(handle)
                handle.close()
                
                # Get the IDs safely
                ids = []
                try:
                    if hasattr(search_results, 'get'):
                        ids = search_results.get('IdList', [])
                    elif isinstance(search_results, dict):
                        ids = search_results.get('IdList', [])
                    else:
                        ids = search_results['IdList'] if search_results else []
                except:
                    ids = []
                
                print(f"   Found {len(ids)} article IDs")
                
                if ids:
                    # Fetch the abstracts
                    handle = Entrez.efetch(
                        db="pubmed",
                        id=",".join(ids[:target_per_search]),
                        rettype="abstract",
                        retmode="text"
                    )
                    
                    abstracts_text = handle.read()
                    handle.close()
                    
                    print(f"   ✅ Retrieved {len(abstracts_text)} characters")
                    
                    # Split into individual abstracts (simple approach)
                    abstract_parts = abstracts_text.split('\n\n\n')
                    
                    for i, abstract_text in enumerate(abstract_parts):
                        if len(abstract_text.strip()) > 200 and specialty_docs < articles_per_specialty:
                            # Clean up the text
                            clean_text = abstract_text.strip()
                            clean_text = ' '.join(clean_text.split())  # Normalize whitespace
                            
                            # Create document
                            doc = {
                                "id": f"pubmed_{specialty.lower()}_{total_fetched+1}",
                                "text": clean_text,
                                "specialty": specialty,
                                "source": "PubMed",
                                "confidence": 0.90,
                                "metadata": {
                                    "search_term": search_term,
                                    "fetched_at": datetime.now().isoformat(),
                                    "pmid_batch": ",".join(ids[:3])  # First few IDs
                                }
                            }
                            
                            documents.append(doc)
                            specialty_docs += 1
                            total_fetched += 1
                
                # Be nice to NCBI servers
                time.sleep(0.5)  # Slightly longer delay for more requests
                
            except Exception as e:
                print(f"   ⚠️ Error with '{search_term}': {e}")
                continue
        
        print(f"   ✅ Total {specialty} articles: {specialty_docs}")
        
        # Longer pause between specialties
        time.sleep(1.0)
    
    # Save the results
    if documents:
        output_file = Path("data/pubmed_simple_dataset.json")
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(documents, f, indent=2, ensure_ascii=False)
        
        # Show statistics
        specialty_counts = {}
        total_chars = 0
        
        for doc in documents:
            specialty = doc['specialty']
            specialty_counts[specialty] = specialty_counts.get(specialty, 0) + 1
            total_chars += len(doc['text'])
        
        avg_length = total_chars / len(documents) if documents else 0
        
        print(f"\n🎉 Success!")
        print("=" * 30)
        print(f"📄 Total articles: {len(documents)}")
        print(f"💾 Saved to: {output_file}")
        print(f"📊 Average length: {avg_length:.0f} characters")
        
        print(f"\n📋 By specialty:")
        for specialty, count in specialty_counts.items():
            print(f"   {specialty}: {count}")
        
        # Show a sample
        if documents:
            sample = documents[0]
            print(f"\n📄 Sample abstract ({sample['specialty']}):")
            print(f"   {sample['text'][:200]}...")
        
        print(f"\n✅ Real medical data ready for training!")
        
    else:
        print("❌ No documents retrieved")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 40)

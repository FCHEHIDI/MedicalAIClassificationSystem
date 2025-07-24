"""
Scaled PubMed Data Fetcher - 1500 Samples
========================================
Fetches 1500 medical abstracts (300 per specialty) for robust training.
"""

import json
import sys
from pathlib import Path

print("🏥 Scaling to 1500 Real PubMed Samples")
print("=" * 50)

try:
    from Bio import Entrez
    print("✅ BioPython available")
    
    # Configure email for PubMed API
    Entrez.email = "fareschehidi7@gmail.com"
    
    # Comprehensive search strategies for each specialty with diverse terminology
    search_strategies = {
        "Cardiology": [
            "cardiology[MeSH] OR cardiovascular diseases[MeSH]",
            "myocardial infarction[MeSH] OR heart failure[MeSH]", 
            "coronary artery disease OR cardiac catheterization",
            "arrhythmia OR atrial fibrillation OR pacemaker",
            "echocardiography OR cardiac imaging",
            "heart surgery OR cardiac surgery",
            "hypertension[MeSH] OR cardiac rehabilitation",
            "atherosclerosis OR coronary stenosis",
            "heart valve disease OR cardiac enzymes",
            "chest pain OR angina OR cardiac arrest",
            "heart attack OR cardiovascular disease"
        ],
        "Emergency": [
            "emergency medicine[MeSH] OR emergency service[MeSH]",
            "trauma surgery[MeSH] OR multiple trauma[MeSH]",
            "critical care[MeSH] OR intensive care[MeSH]",
            "sepsis[MeSH] OR septic shock[MeSH]",
            "stroke[MeSH] OR cerebrovascular accident",
            "cardiac arrest OR resuscitation",
            "emergency procedures OR acute care",
            "trauma center OR emergency department",
            "poisoning OR overdose OR emergency surgery",
            "emergency room OR urgent care",
            "accident OR injury OR emergency treatment"
        ],
        "Pulmonology": [
            "pulmonology[MeSH] OR respiratory tract diseases[MeSH]",
            "pneumonia[MeSH] OR lung diseases[MeSH]",
            "asthma[MeSH] OR chronic obstructive pulmonary[MeSH]",
            "lung cancer[MeSH] OR pulmonary neoplasms[MeSH]",
            "pulmonary embolism[MeSH] OR respiratory failure[MeSH]",
            "bronchoscopy OR pulmonary function",
            "respiratory infections OR pneumothorax",
            "lung transplant OR pulmonary fibrosis",
            "sleep apnea OR respiratory therapy",
            "breathing problems OR lung disease",
            "cough OR shortness of breath OR COPD"
        ],
        "Gastroenterology": [
            "gastroenterology[MeSH] OR gastrointestinal diseases[MeSH]",
            "inflammatory bowel diseases[MeSH] OR crohn disease[MeSH]",
            "liver diseases[MeSH] OR hepatitis[MeSH]",
            "endoscopy[MeSH] OR colonoscopy[MeSH]",
            "peptic ulcer[MeSH] OR gastroesophageal reflux[MeSH]",
            "pancreatitis[MeSH] OR pancreatic diseases[MeSH]",
            "gastrointestinal bleeding OR digestive system",
            "colorectal cancer OR gastric cancer",
            "liver transplant OR biliary tract diseases",
            "stomach pain OR abdominal pain",
            "digestive disorders OR bowel disease"
        ],
        "Dermatology": [
            "dermatology[MeSH] OR skin diseases[MeSH]",
            "melanoma[MeSH] OR skin neoplasms[MeSH]",
            "psoriasis[MeSH] OR dermatitis[MeSH]",
            "skin cancer OR basal cell carcinoma",
            "atopic dermatitis OR eczema",
            "dermatopathology OR skin biopsy",
            "acne vulgaris OR dermatologic procedures",
            "wound healing OR skin infection",
            "dermatologic surgery OR cosmetic dermatology",
            "rash OR skin lesion OR skin condition",
            "dermatitis OR skin irritation OR skin care"
        ]
    }
    
    documents = []
    target_per_specialty = 500  # 500 per specialty = 2500 total
    articles_per_search = 60    # ~60 articles per search (500/9 searches)
    
    print(f"🎯 Target: {target_per_specialty * len(search_strategies)} total samples")
    print(f"📊 Strategy: {articles_per_search} articles per search term")
    
    for specialty, search_terms in search_strategies.items():
        print(f"\n🔬 Fetching {specialty} literature...")
        specialty_docs = 0
        
        for i, search_term in enumerate(search_terms):
            if specialty_docs >= target_per_specialty:
                break
                
            print(f"   [{i+1}/{len(search_terms)}] {search_term[:60]}...")
            
            try:
                # Search PubMed
                handle = Entrez.esearch(
                    db="pubmed",
                    term=search_term,
                    retmax=articles_per_search,
                    sort="relevance"
                )
                search_results = Entrez.read(handle)
                handle.close()
                
                if 'IdList' in search_results:
                    ids = search_results['IdList']
                    
                    if ids:
                        # Fetch abstracts in batches
                        batch_size = 10
                        for batch_start in range(0, min(len(ids), articles_per_search), batch_size):
                            batch_end = min(batch_start + batch_size, len(ids))
                            batch_ids = ids[batch_start:batch_end]
                            
                            try:
                                fetch_handle = Entrez.efetch(
                                    db="pubmed",
                                    id=','.join(batch_ids),
                                    rettype="abstract",
                                    retmode="xml"
                                )
                                papers = Entrez.read(fetch_handle)
                                fetch_handle.close()
                                
                                for paper in papers['PubmedArticle']:
                                    if specialty_docs >= target_per_specialty:
                                        break
                                        
                                    try:
                                        medline = paper['MedlineCitation']
                                        article = medline['Article']
                                        
                                        # Extract title
                                        title = str(article.get('ArticleTitle', ''))
                                        
                                        # Extract abstract
                                        abstract_text = ""
                                        if 'Abstract' in article:
                                            abstract_sections = article['Abstract'].get('AbstractText', [])
                                            if isinstance(abstract_sections, list):
                                                abstract_text = ' '.join([str(section) for section in abstract_sections])
                                            else:
                                                abstract_text = str(abstract_sections)
                                        
                                        # Combine title and abstract
                                        full_text = f"{title} {abstract_text}".strip()
                                        
                                        if len(full_text) > 50:  # Ensure meaningful content
                                            documents.append({
                                                'id': f"pubmed_{medline['PMID']}",
                                                'text': full_text,
                                                'specialty': specialty,
                                                'source': 'PubMed',
                                                'confidence': 0.95,
                                                'metadata': {
                                                    'pmid': str(medline['PMID']),
                                                    'search_term': search_term,
                                                    'title': title,
                                                    'abstract_length': len(abstract_text)
                                                }
                                            })
                                            specialty_docs += 1
                                            
                                    except Exception as e:
                                        print(f"      ⚠️  Paper parsing error: {e}")
                                        continue
                                        
                            except Exception as e:
                                print(f"      ❌ Batch fetch error: {e}")
                                continue
                                
                        # Small delay to be respectful to PubMed
                        import time
                        time.sleep(0.5)
                        
                print(f"   ✅ {specialty}: {specialty_docs} documents collected")
                
            except Exception as e:
                print(f"   ❌ Search error for {search_term}: {e}")
                continue
        
        print(f"🎯 {specialty} completed: {specialty_docs}/{target_per_specialty} documents")
    
    # Save the large dataset
    output_file = Path("data/pubmed_large_dataset.json")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(documents, f, indent=2, ensure_ascii=False)
    
    print(f"\n🎉 Large Dataset Created Successfully!")
    print(f"📊 Total documents: {len(documents)}")
    print(f"📁 Saved to: {output_file}")
    
    # Distribution summary
    from collections import Counter
    specialty_counts = Counter([doc['specialty'] for doc in documents])
    print(f"\n📈 Distribution:")
    for specialty, count in specialty_counts.items():
        print(f"   {specialty}: {count} samples")

except ImportError:
    print("❌ BioPython not installed. Installing...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "biopython"])
    print("✅ BioPython installed, please run again")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

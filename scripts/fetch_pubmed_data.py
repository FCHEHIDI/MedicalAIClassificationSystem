"""
PubMed Data Ingestion for Medical Classification
==============================================

Fetch real medical literature from PubMed for each specialty.
This script will create a dataset of medical abstracts categorized
by specialty for training our classification model.
"""

import sys
import json
from pathlib import Path
from datetime import datetime
from typing import List, Dict

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

print("🏥 PubMed Medical Data Ingestion")
print("=" * 50)

try:
    from Bio import Entrez
    import time
    
    # Configure Entrez with your email
    Entrez.email = "fareschehidi7@gmail.com"
    
    print("✅ BioPython configured with your email")
    
    def parse_medline(medline_text: str):
        """Parse MEDLINE format text to extract fields."""
        lines = medline_text.split('\n')
        title = ""
        abstract = ""
        authors = ""
        journal = ""
        
        for line in lines:
            if line.startswith('TI  - '):
                title = line[6:].strip()
            elif line.startswith('AB  - '):
                abstract = line[6:].strip()
            elif line.startswith('AU  - '):
                if authors:
                    authors += "; " + line[6:].strip()
                else:
                    authors = line[6:].strip()
            elif line.startswith('TA  - '):
                journal = line[6:].strip()
        
        return title, abstract, authors, journal
    
    # Medical specialties and their PubMed search terms
    specialty_searches = {
        "Cardiology": [
            "cardiology[MeSH] OR cardiac[Title/Abstract]",
            "heart failure[MeSH] OR myocardial infarction[MeSH]",
            "coronary artery disease[MeSH] OR arrhythmia[MeSH]"
        ],
        "Emergency": [
            "emergency medicine[MeSH] OR trauma[Title/Abstract]", 
            "critical care[MeSH] OR emergency department[Title/Abstract]",
            "resuscitation[MeSH] OR acute care[Title/Abstract]"
        ],
        "Pulmonology": [
            "pulmonology[MeSH] OR respiratory[Title/Abstract]",
            "lung diseases[MeSH] OR asthma[MeSH]",
            "pneumonia[MeSH] OR COPD[Title/Abstract]"
        ],
        "Gastroenterology": [
            "gastroenterology[MeSH] OR digestive system[MeSH]",
            "liver diseases[MeSH] OR inflammatory bowel disease[MeSH]",
            "endoscopy[MeSH] OR gastrointestinal[Title/Abstract]"
        ],
        "Dermatology": [
            "dermatology[MeSH] OR skin diseases[MeSH]",
            "melanoma[MeSH] OR dermatitis[MeSH]",
            "skin cancer[Title/Abstract] OR dermatological[Title/Abstract]"
        ]
    }
    
    # Fetch articles for each specialty
    medical_documents = []
    docs_per_specialty = 10  # Adjust as needed
    
    for specialty, search_terms in specialty_searches.items():
        print(f"\n📚 Fetching {specialty} articles...")
        
        specialty_docs = 0
        for i, search_term in enumerate(search_terms):
            if specialty_docs >= docs_per_specialty:
                break
                
            try:
                print(f"   🔍 Search {i+1}: {search_term[:50]}...")
                
                # Search PubMed
                handle = Entrez.esearch(
                    db="pubmed",
                    term=search_term,
                    retmax=docs_per_specialty // len(search_terms) + 2,
                    sort="relevance"
                )
                search_results = Entrez.read(handle)
                handle.close()
                
                article_ids = search_results['IdList'] if 'IdList' in search_results else []
                print(f"      Found {len(article_ids)} articles")
                
                if not article_ids:
                    continue
                
                # Fetch abstracts
                for article_id in article_ids[:docs_per_specialty // len(search_terms)]:
                    if specialty_docs >= docs_per_specialty:
                        break
                        
                    try:
                        # Get detailed info
                        handle = Entrez.efetch(
                            db="pubmed",
                            id=article_id,
                            rettype="medline",
                            retmode="text"
                        )
                        
                        medline_text = handle.read()
                        handle.close()
                        
                        # Parse the MEDLINE format
                        title, abstract, authors, journal = parse_medline(medline_text)
                        
                        if title and abstract and len(abstract) > 100:
                            # Create document
                            full_text = f"{title} {abstract}"
                            
                            document = {
                                "id": f"pubmed_{article_id}",
                                "text": full_text,
                                "specialty": specialty,
                                "source": "PubMed",
                                "confidence": 0.95,  # High confidence for PubMed data
                                "metadata": {
                                    "pmid": article_id,
                                    "title": title,
                                    "abstract": abstract,
                                    "authors": authors,
                                    "journal": journal,
                                    "search_term": search_term,
                                    "fetched_at": datetime.now().isoformat()
                                }
                            }
                            
                            medical_documents.append(document)
                            specialty_docs += 1
                            
                            print(f"      ✅ Added article {article_id} ({len(full_text)} chars)")
                        
                        # Rate limiting - be nice to NCBI
                        time.sleep(0.1)
                        
                    except Exception as e:
                        print(f"      ⚠️ Error fetching article {article_id}: {e}")
                        continue
                
                # Brief pause between searches
                time.sleep(0.2)
                
            except Exception as e:
                print(f"   ❌ Search failed: {e}")
                continue
        
        print(f"   ✅ Collected {specialty_docs} {specialty} articles")
    
    def parse_medline(medline_text: str):
        """Parse MEDLINE format text to extract fields."""
        lines = medline_text.split('\n')
        title = ""
        abstract = ""
        authors = ""
        journal = ""
        
        for line in lines:
            if line.startswith('TI  - '):
                title = line[6:].strip()
            elif line.startswith('AB  - '):
                abstract = line[6:].strip()
            elif line.startswith('AU  - '):
                if authors:
                    authors += "; " + line[6:].strip()
                else:
                    authors = line[6:].strip()
            elif line.startswith('TA  - '):
                journal = line[6:].strip()
        
        return title, abstract, authors, journal
    
    # Save the dataset
    output_file = Path("data/pubmed_medical_dataset.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(medical_documents, f, indent=2, ensure_ascii=False)
    
    # Generate statistics
    specialty_counts = {}
    total_chars = 0
    
    for doc in medical_documents:
        specialty = doc['specialty']
        specialty_counts[specialty] = specialty_counts.get(specialty, 0) + 1
        total_chars += len(doc['text'])
    
    print(f"\n🎉 PubMed Data Ingestion Complete!")
    print("=" * 50)
    print(f"📄 Total Articles: {len(medical_documents)}")
    print(f"💾 Saved to: {output_file}")
    print(f"📊 Average length: {total_chars / len(medical_documents):.0f} characters")
    
    print(f"\n📋 Specialty Distribution:")
    for specialty, count in specialty_counts.items():
        print(f"   {specialty}: {count} articles")
    
    print(f"\n✅ Ready for model training!")
    print(f"🚀 Next: Run training pipeline with this data")
    
except ImportError as e:
    print(f"❌ Missing dependency: {e}")
    print("💡 Install with: pip install biopython")
    
except Exception as e:
    print(f"❌ Error during ingestion: {e}")
    print(f"🔍 Error type: {type(e).__name__}")

print("\n" + "=" * 50)

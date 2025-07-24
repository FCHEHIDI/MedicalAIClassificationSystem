"""
PubMed Connection Test
=====================

Simple test to verify PubMed API connectivity with your email.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

print("🧪 Testing PubMed Connection")
print("=" * 40)

try:
    print("1. Testing Bio imports...")
    from Bio import Entrez
    print("   ✅ BioPython imported successfully")
    
    print("2. Setting up Entrez with your email...")
    Entrez.email = "fareschehidi7@gmail.com"
    print("   ✅ Email configured")
    
    print("3. Testing PubMed search...")
    # Search for a few cardiology articles
    handle = Entrez.esearch(
        db="pubmed",
        term="cardiology[MeSH] AND heart failure",
        retmax=5
    )
    search_results = Entrez.read(handle)
    handle.close()
    
    print(f"   ✅ Found {len(search_results['IdList'])} articles")
    print(f"   📋 Article IDs: {search_results['IdList'][:3]}...")
    
    if search_results['IdList']:
        print("4. Testing article fetch...")
        # Get one abstract
        handle = Entrez.efetch(
            db="pubmed",
            id=search_results['IdList'][0],
            rettype="abstract",
            retmode="text"
        )
        
        abstract = handle.read()
        handle.close()
        
        print(f"   ✅ Fetched abstract ({len(abstract)} characters)")
        print(f"   📄 Preview: {abstract[:200]}...")
    
    print("\n🎉 PubMed Connection Test PASSED!")
    print("✅ Ready to ingest medical literature data")
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    
except Exception as e:
    print(f"❌ Connection test failed: {e}")
    print("💡 This might be a network issue or API rate limiting")

print("\n" + "="*40)

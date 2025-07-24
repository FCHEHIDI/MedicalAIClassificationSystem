"""
Data Ingestion Module
====================

Handles ingestion from multiple medical data sources:
- MIMIC-III Clinical Database
- PubMed Medical Literature 
- Synthetic medical text generation
- Custom datasets

All data sources are processed to extract medical texts and map them
to the 5 target specialties: Cardiology, Emergency, Pulmonology, 
Gastroenterology, and Dermatology.
"""

import json
import logging
import pandas as pd
import requests
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Generator
from dataclasses import dataclass
from datetime import datetime

import nltk
from Bio import Entrez
from sqlalchemy import create_engine, text
import openai

from ..config import settings
from ..utils.logging import get_logger


logger = get_logger(__name__)


@dataclass
class MedicalDocument:
    """Represents a medical document with classification metadata."""
    
    id: str
    text: str
    specialty: str
    source: str
    confidence: Optional[float] = None
    keywords: Optional[List[str]] = None
    metadata: Optional[Dict] = None
    created_at: Optional[datetime] = None


class BaseDataLoader(ABC):
    """Abstract base class for medical data loaders."""
    
    def __init__(self, specialty_mapping: Optional[Dict[str, str]] = None):
        """
        Initialize data loader.
        
        Args:
            specialty_mapping: Custom mapping from source categories to target specialties
        """
        self.specialty_mapping = specialty_mapping or self._get_default_mapping()
        self.target_specialties = settings.model.medical_specialties
        
    def _get_default_mapping(self) -> Dict[str, str]:
        """Get default specialty mapping for this data source."""
        return {}
        
    @abstractmethod
    def load_data(self, limit: Optional[int] = None) -> Generator[MedicalDocument, None, None]:
        """Load medical documents from the data source."""
        pass
    
    def _map_specialty(self, original_specialty: str) -> Optional[str]:
        """Map original specialty to target specialty."""
        mapped = self.specialty_mapping.get(original_specialty.lower())
        return mapped if mapped in self.target_specialties else None


class MIMICDataLoader(BaseDataLoader):
    """
    Loads clinical notes from MIMIC-III database.
    
    MIMIC-III is a large, freely-available database comprising deidentified 
    health-related data from patients admitted to critical care units.
    """
    
    def __init__(self, connection_string: str, specialty_mapping: Optional[Dict[str, str]] = None):
        """
        Initialize MIMIC data loader.
        
        Args:
            connection_string: Database connection string for MIMIC-III
            specialty_mapping: Custom specialty mapping
        """
        super().__init__(specialty_mapping)
        self.connection_string = connection_string
        self.engine = create_engine(connection_string)
        
    def _get_default_mapping(self) -> Dict[str, str]:
        """Default mapping from MIMIC specialties to target specialties."""
        return {
            # Cardiology mappings
            'cardiology': 'Cardiology',
            'cardiac surgery': 'Cardiology', 
            'cardiovascular': 'Cardiology',
            'ccu': 'Cardiology',  # Cardiac Care Unit
            
            # Emergency mappings
            'emergency': 'Emergency',
            'trauma': 'Emergency',
            'emergency department': 'Emergency',
            'acute care': 'Emergency',
            'intensive care': 'Emergency',
            
            # Pulmonology mappings
            'pulmonology': 'Pulmonology',
            'respiratory': 'Pulmonology', 
            'pneumonia': 'Pulmonology',
            'copd': 'Pulmonology',
            'lung': 'Pulmonology',
            
            # Gastroenterology mappings
            'gastroenterology': 'Gastroenterology',
            'gi': 'Gastroenterology',
            'hepatology': 'Gastroenterology',
            'gastro': 'Gastroenterology',
            
            # Dermatology mappings  
            'dermatology': 'Dermatology',
            'skin': 'Dermatology',
            'dermatologic': 'Dermatology',
        }
    
    def load_data(self, limit: Optional[int] = None) -> Generator[MedicalDocument, None, None]:
        """
        Load clinical notes from MIMIC-III database.
        
        Args:
            limit: Maximum number of documents to load
            
        Yields:
            MedicalDocument: Clinical notes with specialty classification
        """
        logger.info(f"Loading MIMIC-III data (limit: {limit})")
        
        # Query to get clinical notes with service information
        query = """
        SELECT DISTINCT
            n.row_id,
            n.text,
            s.curr_service,
            s.prev_service,
            a.admission_type,
            a.admission_location,
            a.discharge_location,
            p.gender,
            p.dob,
            a.admittime,
            a.dischtime
        FROM noteevents n
        JOIN admissions a ON n.hadm_id = a.hadm_id  
        JOIN services s ON n.hadm_id = s.hadm_id
        JOIN patients p ON n.subject_id = p.subject_id
        WHERE n.category IN ('Discharge summary', 'Physician notes', 'Nursing notes')
        AND n.text IS NOT NULL
        AND LENGTH(n.text) > 100
        """
        
        if limit:
            query += f" LIMIT {limit}"
            
        try:
            with self.engine.connect() as conn:
                result = conn.execute(text(query))
                
                for row in result:
                    # Determine specialty from service information
                    specialty = self._determine_specialty(
                        row.curr_service, 
                        row.prev_service,
                        row.admission_type,
                        row.text
                    )
                    
                    if specialty:
                        doc = MedicalDocument(
                            id=f"mimic_{row.row_id}",
                            text=self._clean_text(row.text),
                            specialty=specialty,
                            source="MIMIC-III",
                            metadata={
                                "service": row.curr_service,
                                "admission_type": row.admission_type,
                                "gender": row.gender,
                                "admit_time": row.admittime,
                                "discharge_time": row.dischtime
                            },
                            created_at=datetime.now()
                        )
                        
                        yield doc
                        
        except Exception as e:
            logger.error(f"Error loading MIMIC data: {e}")
            raise
    
    def _determine_specialty(self, curr_service: str, prev_service: str, 
                           admission_type: str, text: str) -> Optional[str]:
        """Determine specialty from service and text content."""
        
        # First try current service
        if curr_service:
            mapped = self._map_specialty(curr_service)
            if mapped:
                return mapped
                
        # Then try previous service
        if prev_service:
            mapped = self._map_specialty(prev_service) 
            if mapped:
                return mapped
                
        # Finally, use text analysis for keyword matching
        return self._classify_by_keywords(text)
    
    def _classify_by_keywords(self, text: str) -> Optional[str]:
        """Classify specialty based on medical keywords in text."""
        text_lower = text.lower()
        
        # Keyword patterns for each specialty
        specialty_keywords = {
            'Cardiology': [
                'cardiac', 'heart', 'ecg', 'ekg', 'myocardial', 'coronary',
                'arrhythmia', 'hypertension', 'chest pain', 'angina', 'stent',
                'catheterization', 'troponin', 'cardiology'
            ],
            'Emergency': [
                'trauma', 'emergency', 'acute', 'urgent', 'resuscitation',
                'code blue', 'shock', 'overdose', 'accident', 'injury'
            ],
            'Pulmonology': [
                'respiratory', 'lung', 'copd', 'asthma', 'pneumonia', 'dyspnea',
                'bronchitis', 'ventilator', 'oxygen', 'shortness of breath'
            ],
            'Gastroenterology': [
                'gastro', 'abdominal', 'liver', 'stomach', 'intestinal', 'gi',
                'endoscopy', 'colonoscopy', 'hepatitis', 'cirrhosis'
            ],
            'Dermatology': [
                'skin', 'rash', 'lesion', 'dermatitis', 'melanoma', 'biopsy',
                'dermatology', 'wound', 'ulcer'
            ]
        }
        
        # Count keyword matches for each specialty
        specialty_scores = {}
        for specialty, keywords in specialty_keywords.items():
            score = sum(1 for keyword in keywords if keyword in text_lower)
            if score > 0:
                specialty_scores[specialty] = score
        
        # Return specialty with highest score (minimum 2 matches)
        if specialty_scores:
            best_specialty = max(specialty_scores.items(), key=lambda x: x[1])
            if best_specialty[1] >= 2:
                return best_specialty[0]
                
        return None
    
    def _clean_text(self, text: str) -> str:
        """Clean and preprocess medical text."""
        # Remove PHI markers and clean up formatting
        text = text.replace('[**', '').replace('**]', '')
        text = ' '.join(text.split())  # Normalize whitespace
        return text[:10000]  # Truncate very long texts


class PubMedDataLoader(BaseDataLoader):
    """
    Loads medical abstracts from PubMed database.
    
    Uses the Entrez API to search and retrieve medical literature
    abstracts classified by medical specialty.
    """
    
    def __init__(self, email: str, api_key: Optional[str] = None, 
                 specialty_mapping: Optional[Dict[str, str]] = None):
        """
        Initialize PubMed data loader.
        
        Args:
            email: Email for Entrez API access
            api_key: Optional API key for higher rate limits
            specialty_mapping: Custom specialty mapping
        """
        super().__init__(specialty_mapping)
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            
    def load_data(self, limit: Optional[int] = None) -> Generator[MedicalDocument, None, None]:
        """
        Load medical abstracts from PubMed.
        
        Args:
            limit: Maximum number of documents per specialty
            
        Yields:
            MedicalDocument: Medical abstracts with specialty classification
        """
        logger.info(f"Loading PubMed data (limit: {limit})")
        
        # Search terms for each specialty
        search_terms = {
            'Cardiology': [
                'cardiology[MeSH] OR cardiac[Title/Abstract] OR heart disease[Title/Abstract]',
                'myocardial infarction[MeSH] OR coronary artery disease[MeSH]',
                'arrhythmia[MeSH] OR heart failure[MeSH]'
            ],
            'Emergency': [
                'emergency medicine[MeSH] OR trauma[Title/Abstract]',
                'critical care[MeSH] OR intensive care[Title/Abstract]', 
                'emergency department[Title/Abstract] OR acute care[Title/Abstract]'
            ],
            'Pulmonology': [
                'pulmonology[MeSH] OR respiratory system[MeSH]',
                'lung diseases[MeSH] OR asthma[MeSH] OR copd[Title/Abstract]',
                'pneumonia[MeSH] OR respiratory failure[Title/Abstract]'
            ],
            'Gastroenterology': [
                'gastroenterology[MeSH] OR digestive system[MeSH]',
                'liver diseases[MeSH] OR inflammatory bowel disease[MeSH]',
                'endoscopy[MeSH] OR gastrointestinal[Title/Abstract]'
            ],
            'Dermatology': [
                'dermatology[MeSH] OR skin diseases[MeSH]',
                'melanoma[MeSH] OR dermatitis[MeSH]',
                'skin cancer[Title/Abstract] OR skin lesions[Title/Abstract]'
            ]
        }
        
        docs_per_specialty = limit // len(search_terms) if limit else 100
        
        for specialty, terms in search_terms.items():
            logger.info(f"Searching PubMed for {specialty} articles")
            
            for term in terms:
                try:
                    # Search PubMed
                    handle = Entrez.esearch(
                        db="pubmed",
                        term=term,
                        retmax=docs_per_specialty // len(terms),
                        sort="relevance"
                    )
                    search_results = Entrez.read(handle)
                    handle.close()
                    
                    if not search_results["IdList"]:
                        continue
                    
                    # Fetch abstracts
                    ids = search_results["IdList"]
                    handle = Entrez.efetch(
                        db="pubmed",
                        id=ids,
                        rettype="medline",
                        retmode="text"
                    )
                    
                    abstracts = handle.read()
                    handle.close()
                    
                    # Parse abstracts
                    for doc in self._parse_pubmed_abstracts(abstracts, specialty):
                        yield doc
                        
                    # Rate limiting
                    time.sleep(0.1)
                    
                except Exception as e:
                    logger.error(f"Error fetching PubMed data for {specialty}: {e}")
                    continue
    
    def _parse_pubmed_abstracts(self, medline_text: str, specialty: str) -> Generator[MedicalDocument, None, None]:
        """Parse PubMed MEDLINE format abstracts."""
        
        articles = medline_text.split('\n\n')
        
        for article in articles:
            if not article.strip():
                continue
                
            lines = article.split('\n')
            pmid = None
            title = ""
            abstract = ""
            authors = ""
            journal = ""
            
            for line in lines:
                if line.startswith('PMID- '):
                    pmid = line[6:].strip()
                elif line.startswith('TI  - '):
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
            
            if pmid and (title or abstract):
                full_text = f"{title} {abstract}".strip()
                
                if len(full_text) > 50:  # Minimum text length
                    yield MedicalDocument(
                        id=f"pubmed_{pmid}",
                        text=full_text,
                        specialty=specialty,
                        source="PubMed",
                        metadata={
                            "pmid": pmid,
                            "title": title,
                            "authors": authors,
                            "journal": journal
                        },
                        created_at=datetime.now()
                    )


class SyntheticDataGenerator(BaseDataLoader):
    """
    Generates synthetic medical documents using templates and patterns.
    
    Creates realistic medical text data for training when real data
    is limited or for augmentation purposes.
    """
    
    def __init__(self, use_openai: bool = False, openai_api_key: Optional[str] = None):
        """
        Initialize synthetic data generator.
        
        Args:
            use_openai: Whether to use OpenAI API for generation
            openai_api_key: OpenAI API key if using GPT models
        """
        super().__init__()
        self.use_openai = use_openai
        if use_openai and openai_api_key:
            openai.api_key = openai_api_key
        
        # Load medical terminology and templates
        self._load_medical_vocabulary()
        
    def _load_medical_vocabulary(self):
        """Load medical vocabulary and templates for each specialty."""
        
        self.specialty_templates = {
            'Cardiology': {
                'symptoms': ['chest pain', 'shortness of breath', 'palpitations', 'syncope', 'fatigue'],
                'conditions': ['myocardial infarction', 'heart failure', 'arrhythmia', 'angina', 'hypertension'],
                'procedures': ['ECG', 'echocardiogram', 'catheterization', 'angioplasty', 'bypass surgery'],
                'medications': ['aspirin', 'beta-blockers', 'ACE inhibitors', 'statins', 'nitroglycerin'],
                'templates': [
                    "Patient presents with {symptom}. {procedure} shows {finding}. Diagnosed with {condition}. Started on {medication}.",
                    "{age}-year-old {gender} with history of {condition} presents with {symptom}. {procedure} performed. {treatment} recommended.",
                    "Acute {condition} with {symptom} and {finding}. {procedure} shows {result}. {medication} initiated."
                ]
            },
            'Emergency': {
                'symptoms': ['trauma', 'severe pain', 'unconsciousness', 'bleeding', 'difficulty breathing'],
                'conditions': ['fracture', 'laceration', 'overdose', 'shock', 'cardiac arrest'],
                'procedures': ['CT scan', 'X-ray', 'intubation', 'IV access', 'wound repair'],
                'medications': ['morphine', 'epinephrine', 'fluids', 'antibiotics', 'tetanus shot'],
                'templates': [
                    "Patient brought to ED with {symptom}. Vitals show {finding}. {procedure} performed. Diagnosis: {condition}.",
                    "Trauma patient with {condition}. {symptom} noted. Emergency {procedure} completed. {medication} administered.",
                    "Critical patient presenting with {symptom}. {procedure} shows {finding}. Immediate {treatment} required."
                ]
            },
            'Pulmonology': {
                'symptoms': ['cough', 'dyspnea', 'wheezing', 'chest tightness', 'sputum production'],
                'conditions': ['asthma', 'COPD', 'pneumonia', 'pulmonary embolism', 'lung cancer'],
                'procedures': ['chest X-ray', 'CT chest', 'spirometry', 'bronchoscopy', 'arterial blood gas'],
                'medications': ['albuterol', 'steroids', 'antibiotics', 'oxygen therapy', 'bronchodilators'],
                'templates': [
                    "Patient with {condition} presents with {symptom}. {procedure} shows {finding}. {medication} prescribed.",
                    "{age}-year-old smoker with {symptom} and {finding}. {procedure} confirms {condition}. {treatment} started.",
                    "Respiratory failure with {symptom}. {procedure} shows {result}. {medication} and {treatment} initiated."
                ]
            },
            'Gastroenterology': {
                'symptoms': ['abdominal pain', 'nausea', 'vomiting', 'diarrhea', 'bloating'],
                'conditions': ['gastritis', 'peptic ulcer', 'IBD', 'hepatitis', 'gallstones'],
                'procedures': ['endoscopy', 'colonoscopy', 'CT abdomen', 'liver biopsy', 'ERCP'],
                'medications': ['proton pump inhibitors', 'antacids', 'antibiotics', 'steroids', 'anti-emetics'],
                'templates': [
                    "Patient with {symptom} for {duration}. {procedure} reveals {finding}. Diagnosis: {condition}. {medication} started.",
                    "{age}-year-old with {condition} presenting with {symptom}. {procedure} shows {result}. {treatment} recommended.",
                    "GI bleeding with {symptom}. Emergency {procedure} performed. {finding} noted. {medication} given."
                ]
            },
            'Dermatology': {
                'symptoms': ['rash', 'itching', 'lesion', 'discoloration', 'scaling'],
                'conditions': ['eczema', 'psoriasis', 'melanoma', 'basal cell carcinoma', 'dermatitis'],
                'procedures': ['biopsy', 'dermoscopy', 'patch testing', 'excision', 'cryotherapy'],
                'medications': ['topical steroids', 'antihistamines', 'antibiotics', 'retinoids', 'moisturizers'],
                'templates': [
                    "Patient presents with {symptom} on {location}. {procedure} performed. Diagnosis: {condition}. {medication} prescribed.",
                    "New {symptom} noted. {procedure} shows {finding}. Suspicious for {condition}. {treatment} recommended.",
                    "{age}-year-old with {condition} and {symptom}. {procedure} confirms {result}. {medication} started."
                ]
            }
        }
    
    def load_data(self, limit: Optional[int] = None) -> Generator[MedicalDocument, None, None]:
        """
        Generate synthetic medical documents.
        
        Args:
            limit: Total number of documents to generate
            
        Yields:
            MedicalDocument: Synthetic medical documents
        """
        logger.info(f"Generating synthetic medical data (limit: {limit})")
        
        docs_per_specialty = (limit or 500) // len(self.target_specialties)
        
        for specialty in self.target_specialties:
            logger.info(f"Generating {docs_per_specialty} synthetic {specialty} documents")
            
            for i in range(docs_per_specialty):
                if self.use_openai:
                    doc = self._generate_openai_document(specialty, i)
                else:
                    doc = self._generate_template_document(specialty, i)
                    
                if doc:
                    yield doc
    
    def _generate_template_document(self, specialty: str, doc_id: int) -> MedicalDocument:
        """Generate document using templates and vocabulary."""
        import random
        
        template_data = self.specialty_templates[specialty]
        template = random.choice(template_data['templates'])
        
        # Fill template with random medical terms
        text = template.format(
            symptom=random.choice(template_data['symptoms']),
            condition=random.choice(template_data['conditions']),
            procedure=random.choice(template_data['procedures']),
            medication=random.choice(template_data['medications']),
            age=random.randint(25, 85),
            gender=random.choice(['male', 'female']),
            finding=f"{random.choice(['elevated', 'abnormal', 'positive', 'negative'])} {random.choice(['levels', 'results', 'findings'])}",
            result=random.choice(['consistent with', 'suggestive of', 'shows evidence of']),
            treatment=random.choice(['therapy', 'management', 'intervention', 'monitoring']),
            duration=random.choice(['2 days', '1 week', '3 months', '6 months']),
            location=random.choice(['arm', 'leg', 'face', 'back', 'chest'])
        )
        
        return MedicalDocument(
            id=f"synthetic_{specialty.lower()}_{doc_id}",
            text=text,
            specialty=specialty,
            source="Synthetic",
            confidence=0.95,
            created_at=datetime.now()
        )
    
    def _generate_openai_document(self, specialty: str, doc_id: int) -> Optional[MedicalDocument]:
        """Generate document using OpenAI API."""
        try:
            prompt = f"""Generate a realistic medical case note for {specialty}. 
            Include patient presentation, examination findings, diagnosis, and treatment plan. 
            Make it sound like a real medical professional wrote it. 
            Keep it concise (2-3 sentences). Do not include patient names or identifying information."""
            
            response = openai.Completion.create(
                engine="text-davinci-003",
                prompt=prompt,
                max_tokens=200,
                temperature=0.8
            )
            
            text = response.choices[0].text.strip()
            
            if len(text) > 50:
                return MedicalDocument(
                    id=f"openai_{specialty.lower()}_{doc_id}",
                    text=text,
                    specialty=specialty,
                    source="OpenAI-Synthetic",
                    confidence=0.90,
                    created_at=datetime.now()
                )
                
        except Exception as e:
            logger.error(f"Error generating OpenAI document: {e}")
            
        return None


class DataIngestionPipeline:
    """
    Orchestrates data ingestion from multiple sources.
    
    Manages the complete data ingestion pipeline including:
    - Loading from multiple data sources
    - Data validation and quality checks
    - Storage to database
    - Progress tracking and logging
    """
    
    def __init__(self, database_url: str):
        """
        Initialize data ingestion pipeline.
        
        Args:
            database_url: Database connection string
        """
        self.database_url = database_url
        self.loaders: List[BaseDataLoader] = []
        
    def add_loader(self, loader: BaseDataLoader):
        """Add a data loader to the pipeline."""
        self.loaders.append(loader)
        
    def run_pipeline(self, output_path: Optional[Path] = None, 
                    limit_per_source: Optional[int] = None) -> Dict[str, int]:
        """
        Run the complete data ingestion pipeline.
        
        Args:
            output_path: Optional path to save processed data
            limit_per_source: Limit documents per data source
            
        Returns:
            Dict with ingestion statistics
        """
        logger.info("Starting data ingestion pipeline")
        
        stats = {
            'total_documents': 0,
            'by_specialty': {spec: 0 for spec in settings.model.medical_specialties},
            'by_source': {},
            'errors': 0
        }
        
        all_documents = []
        
        for loader in self.loaders:
            loader_name = loader.__class__.__name__
            logger.info(f"Running loader: {loader_name}")
            
            loader_docs = 0
            try:
                for doc in loader.load_data(limit=limit_per_source):
                    all_documents.append(doc)
                    loader_docs += 1
                    stats['total_documents'] += 1
                    stats['by_specialty'][doc.specialty] += 1
                    
                stats['by_source'][loader_name] = loader_docs
                
            except Exception as e:
                logger.error(f"Error in {loader_name}: {e}")
                stats['errors'] += 1
        
        # Save to file if requested
        if output_path:
            self._save_documents(all_documents, output_path)
            
        # Save to database
        self._save_to_database(all_documents)
        
        logger.info(f"Pipeline complete. Processed {stats['total_documents']} documents")
        return stats
    
    def _save_documents(self, documents: List[MedicalDocument], output_path: Path):
        """Save documents to JSON file."""
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        doc_dicts = []
        for doc in documents:
            doc_dict = {
                'id': doc.id,
                'text': doc.text,
                'specialty': doc.specialty,
                'source': doc.source,
                'confidence': doc.confidence,
                'keywords': doc.keywords,
                'metadata': doc.metadata,
                'created_at': doc.created_at.isoformat() if doc.created_at else None
            }
            doc_dicts.append(doc_dict)
            
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(doc_dicts, f, indent=2, ensure_ascii=False)
            
        logger.info(f"Saved {len(documents)} documents to {output_path}")
    
    def _save_to_database(self, documents: List[MedicalDocument]):
        """Save documents to database."""
        # This will be implemented when we create the storage module
        logger.info(f"Would save {len(documents)} documents to database")

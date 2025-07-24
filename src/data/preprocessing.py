"""
Data Preprocessing Module
========================

Handles preprocessing and feature extraction for medical text data.
Includes medical-specific text cleaning, feature engineering, and 
data validation tailored for healthcare documentation.
"""

import re
import string
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
from pathlib import Path

import nltk
import spacy
from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.preprocessing import LabelEncoder, StandardScaler
from textstat import flesch_reading_ease, automated_readability_index
import pandas as pd

from ..config import settings
from ..utils.logging import get_logger
from .ingestion import MedicalDocument

logger = get_logger(__name__)


@dataclass
class ProcessedFeatures:
    """Container for processed medical text features."""
    
    text_features: np.ndarray
    text_metadata: Dict[str, Any]
    labels: np.ndarray
    feature_names: List[str]
    document_ids: List[str]


class MedicalTextPreprocessor:
    """
    Preprocesses medical text with healthcare-specific considerations.
    
    Features:
    - Medical abbreviation expansion
    - PHI removal/masking
    - Medical terminology normalization
    - Clinical text cleaning
    """
    
    def __init__(self, preserve_medical_terms: bool = True):
        """
        Initialize medical text preprocessor.
        
        Args:
            preserve_medical_terms: Whether to preserve medical terminology during cleaning
        """
        self.preserve_medical_terms = preserve_medical_terms
        self.nlp = None
        self._load_medical_resources()
        
    def _load_medical_resources(self):
        """Load medical dictionaries and resources."""
        
        # Common medical abbreviations
        self.medical_abbreviations = {
            # Vital signs and measurements
            'bp': 'blood pressure',
            'hr': 'heart rate', 
            'rr': 'respiratory rate',
            'temp': 'temperature',
            'o2 sat': 'oxygen saturation',
            'bmi': 'body mass index',
            'wbc': 'white blood cell',
            'rbc': 'red blood cell',
            'hgb': 'hemoglobin',
            'hct': 'hematocrit',
            
            # Medical conditions
            'mi': 'myocardial infarction',
            'chf': 'congestive heart failure',
            'copd': 'chronic obstructive pulmonary disease',
            'dm': 'diabetes mellitus',
            'htn': 'hypertension',
            'cad': 'coronary artery disease',
            'afib': 'atrial fibrillation',
            'dvt': 'deep vein thrombosis',
            'pe': 'pulmonary embolism',
            'uti': 'urinary tract infection',
            
            # Procedures and tests
            'ecg': 'electrocardiogram',
            'ekg': 'electrocardiogram',
            'cbc': 'complete blood count',
            'bmp': 'basic metabolic panel',
            'cmp': 'comprehensive metabolic panel',
            'ct': 'computed tomography',
            'mri': 'magnetic resonance imaging',
            'cxr': 'chest x-ray',
            'echo': 'echocardiogram',
            'eeg': 'electroencephalogram',
            
            # Medications
            'asa': 'aspirin',
            'ace i': 'ace inhibitor',
            'arb': 'angiotensin receptor blocker',
            'ppi': 'proton pump inhibitor',
            'nsaid': 'nonsteroidal anti-inflammatory drug',
            
            # Units and timing
            'bid': 'twice daily',
            'tid': 'three times daily',
            'qid': 'four times daily',
            'qd': 'once daily',
            'prn': 'as needed',
            'po': 'by mouth',
            'iv': 'intravenous',
            'im': 'intramuscular',
            'mg': 'milligram',
            'mcg': 'microgram',
            'ml': 'milliliter',
            'l': 'liter'
        }
        
        # Medical specialty keywords for context preservation
        self.specialty_keywords = {
            'Cardiology': [
                'cardiac', 'cardio', 'heart', 'coronary', 'myocardial', 'arrhythmia',
                'ecg', 'ekg', 'echo', 'catheter', 'stent', 'bypass', 'valve',
                'angina', 'infarction', 'ischemia', 'pericardial', 'atrial', 'ventricular'
            ],
            'Emergency': [
                'trauma', 'emergency', 'urgent', 'acute', 'critical', 'shock',
                'resuscitation', 'code', 'intubation', 'ventilation', 'iv', 'fluids',
                'stabilize', 'triage', 'ambulance', 'accident', 'injury'
            ],
            'Pulmonology': [
                'lung', 'pulmonary', 'respiratory', 'breath', 'oxygen', 'ventilator',
                'asthma', 'copd', 'pneumonia', 'bronchitis', 'emphysema', 'fibrosis',
                'pleural', 'thorax', 'dyspnea', 'cough', 'sputum', 'wheeze'
            ],
            'Gastroenterology': [
                'gastro', 'gi', 'digestive', 'abdominal', 'stomach', 'intestinal',
                'liver', 'hepatic', 'pancreatic', 'gallbladder', 'colon', 'rectal',
                'endoscopy', 'colonoscopy', 'biopsy', 'ulcer', 'reflux', 'bowel'
            ],
            'Dermatology': [
                'skin', 'dermal', 'cutaneous', 'lesion', 'rash', 'biopsy',
                'melanoma', 'carcinoma', 'dermatitis', 'eczema', 'psoriasis',
                'wound', 'ulcer', 'pigmented', 'mole', 'excision', 'cryotherapy'
            ]
        }
        
        # Load spaCy model for NLP
        try:
            self.nlp = spacy.load("en_core_web_sm")
        except OSError:
            logger.warning("spaCy model 'en_core_web_sm' not found. Some features may be limited.")
            self.nlp = None
    
    def preprocess_document(self, document: MedicalDocument) -> Dict[str, Any]:
        """
        Preprocess a medical document.
        
        Args:
            document: Medical document to preprocess
            
        Returns:
            Dict containing processed text and metadata
        """
        text = document.text
        
        # Step 1: Remove/mask PHI (Protected Health Information)
        text = self._remove_phi(text)
        
        # Step 2: Expand medical abbreviations
        text = self._expand_abbreviations(text)
        
        # Step 3: Clean and normalize text
        text = self._clean_text(text)
        
        # Step 4: Extract medical entities and concepts
        entities = self._extract_medical_entities(text)
        
        # Step 5: Calculate text statistics
        stats = self._calculate_text_stats(text)
        
        return {
            'processed_text': text,
            'original_length': len(document.text),
            'processed_length': len(text),
            'medical_entities': entities,
            'text_stats': stats,
            'specialty_keywords': self._count_specialty_keywords(text),
            'document_id': document.id,
            'specialty': document.specialty
        }
    
    def _remove_phi(self, text: str) -> str:
        """Remove or mask Protected Health Information."""
        
        # Remove PHI markers (MIMIC-III style)
        text = re.sub(r'\[\*\*[^*]*\*\*\]', '[PHI]', text)
        
        # Remove dates (basic pattern)
        text = re.sub(r'\d{1,2}[/-]\d{1,2}[/-]\d{2,4}', '[DATE]', text)
        text = re.sub(r'\d{1,2}/\d{1,2}', '[DATE]', text)
        
        # Remove phone numbers
        text = re.sub(r'\b\d{3}[-.]?\d{3}[-.]?\d{4}\b', '[PHONE]', text)
        
        # Remove SSN patterns
        text = re.sub(r'\b\d{3}-\d{2}-\d{4}\b', '[SSN]', text)
        
        # Remove email addresses
        text = re.sub(r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b', '[EMAIL]', text)
        
        return text
    
    def _expand_abbreviations(self, text: str) -> str:
        """Expand common medical abbreviations."""
        
        # Convert to lowercase for matching, preserve original case
        text_lower = text.lower()
        
        for abbrev, expansion in self.medical_abbreviations.items():
            # Use word boundaries to avoid partial matches
            pattern = r'\b' + re.escape(abbrev) + r'\b'
            text_lower = re.sub(pattern, expansion, text_lower)
        
        return text_lower
    
    def _clean_text(self, text: str) -> str:
        """Clean and normalize medical text."""
        
        # Remove extra whitespace
        text = re.sub(r'\s+', ' ', text)
        
        # Remove special characters but preserve medical notation
        if not self.preserve_medical_terms:
            # Remove punctuation except hyphens in medical terms
            text = re.sub(r'[^\w\s-]', ' ', text)
        
        # Normalize case
        text = text.lower()
        
        # Remove very short words (but preserve important medical abbreviations)
        words = text.split()
        filtered_words = []
        
        for word in words:
            # Keep words that are:
            # - Length >= 3
            # - Known medical abbreviations
            # - Numbers (dosages, measurements)
            if (len(word) >= 3 or 
                word in self.medical_abbreviations or 
                re.match(r'\d+', word)):
                filtered_words.append(word)
        
        return ' '.join(filtered_words)
    
    def _extract_medical_entities(self, text: str) -> Dict[str, List[str]]:
        """Extract medical entities using NLP."""
        
        entities = {
            'medications': [],
            'conditions': [],
            'procedures': [],
            'anatomy': [],
            'symptoms': []
        }
        
        if self.nlp is None:
            return entities
        
        doc = self.nlp(text)
        
        # Extract named entities
        for ent in doc.ents:
            if ent.label_ in ['PERSON', 'GPE']:  # Filter out potential PHI
                continue
                
            entity_text = ent.text.lower()
            
            # Classify entities based on context and patterns
            if any(med in entity_text for med in ['mg', 'ml', 'tablet', 'dose']):
                entities['medications'].append(entity_text)
            elif any(cond in entity_text for cond in ['itis', 'osis', 'emia', 'pathy']):
                entities['conditions'].append(entity_text)
            elif any(proc in entity_text for proc in ['scopy', 'tomy', 'graphy', 'plasty']):
                entities['procedures'].append(entity_text)
            elif any(anat in entity_text for anat in ['heart', 'lung', 'liver', 'brain']):
                entities['anatomy'].append(entity_text)
        
        # Remove duplicates
        for key in entities:
            entities[key] = list(set(entities[key]))
        
        return entities
    
    def _calculate_text_stats(self, text: str) -> Dict[str, float]:
        """Calculate text readability and complexity statistics."""
        
        if not text.strip():
            return {}
        
        stats = {
            'word_count': len(text.split()),
            'sentence_count': len(re.split(r'[.!?]+', text)),
            'avg_word_length': np.mean([len(word) for word in text.split()]),
            'unique_words': len(set(text.split())),
        }
        
        try:
            stats['flesch_reading_ease'] = flesch_reading_ease(text)
            stats['automated_readability_index'] = automated_readability_index(text)
        except:
            stats['flesch_reading_ease'] = 0
            stats['automated_readability_index'] = 0
        
        if stats['word_count'] > 0:
            stats['lexical_diversity'] = stats['unique_words'] / stats['word_count']
        else:
            stats['lexical_diversity'] = 0
        
        return stats
    
    def _count_specialty_keywords(self, text: str) -> Dict[str, int]:
        """Count keywords for each medical specialty."""
        
        keyword_counts = {}
        
        for specialty, keywords in self.specialty_keywords.items():
            count = sum(1 for keyword in keywords if keyword in text.lower())
            keyword_counts[specialty] = count
            
        return keyword_counts


class FeatureExtractor:
    """
    Extracts machine learning features from processed medical text.
    
    Supports multiple feature extraction methods:
    - TF-IDF vectorization
    - N-gram features
    - Medical entity features
    - Text statistics features
    """
    
    def __init__(self, 
                 max_features: int = 10000,
                 ngram_range: Tuple[int, int] = (1, 3),
                 min_df: int = 2,
                 max_df: float = 0.8):
        """
        Initialize feature extractor.
        
        Args:
            max_features: Maximum number of features for TF-IDF
            ngram_range: Range of n-grams to extract
            min_df: Minimum document frequency for features
            max_df: Maximum document frequency for features
        """
        self.max_features = max_features
        self.ngram_range = ngram_range
        self.min_df = min_df
        self.max_df = max_df
        
        # Initialize vectorizers
        self.tfidf_vectorizer = TfidfVectorizer(
            max_features=max_features,
            ngram_range=ngram_range,
            min_df=min_df,
            max_df=max_df,
            stop_words='english'
        )
        
        self.count_vectorizer = CountVectorizer(
            max_features=max_features // 2,
            ngram_range=(1, 2),
            min_df=min_df,
            max_df=max_df,
            stop_words='english'
        )
        
        self.label_encoder = LabelEncoder()
        self.scaler = StandardScaler()
        
        self.is_fitted = False
    
    def fit_transform(self, documents: List[Dict[str, Any]]) -> ProcessedFeatures:
        """
        Fit extractors and transform documents to features.
        
        Args:
            documents: List of preprocessed documents
            
        Returns:
            ProcessedFeatures object with extracted features
        """
        logger.info(f"Extracting features from {len(documents)} documents")
        
        # Extract text and labels
        texts = [doc['processed_text'] for doc in documents]
        labels = [doc['specialty'] for doc in documents]
        doc_ids = [doc['document_id'] for doc in documents]
        
        # Fit and transform TF-IDF features
        tfidf_features = self.tfidf_vectorizer.fit_transform(texts)
        
        # Extract additional features
        additional_features = self._extract_additional_features(documents)
        
        # Combine features
        if additional_features.shape[1] > 0:
            # Scale additional features
            additional_features = self.scaler.fit_transform(additional_features)
            
            # Combine TF-IDF and additional features
            from scipy.sparse import hstack, csr_matrix
            additional_sparse = csr_matrix(additional_features)
            combined_features = hstack([tfidf_features, additional_sparse])
        else:
            combined_features = tfidf_features
        
        # Encode labels
        encoded_labels = self.label_encoder.fit_transform(labels)
        
        # Get feature names
        feature_names = list(self.tfidf_vectorizer.get_feature_names_out())
        if additional_features.shape[1] > 0:
            feature_names.extend(self._get_additional_feature_names())
        
        self.is_fitted = True
        
        return ProcessedFeatures(
            text_features=combined_features,
            text_metadata={
                'feature_names': feature_names,
                'label_classes': self.label_encoder.classes_.tolist(),
                'n_samples': len(documents),
                'n_features': combined_features.shape[1]
            },
            labels=encoded_labels,
            feature_names=feature_names,
            document_ids=doc_ids
        )
    
    def transform(self, documents: List[Dict[str, Any]]) -> ProcessedFeatures:
        """
        Transform new documents using fitted extractors.
        
        Args:
            documents: List of preprocessed documents
            
        Returns:
            ProcessedFeatures object with extracted features
        """
        if not self.is_fitted:
            raise ValueError("FeatureExtractor must be fitted before transform")
        
        # Extract text and labels
        texts = [doc['processed_text'] for doc in documents]
        labels = [doc['specialty'] for doc in documents]
        doc_ids = [doc['document_id'] for doc in documents]
        
        # Transform TF-IDF features
        tfidf_features = self.tfidf_vectorizer.transform(texts)
        
        # Extract additional features
        additional_features = self._extract_additional_features(documents)
        
        # Combine features
        if additional_features.shape[1] > 0:
            # Scale additional features
            additional_features = self.scaler.transform(additional_features)
            
            # Combine TF-IDF and additional features
            from scipy.sparse import hstack, csr_matrix
            additional_sparse = csr_matrix(additional_features)
            combined_features = hstack([tfidf_features, additional_sparse])
        else:
            combined_features = tfidf_features
        
        # Encode labels
        encoded_labels = self.label_encoder.transform(labels)
        
        return ProcessedFeatures(
            text_features=combined_features,
            text_metadata={
                'n_samples': len(documents),
                'n_features': combined_features.shape[1]
            },
            labels=encoded_labels,
            feature_names=self.tfidf_vectorizer.get_feature_names_out().tolist(),
            document_ids=doc_ids
        )
    
    def _extract_additional_features(self, documents: List[Dict[str, Any]]) -> np.ndarray:
        """Extract additional non-TF-IDF features."""
        
        features = []
        
        for doc in documents:
            doc_features = []
            
            # Text statistics
            stats = doc.get('text_stats', {})
            doc_features.extend([
                stats.get('word_count', 0),
                stats.get('sentence_count', 0),
                stats.get('avg_word_length', 0),
                stats.get('unique_words', 0),
                stats.get('lexical_diversity', 0),
                stats.get('flesch_reading_ease', 0),
                stats.get('automated_readability_index', 0)
            ])
            
            # Specialty keyword counts
            keyword_counts = doc.get('specialty_keywords', {})
            for specialty in settings.model.medical_specialties:
                doc_features.append(keyword_counts.get(specialty, 0))
            
            # Medical entity counts
            entities = doc.get('medical_entities', {})
            doc_features.extend([
                len(entities.get('medications', [])),
                len(entities.get('conditions', [])),
                len(entities.get('procedures', [])),
                len(entities.get('anatomy', [])),
                len(entities.get('symptoms', []))
            ])
            
            features.append(doc_features)
        
        return np.array(features)
    
    def _get_additional_feature_names(self) -> List[str]:
        """Get names for additional features."""
        
        names = [
            'word_count', 'sentence_count', 'avg_word_length', 
            'unique_words', 'lexical_diversity', 'flesch_reading_ease',
            'automated_readability_index'
        ]
        
        # Add specialty keyword features
        for specialty in settings.model.medical_specialties:
            names.append(f'keywords_{specialty.lower()}')
        
        # Add medical entity features
        names.extend([
            'entity_medications', 'entity_conditions', 'entity_procedures',
            'entity_anatomy', 'entity_symptoms'
        ])
        
        return names


class DataValidator:
    """
    Validates data quality and consistency for medical documents.
    
    Performs checks for:
    - Text quality and completeness
    - Label consistency
    - Medical terminology validity
    - Data distribution balance
    """
    
    def __init__(self):
        """Initialize data validator."""
        self.validation_rules = self._setup_validation_rules()
    
    def _setup_validation_rules(self) -> Dict[str, Dict]:
        """Setup validation rules for medical data."""
        
        return {
            'text_quality': {
                'min_length': 10,
                'max_length': 50000,
                'min_words': 5,
                'max_repeated_chars': 10,
                'required_chars': string.ascii_letters
            },
            'label_consistency': {
                'valid_specialties': settings.model.medical_specialties,
                'min_samples_per_class': 10
            },
            'medical_content': {
                'min_medical_terms': 1,
                'max_phi_markers': 5,
                'required_medical_patterns': [
                    r'\b(patient|diagnosis|treatment|symptoms?|condition)\b'
                ]
            }
        }
    
    def validate_documents(self, documents: List[MedicalDocument]) -> Dict[str, Any]:
        """
        Validate a collection of medical documents.
        
        Args:
            documents: List of medical documents to validate
            
        Returns:
            Dict with validation results and statistics
        """
        logger.info(f"Validating {len(documents)} medical documents")
        
        validation_results = {
            'total_documents': len(documents),
            'valid_documents': 0,
            'validation_errors': [],
            'warnings': [],
            'statistics': {},
            'quality_score': 0.0
        }
        
        valid_docs = []
        
        for doc in documents:
            doc_validation = self._validate_single_document(doc)
            
            if doc_validation['is_valid']:
                valid_docs.append(doc)
                validation_results['valid_documents'] += 1
            else:
                validation_results['validation_errors'].extend([
                    f"Doc {doc.id}: {error}" for error in doc_validation['errors']
                ])
            
            if doc_validation['warnings']:
                validation_results['warnings'].extend([
                    f"Doc {doc.id}: {warning}" for warning in doc_validation['warnings']
                ])
        
        # Calculate statistics
        validation_results['statistics'] = self._calculate_validation_statistics(valid_docs)
        
        # Calculate overall quality score
        if len(documents) > 0:
            validation_results['quality_score'] = (
                validation_results['valid_documents'] / len(documents)
            )
        
        logger.info(f"Validation complete. {validation_results['valid_documents']}/{len(documents)} documents valid")
        
        return validation_results
    
    def _validate_single_document(self, doc: MedicalDocument) -> Dict[str, Any]:
        """Validate a single medical document."""
        
        result = {
            'is_valid': True,
            'errors': [],
            'warnings': []
        }
        
        # Text quality checks
        text_validation = self._validate_text_quality(doc.text)
        result['errors'].extend(text_validation['errors'])
        result['warnings'].extend(text_validation['warnings'])
        
        # Label consistency checks
        label_validation = self._validate_label(doc.specialty)
        result['errors'].extend(label_validation['errors'])
        
        # Medical content checks
        content_validation = self._validate_medical_content(doc.text)
        result['warnings'].extend(content_validation['warnings'])
        
        # Set overall validity
        result['is_valid'] = len(result['errors']) == 0
        
        return result
    
    def _validate_text_quality(self, text: str) -> Dict[str, List[str]]:
        """Validate text quality."""
        
        errors = []
        warnings = []
        rules = self.validation_rules['text_quality']
        
        # Length checks
        if len(text) < rules['min_length']:
            errors.append(f"Text too short: {len(text)} < {rules['min_length']}")
        elif len(text) > rules['max_length']:
            warnings.append(f"Text very long: {len(text)} > {rules['max_length']}")
        
        # Word count checks
        words = text.split()
        if len(words) < rules['min_words']:
            errors.append(f"Too few words: {len(words)} < {rules['min_words']}")
        
        # Character repetition checks
        for char in text:
            if text.count(char * rules['max_repeated_chars']) > 0:
                warnings.append(f"Excessive character repetition detected: {char}")
                break
        
        # Required character checks
        if not any(c in text for c in rules['required_chars']):
            errors.append("Text contains no alphabetic characters")
        
        return {'errors': errors, 'warnings': warnings}
    
    def _validate_label(self, specialty: str) -> Dict[str, List[str]]:
        """Validate specialty label."""
        
        errors = []
        rules = self.validation_rules['label_consistency']
        
        if specialty not in rules['valid_specialties']:
            errors.append(f"Invalid specialty: {specialty}")
        
        return {'errors': errors}
    
    def _validate_medical_content(self, text: str) -> Dict[str, List[str]]:
        """Validate medical content quality."""
        
        warnings = []
        rules = self.validation_rules['medical_content']
        
        # Check for medical terminology
        medical_pattern_found = False
        for pattern in rules['required_medical_patterns']:
            if re.search(pattern, text, re.IGNORECASE):
                medical_pattern_found = True
                break
        
        if not medical_pattern_found:
            warnings.append("No medical terminology patterns detected")
        
        # Check for excessive PHI markers
        phi_count = len(re.findall(r'\[PHI\]|\[\*\*[^*]*\*\*\]', text))
        if phi_count > rules['max_phi_markers']:
            warnings.append(f"High number of PHI markers: {phi_count}")
        
        return {'warnings': warnings}
    
    def _calculate_validation_statistics(self, documents: List[MedicalDocument]) -> Dict[str, Any]:
        """Calculate statistics for validated documents."""
        
        if not documents:
            return {}
        
        # Specialty distribution
        specialty_counts = {}
        for doc in documents:
            specialty_counts[doc.specialty] = specialty_counts.get(doc.specialty, 0) + 1
        
        # Text length statistics
        text_lengths = [len(doc.text) for doc in documents]
        
        # Source distribution
        source_counts = {}
        for doc in documents:
            source_counts[doc.source] = source_counts.get(doc.source, 0) + 1
        
        return {
            'specialty_distribution': specialty_counts,
            'source_distribution': source_counts,
            'text_length_stats': {
                'mean': np.mean(text_lengths),
                'median': np.median(text_lengths),
                'min': np.min(text_lengths),
                'max': np.max(text_lengths),
                'std': np.std(text_lengths)
            },
            'total_words': sum(len(doc.text.split()) for doc in documents),
            'avg_words_per_doc': np.mean([len(doc.text.split()) for doc in documents])
        }

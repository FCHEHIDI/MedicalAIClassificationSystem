"""
Data Pipeline for Medical Classification Engine
==============================================

This module handles data ingestion, processing, and preparation for the medical
document classification system. It supports multiple data sources including:

- MIMIC-III clinical notes
- PubMed medical abstracts  
- Synthetic medical text generation
- Custom medical datasets

The pipeline ensures consistent preprocessing and feature extraction across
all data sources while maintaining HIPAA compliance principles.
"""

from .ingestion import (
    MIMICDataLoader,
    PubMedDataLoader,
    SyntheticDataGenerator,
    DataIngestionPipeline
)

from .preprocessing import (
    MedicalTextPreprocessor,
    FeatureExtractor,
    DataValidator
)

from .storage import (
    DatabaseManager,
    FeatureStore,
    DataWarehouse
)

__all__ = [
    # Data ingestion
    "MIMICDataLoader",
    "PubMedDataLoader", 
    "SyntheticDataGenerator",
    "DataIngestionPipeline",
    
    # Preprocessing
    "MedicalTextPreprocessor",
    "FeatureExtractor",
    "DataValidator",
    
    # Storage
    "DatabaseManager",
    "FeatureStore", 
    "DataWarehouse"
]

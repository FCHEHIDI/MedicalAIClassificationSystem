"""
Data Storage Module
==================

Handles database operations, feature store management, and data warehousing
for the medical classification system. Provides abstraction layers for
different storage backends and ensures data consistency.
"""

import json
import logging
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple
from sqlalchemy import (
    create_engine, Column, Integer, String, Text, DateTime, Float, JSON,
    Boolean, ForeignKey, Index, UniqueConstraint
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship, Session
from sqlalchemy.dialects.postgresql import UUID
import pandas as pd
import numpy as np
from pathlib import Path

from ..config import settings
from ..utils.logging import get_logger
from .ingestion import MedicalDocument
from .preprocessing import ProcessedFeatures

logger = get_logger(__name__)

# Database models
Base = declarative_base()


class MedicalDocumentORM(Base):
    """Database model for medical documents."""
    
    __tablename__ = 'medical_documents'
    
    id = Column(String, primary_key=True)
    text = Column(Text, nullable=False)
    specialty = Column(String(50), nullable=False, index=True)
    source = Column(String(50), nullable=False, index=True)
    confidence = Column(Float, nullable=True)
    keywords = Column(JSON, nullable=True)
    metadata = Column(JSON, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, index=True)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    features = relationship("DocumentFeatureORM", back_populates="document", cascade="all, delete-orphan")
    predictions = relationship("PredictionORM", back_populates="document", cascade="all, delete-orphan")
    
    # Indexes for performance
    __table_args__ = (
        Index('idx_specialty_source', 'specialty', 'source'),
        Index('idx_created_at', 'created_at'),
    )


class DocumentFeatureORM(Base):
    """Database model for extracted document features."""
    
    __tablename__ = 'document_features'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    document_id = Column(String, ForeignKey('medical_documents.id'), nullable=False, index=True)
    feature_version = Column(String(20), nullable=False, default='v1.0')
    
    # Text statistics
    word_count = Column(Integer, nullable=True)
    sentence_count = Column(Integer, nullable=True)
    avg_word_length = Column(Float, nullable=True)
    unique_words = Column(Integer, nullable=True)
    lexical_diversity = Column(Float, nullable=True)
    flesch_reading_ease = Column(Float, nullable=True)
    automated_readability_index = Column(Float, nullable=True)
    
    # Medical entity counts
    entity_medications = Column(Integer, default=0)
    entity_conditions = Column(Integer, default=0)
    entity_procedures = Column(Integer, default=0)
    entity_anatomy = Column(Integer, default=0)
    entity_symptoms = Column(Integer, default=0)
    
    # Specialty keyword counts (JSON for flexibility)
    specialty_keywords = Column(JSON, nullable=True)
    medical_entities = Column(JSON, nullable=True)
    
    # Feature vectors (stored as JSON for simplicity)
    tfidf_features = Column(JSON, nullable=True)  # Sparse representation
    additional_features = Column(JSON, nullable=True)
    
    created_at = Column(DateTime, default=datetime.utcnow)
    
    # Relationships
    document = relationship("MedicalDocumentORM", back_populates="features")
    
    __table_args__ = (
        UniqueConstraint('document_id', 'feature_version', name='_doc_feature_version_uc'),
    )


class ModelORM(Base):
    """Database model for trained models."""
    
    __tablename__ = 'models'
    
    id = Column(String, primary_key=True)
    name = Column(String(100), nullable=False)
    version = Column(String(20), nullable=False)
    model_type = Column(String(50), nullable=False)
    hyperparameters = Column(JSON, nullable=True)
    training_data_info = Column(JSON, nullable=True)
    performance_metrics = Column(JSON, nullable=True)
    feature_version = Column(String(20), nullable=False)
    model_path = Column(String(500), nullable=True)  # Path to serialized model
    is_active = Column(Boolean, default=False, index=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    
    # Relationships
    predictions = relationship("PredictionORM", back_populates="model")
    
    __table_args__ = (
        Index('idx_model_active', 'is_active'),
        Index('idx_model_created', 'created_at'),
    )


class PredictionORM(Base):
    """Database model for model predictions."""
    
    __tablename__ = 'predictions'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    document_id = Column(String, ForeignKey('medical_documents.id'), nullable=False, index=True)
    model_id = Column(String, ForeignKey('models.id'), nullable=False, index=True)
    predicted_specialty = Column(String(50), nullable=False)
    confidence_scores = Column(JSON, nullable=True)  # Probabilities for all classes
    prediction_time = Column(DateTime, default=datetime.utcnow, index=True)
    is_correct = Column(Boolean, nullable=True)  # For evaluation
    
    # Relationships
    document = relationship("MedicalDocumentORM", back_populates="predictions")
    model = relationship("ModelORM", back_populates="predictions")


class DatasetORM(Base):
    """Database model for dataset management."""
    
    __tablename__ = 'datasets'
    
    id = Column(String, primary_key=True)
    name = Column(String(100), nullable=False)
    description = Column(Text, nullable=True)
    version = Column(String(20), nullable=False)
    split_type = Column(String(20), nullable=False)  # train, validation, test
    document_count = Column(Integer, nullable=False)
    specialty_distribution = Column(JSON, nullable=True)
    source_distribution = Column(JSON, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    
    __table_args__ = (
        Index('idx_dataset_split', 'split_type'),
    )


class DatabaseManager:
    """
    Manages database connections and operations.
    
    Provides high-level interface for storing and retrieving
    medical documents, features, and model artifacts.
    """
    
    def __init__(self, database_url: Optional[str] = None):
        """
        Initialize database manager.
        
        Args:
            database_url: Database connection URL
        """
        self.database_url = database_url or settings.database.url
        self.engine = create_engine(self.database_url)
        self.SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=self.engine)
        
    def create_tables(self):
        """Create all database tables."""
        logger.info("Creating database tables")
        Base.metadata.create_all(bind=self.engine)
        
    def get_session(self) -> Session:
        """Get database session."""
        return self.SessionLocal()
    
    def store_documents(self, documents: List[MedicalDocument]) -> int:
        """
        Store medical documents in database.
        
        Args:
            documents: List of medical documents to store
            
        Returns:
            Number of documents stored
        """
        logger.info(f"Storing {len(documents)} documents to database")
        
        with self.get_session() as session:
            stored_count = 0
            
            for doc in documents:
                # Check if document already exists
                existing = session.query(MedicalDocumentORM).filter_by(id=doc.id).first()
                
                if existing:
                    logger.debug(f"Document {doc.id} already exists, skipping")
                    continue
                
                # Create new document record
                doc_orm = MedicalDocumentORM(
                    id=doc.id,
                    text=doc.text,
                    specialty=doc.specialty,
                    source=doc.source,
                    confidence=doc.confidence,
                    keywords=doc.keywords,
                    metadata=doc.metadata,
                    created_at=doc.created_at or datetime.utcnow()
                )
                
                session.add(doc_orm)
                stored_count += 1
            
            session.commit()
            logger.info(f"Stored {stored_count} new documents")
            
        return stored_count
    
    def get_documents(self, 
                     specialty: Optional[str] = None,
                     source: Optional[str] = None,
                     limit: Optional[int] = None,
                     offset: int = 0) -> List[MedicalDocument]:
        """
        Retrieve medical documents from database.
        
        Args:
            specialty: Filter by specialty
            source: Filter by source
            limit: Maximum number of documents
            offset: Number of documents to skip
            
        Returns:
            List of medical documents
        """
        with self.get_session() as session:
            query = session.query(MedicalDocumentORM)
            
            if specialty:
                query = query.filter_by(specialty=specialty)
            if source:
                query = query.filter_by(source=source)
            
            if offset:
                query = query.offset(offset)
            if limit:
                query = query.limit(limit)
            
            doc_orms = query.all()
            
            # Convert to MedicalDocument objects
            documents = []
            for doc_orm in doc_orms:
                doc = MedicalDocument(
                    id=doc_orm.id,
                    text=doc_orm.text,
                    specialty=doc_orm.specialty,
                    source=doc_orm.source,
                    confidence=doc_orm.confidence,
                    keywords=doc_orm.keywords,
                    metadata=doc_orm.metadata,
                    created_at=doc_orm.created_at
                )
                documents.append(doc)
            
            return documents
    
    def store_features(self, features: ProcessedFeatures, feature_version: str = "v1.0") -> int:
        """
        Store extracted features in database.
        
        Args:
            features: Processed features object
            feature_version: Version identifier for features
            
        Returns:
            Number of feature records stored
        """
        logger.info(f"Storing features for {len(features.document_ids)} documents")
        
        with self.get_session() as session:
            stored_count = 0
            
            for i, doc_id in enumerate(features.document_ids):
                # Check if features already exist for this version
                existing = session.query(DocumentFeatureORM).filter_by(
                    document_id=doc_id, 
                    feature_version=feature_version
                ).first()
                
                if existing:
                    logger.debug(f"Features for document {doc_id} v{feature_version} already exist")
                    continue
                
                # Extract sparse TF-IDF features (store indices and values)
                doc_features = features.text_features[i]
                if hasattr(doc_features, 'toarray'):
                    # Sparse matrix - store non-zero elements
                    nonzero_indices = doc_features.nonzero()[1]
                    nonzero_values = doc_features.data
                    tfidf_data = {
                        'indices': nonzero_indices.tolist(),
                        'values': nonzero_values.tolist(),
                        'shape': doc_features.shape
                    }
                else:
                    # Dense array
                    tfidf_data = doc_features.tolist()
                
                # Create feature record
                feature_orm = DocumentFeatureORM(
                    document_id=doc_id,
                    feature_version=feature_version,
                    tfidf_features=tfidf_data
                )
                
                session.add(feature_orm)
                stored_count += 1
            
            session.commit()
            logger.info(f"Stored features for {stored_count} documents")
            
        return stored_count
    
    def get_dataset_statistics(self) -> Dict[str, Any]:
        """Get comprehensive dataset statistics."""
        
        with self.get_session() as session:
            # Total documents
            total_docs = session.query(MedicalDocumentORM).count()
            
            # Specialty distribution
            specialty_counts = {}
            for specialty in settings.model.medical_specialties:
                count = session.query(MedicalDocumentORM).filter_by(specialty=specialty).count()
                specialty_counts[specialty] = count
            
            # Source distribution
            source_results = session.query(
                MedicalDocumentORM.source,
                session.query().count().label('count')
            ).group_by(MedicalDocumentORM.source).all()
            
            source_counts = {source: count for source, count in source_results}
            
            # Recent document counts
            from sqlalchemy import func
            recent_count = session.query(MedicalDocumentORM).filter(
                MedicalDocumentORM.created_at >= func.now() - func.interval('7 days')
            ).count()
            
            return {
                'total_documents': total_docs,
                'specialty_distribution': specialty_counts,
                'source_distribution': source_counts,
                'recent_documents': recent_count,
                'total_features': session.query(DocumentFeatureORM).count(),
                'total_models': session.query(ModelORM).count(),
                'total_predictions': session.query(PredictionORM).count()
            }


class FeatureStore:
    """
    Manages feature storage and retrieval for ML pipelines.
    
    Provides versioned feature storage with efficient retrieval
    for training and inference pipelines.
    """
    
    def __init__(self, database_manager: DatabaseManager):
        """
        Initialize feature store.
        
        Args:
            database_manager: Database manager instance
        """
        self.db_manager = database_manager
        
    def store_features(self, features: ProcessedFeatures, version: str = "v1.0") -> bool:
        """Store processed features with version control."""
        
        try:
            self.db_manager.store_features(features, version)
            return True
        except Exception as e:
            logger.error(f"Error storing features: {e}")
            return False
    
    def get_features(self, 
                    document_ids: Optional[List[str]] = None,
                    version: str = "v1.0",
                    specialty: Optional[str] = None) -> Optional[ProcessedFeatures]:
        """
        Retrieve features for specified documents.
        
        Args:
            document_ids: Specific document IDs to retrieve
            version: Feature version to retrieve
            specialty: Filter by specialty
            
        Returns:
            ProcessedFeatures object or None
        """
        with self.db_manager.get_session() as session:
            query = session.query(DocumentFeatureORM, MedicalDocumentORM).join(
                MedicalDocumentORM, DocumentFeatureORM.document_id == MedicalDocumentORM.id
            ).filter(DocumentFeatureORM.feature_version == version)
            
            if document_ids:
                query = query.filter(DocumentFeatureORM.document_id.in_(document_ids))
            if specialty:
                query = query.filter(MedicalDocumentORM.specialty == specialty)
            
            results = query.all()
            
            if not results:
                return None
            
            # Reconstruct features
            feature_data = []
            labels = []
            doc_ids = []
            
            for feature_orm, doc_orm in results:
                # Reconstruct TF-IDF features
                tfidf_data = feature_orm.tfidf_features
                if isinstance(tfidf_data, dict) and 'indices' in tfidf_data:
                    # Sparse representation - would need scipy to reconstruct
                    # For now, store as dense
                    feature_data.append(tfidf_data['values'])
                else:
                    feature_data.append(tfidf_data)
                
                # Get label encoding
                labels.append(doc_orm.specialty)
                doc_ids.append(doc_orm.id)
            
            return ProcessedFeatures(
                text_features=np.array(feature_data),
                text_metadata={'retrieved_from_db': True},
                labels=np.array(labels),
                feature_names=[],  # Would need to store separately
                document_ids=doc_ids
            )
    
    def get_feature_versions(self) -> List[str]:
        """Get list of available feature versions."""
        
        with self.db_manager.get_session() as session:
            versions = session.query(DocumentFeatureORM.feature_version).distinct().all()
            return [v[0] for v in versions]


class DataWarehouse:
    """
    Manages data warehousing operations for analytics and reporting.
    
    Provides aggregated views and analytics for the medical
    classification system.
    """
    
    def __init__(self, database_manager: DatabaseManager):
        """
        Initialize data warehouse.
        
        Args:
            database_manager: Database manager instance
        """
        self.db_manager = database_manager
    
    def create_specialty_summary(self) -> pd.DataFrame:
        """Create specialty summary statistics."""
        
        with self.db_manager.get_session() as session:
            # Query document counts by specialty
            from sqlalchemy import func
            
            results = session.query(
                MedicalDocumentORM.specialty,
                func.count(MedicalDocumentORM.id).label('document_count'),
                func.avg(func.length(MedicalDocumentORM.text)).label('avg_text_length'),
                func.min(MedicalDocumentORM.created_at).label('first_document'),
                func.max(MedicalDocumentORM.created_at).label('latest_document')
            ).group_by(MedicalDocumentORM.specialty).all()
            
            # Convert to DataFrame
            df = pd.DataFrame(results, columns=[
                'specialty', 'document_count', 'avg_text_length', 
                'first_document', 'latest_document'
            ])
            
            return df
    
    def create_source_analysis(self) -> pd.DataFrame:
        """Create data source analysis."""
        
        with self.db_manager.get_session() as session:
            from sqlalchemy import func
            
            results = session.query(
                MedicalDocumentORM.source,
                MedicalDocumentORM.specialty,
                func.count(MedicalDocumentORM.id).label('document_count')
            ).group_by(
                MedicalDocumentORM.source,
                MedicalDocumentORM.specialty
            ).all()
            
            df = pd.DataFrame(results, columns=['source', 'specialty', 'document_count'])
            
            # Pivot to get source-specialty matrix
            pivot_df = df.pivot(index='source', columns='specialty', values='document_count')
            pivot_df = pivot_df.fillna(0)
            
            return pivot_df
    
    def get_model_performance_trends(self) -> pd.DataFrame:
        """Get model performance trends over time."""
        
        with self.db_manager.get_session() as session:
            from sqlalchemy import func, case
            
            # Calculate daily accuracy trends
            results = session.query(
                func.date(PredictionORM.prediction_time).label('date'),
                PredictionORM.model_id,
                func.count(PredictionORM.id).label('total_predictions'),
                func.sum(case((PredictionORM.is_correct == True, 1), else_=0)).label('correct_predictions')
            ).filter(
                PredictionORM.is_correct.isnot(None)
            ).group_by(
                func.date(PredictionORM.prediction_time),
                PredictionORM.model_id
            ).all()
            
            df = pd.DataFrame(results, columns=[
                'date', 'model_id', 'total_predictions', 'correct_predictions'
            ])
            
            if not df.empty:
                df['accuracy'] = df['correct_predictions'] / df['total_predictions']
            
            return df

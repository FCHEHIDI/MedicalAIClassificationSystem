"""
Medical Data Pipeline Script
===========================

Complete data pipeline for ingesting, processing, and storing medical data
from multiple sources. This script demonstrates the full workflow from
raw data to ML-ready features.

Usage:
    python -m src.data.pipeline --source mimic --limit 1000
    python -m src.data.pipeline --source pubmed --email your@email.com
    python -m src.data.pipeline --source synthetic --count 500
    python -m src.data.pipeline --source all --output data/processed/medical_dataset.json
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional
import json
from datetime import datetime

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import settings
from utils.logging import setup_logging, get_logger
from data.ingestion import (
    MIMICDataLoader,
    PubMedDataLoader, 
    SyntheticDataGenerator,
    DataIngestionPipeline
)
from data.preprocessing import MedicalTextPreprocessor, FeatureExtractor, DataValidator
from data.storage import DatabaseManager, FeatureStore, DataWarehouse

logger = get_logger(__name__)


def create_mimic_loader(connection_string: Optional[str] = None) -> MIMICDataLoader:
    """Create MIMIC data loader with connection."""
    if not connection_string:
        # Use default database URL (you would need MIMIC-III access)
        connection_string = settings.database.url.replace('medical_classifier', 'mimic')
        logger.warning("Using default MIMIC connection. Ensure you have proper MIMIC-III access.")
    
    return MIMICDataLoader(connection_string)


def create_pubmed_loader(email: str, api_key: Optional[str] = None) -> PubMedDataLoader:
    """Create PubMed data loader."""
    return PubMedDataLoader(email=email, api_key=api_key)


def create_synthetic_loader(use_openai: bool = False, api_key: Optional[str] = None) -> SyntheticDataGenerator:
    """Create synthetic data generator."""
    return SyntheticDataGenerator(use_openai=use_openai, openai_api_key=api_key)


def run_data_pipeline(args) -> Dict[str, any]:
    """
    Run the complete data pipeline.
    
    Args:
        args: Command line arguments
        
    Returns:
        Pipeline statistics and results
    """
    logger.info("Starting Medical Data Pipeline")
    
    # Initialize components
    pipeline = DataIngestionPipeline(settings.database.url)
    preprocessor = MedicalTextPreprocessor(preserve_medical_terms=True)
    feature_extractor = FeatureExtractor(
        max_features=10000,
        ngram_range=(1, 3),
        min_df=2,
        max_df=0.8
    )
    validator = DataValidator()
    
    # Add data loaders based on source
    if args.source in ['mimic', 'all']:
        logger.info("Adding MIMIC-III data loader")
        try:
            mimic_loader = create_mimic_loader(args.mimic_connection)
            pipeline.add_loader(mimic_loader)
        except Exception as e:
            logger.error(f"Failed to add MIMIC loader: {e}")
            if args.source == 'mimic':
                return {'error': f'MIMIC loader failed: {e}'}
    
    if args.source in ['pubmed', 'all']:
        logger.info("Adding PubMed data loader")
        if not args.email:
            logger.error("Email required for PubMed access")
            if args.source == 'pubmed':
                return {'error': 'Email required for PubMed'}
        else:
            pubmed_loader = create_pubmed_loader(args.email, args.pubmed_api_key)
            pipeline.add_loader(pubmed_loader)
    
    if args.source in ['synthetic', 'all']:
        logger.info("Adding synthetic data generator")
        synthetic_loader = create_synthetic_loader(
            use_openai=args.use_openai,
            api_key=args.openai_api_key
        )
        pipeline.add_loader(synthetic_loader)
    
    # Run ingestion pipeline
    logger.info("Running data ingestion pipeline")
    ingestion_stats = pipeline.run_pipeline(
        output_path=Path(args.output) if args.output else None,
        limit_per_source=args.limit
    )
    
    logger.info(f"Ingestion complete: {ingestion_stats}")
    
    # Load documents for processing
    db_manager = DatabaseManager()
    documents = db_manager.get_documents(limit=args.limit * 3 if args.limit else None)
    
    if not documents:
        logger.warning("No documents found for processing")
        return {
            'ingestion_stats': ingestion_stats,
            'processing_stats': {},
            'error': 'No documents to process'
        }
    
    logger.info(f"Processing {len(documents)} documents")
    
    # Validate documents
    logger.info("Validating document quality")
    validation_results = validator.validate_documents(documents)
    
    logger.info(f"Validation results: {validation_results['quality_score']:.2f} quality score")
    
    if validation_results['valid_documents'] == 0:
        logger.error("No valid documents found")
        return {
            'ingestion_stats': ingestion_stats,
            'validation_results': validation_results,
            'error': 'No valid documents'
        }
    
    # Preprocess documents
    logger.info("Preprocessing medical text")
    processed_docs = []
    
    for doc in documents[:validation_results['valid_documents']]:
        try:
            processed_doc = preprocessor.preprocess_document(doc)
            processed_docs.append(processed_doc)
        except Exception as e:
            logger.error(f"Error preprocessing document {doc.id}: {e}")
    
    if not processed_docs:
        logger.error("No documents successfully preprocessed")
        return {
            'ingestion_stats': ingestion_stats,
            'validation_results': validation_results,
            'error': 'Preprocessing failed'
        }
    
    # Extract features
    logger.info("Extracting ML features")
    try:
        features = feature_extractor.fit_transform(processed_docs)
        logger.info(f"Extracted {features.text_features.shape[1]} features from {features.text_features.shape[0]} documents")
    except Exception as e:
        logger.error(f"Feature extraction failed: {e}")
        return {
            'ingestion_stats': ingestion_stats,
            'validation_results': validation_results,
            'error': f'Feature extraction failed: {e}'
        }
    
    # Store features if database available
    if args.store_features:
        logger.info("Storing features to database")
        try:
            db_manager.create_tables()  # Ensure tables exist
            feature_store = FeatureStore(db_manager)
            feature_store.store_features(features, version=args.feature_version)
            logger.info("Features stored successfully")
        except Exception as e:
            logger.error(f"Feature storage failed: {e}")
    
    # Generate summary statistics
    logger.info("Generating pipeline statistics")
    processing_stats = {
        'total_documents_processed': len(processed_docs),
        'features_extracted': features.text_features.shape[1],
        'specialty_distribution': validation_results['statistics'].get('specialty_distribution', {}),
        'average_text_length': validation_results['statistics'].get('text_length_stats', {}).get('mean', 0),
        'feature_extraction_time': datetime.now().isoformat(),
        'feature_version': args.feature_version
    }
    
    # Save processing results
    if args.output:
        output_path = Path(args.output)
        results_path = output_path.parent / f"{output_path.stem}_results.json"
        
        pipeline_results = {
            'pipeline_config': {
                'source': args.source,
                'limit': args.limit,
                'feature_version': args.feature_version,
                'timestamp': datetime.now().isoformat()
            },
            'ingestion_stats': ingestion_stats,
            'validation_results': validation_results,
            'processing_stats': processing_stats,
            'feature_metadata': features.text_metadata
        }
        
        with open(results_path, 'w') as f:
            json.dump(pipeline_results, f, indent=2, default=str)
        
        logger.info(f"Pipeline results saved to {results_path}")
    
    logger.info("Medical Data Pipeline completed successfully!")
    
    return {
        'ingestion_stats': ingestion_stats,
        'validation_results': validation_results, 
        'processing_stats': processing_stats,
        'features': features
    }


def main():
    """Main pipeline execution."""
    
    parser = argparse.ArgumentParser(description='Medical Data Pipeline')
    
    # Data source options
    parser.add_argument('--source', 
                       choices=['mimic', 'pubmed', 'synthetic', 'all'],
                       default='synthetic',
                       help='Data source to process')
    
    # General options
    parser.add_argument('--limit', type=int, default=100,
                       help='Limit documents per source')
    
    parser.add_argument('--output', type=str,
                       help='Output file for processed data')
    
    parser.add_argument('--feature-version', type=str, default='v1.0',
                       help='Feature version identifier')
    
    parser.add_argument('--store-features', action='store_true',
                       help='Store features to database')
    
    # MIMIC options
    parser.add_argument('--mimic-connection', type=str,
                       help='MIMIC database connection string')
    
    # PubMed options  
    parser.add_argument('--email', type=str,
                       help='Email for PubMed API access')
    
    parser.add_argument('--pubmed-api-key', type=str,
                       help='PubMed API key for higher limits')
    
    # Synthetic data options
    parser.add_argument('--use-openai', action='store_true',
                       help='Use OpenAI for synthetic data generation')
    
    parser.add_argument('--openai-api-key', type=str,
                       help='OpenAI API key')
    
    # Logging options
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(level=args.log_level)
    
    # Validate arguments
    if args.source == 'pubmed' and not args.email:
        parser.error("--email required when using PubMed source")
    
    if args.use_openai and not args.openai_api_key:
        parser.error("--openai-api-key required when using OpenAI")
    
    try:
        # Run pipeline
        results = run_data_pipeline(args)
        
        if 'error' in results:
            logger.error(f"Pipeline failed: {results['error']}")
            sys.exit(1)
        
        # Print summary
        print("\n" + "="*60)
        print("MEDICAL DATA PIPELINE SUMMARY")
        print("="*60)
        
        if 'ingestion_stats' in results:
            stats = results['ingestion_stats']
            print(f"📊 Total Documents: {stats['total_documents']}")
            print(f"📚 By Specialty:")
            for specialty, count in stats['by_specialty'].items():
                print(f"   {specialty}: {count}")
            print(f"🏥 By Source:")
            for source, count in stats['by_source'].items():
                print(f"   {source}: {count}")
        
        if 'validation_results' in results:
            val_stats = results['validation_results']
            print(f"✅ Valid Documents: {val_stats['valid_documents']}/{val_stats['total_documents']}")
            print(f"🎯 Quality Score: {val_stats['quality_score']:.2%}")
        
        if 'processing_stats' in results:
            proc_stats = results['processing_stats'] 
            print(f"🔧 Features Extracted: {proc_stats['features_extracted']}")
            print(f"📝 Average Text Length: {proc_stats['average_text_length']:.0f}")
        
        print("="*60)
        print("✅ Pipeline completed successfully!")
        
    except KeyboardInterrupt:
        logger.info("Pipeline interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed with error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

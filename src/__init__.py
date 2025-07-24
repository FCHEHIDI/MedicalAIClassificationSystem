"""
Medical Classification Engine
============================

A production-ready medical document classification system with MLOps capabilities.

This package provides:
- Medical text classification with 5 clinical specialties
- MLflow integration for experiment tracking
- FastAPI REST API for model serving
- Streamlit dashboard for interactive use
- PostgreSQL integration for data management
- Comprehensive monitoring and logging

Modules:
--------
- api: FastAPI REST API endpoints
- data: Data processing and feature engineering
- models: Machine learning models and training
- dashboard: Streamlit web application
- utils: Common utilities and helpers
- config: Configuration management
"""

__version__ = "1.0.0"
__author__ = "Fares"
__email__ = "fareschehidi7@gmail.com"

# Medical specialties supported by the system
MEDICAL_SPECIALTIES = [
    "Cardiology",
    "Emergency",
    "Pulmonology", 
    "Gastroenterology",
    "Dermatology"
]

# Core exports
from .config import settings
from .utils.logging import setup_logging

__all__ = [
    "MEDICAL_SPECIALTIES",
    "settings",
    "setup_logging",
]

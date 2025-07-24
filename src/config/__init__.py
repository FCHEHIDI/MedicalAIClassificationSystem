"""
Configuration Management
========================

Centralized configuration using Pydantic settings with environment variable support.
Follows 12-factor app principles for configuration management.
"""

import os
from typing import List, Optional

from pydantic import Field
from pydantic_settings import BaseSettings


class DatabaseSettings(BaseSettings):
    """Database configuration settings."""
    
    host: str = Field(default="localhost", description="Database host")
    port: int = Field(default=5432, description="Database port")
    name: str = Field(default="medical_classifier", description="Database name")
    user: str = Field(default="meduser", description="Database user")
    password: str = Field(default="medpass123", description="Database password")
    
    @property
    def url(self) -> str:
        """Construct database URL."""
        return f"postgresql://{self.user}:{self.password}@{self.host}:{self.port}/{self.name}"
    
    class Config:
        env_prefix = "DATABASE_"


class MLFlowSettings(BaseSettings):
    """MLflow configuration settings."""
    
    tracking_uri: str = Field(default="http://localhost:5000", description="MLflow tracking URI")
    experiment_name: str = Field(default="medical_classification", description="Default experiment name")
    s3_endpoint_url: Optional[str] = Field(default=None, description="S3 endpoint for artifacts")
    
    class Config:
        env_prefix = "MLFLOW_"


class APISettings(BaseSettings):
    """API configuration settings."""
    
    host: str = Field(default="0.0.0.0", description="API host")
    port: int = Field(default=8000, description="API port")
    debug: bool = Field(default=False, description="Debug mode")
    secret_key: str = Field(default="change-me-in-production", description="JWT secret key")
    algorithm: str = Field(default="HS256", description="JWT algorithm")
    access_token_expire_minutes: int = Field(default=30, description="Token expiration time")
    
    class Config:
        env_prefix = "API_"


class ModelSettings(BaseSettings):
    """Model configuration settings."""
    
    default_model_name: str = Field(default="medical_classifier_v1", description="Default model name")
    registry_path: str = Field(default="models/", description="Model registry path")
    max_sequence_length: int = Field(default=512, description="Maximum input sequence length")
    batch_size: int = Field(default=32, description="Default batch size")
    
    # Medical specialties - ensures consistency across the system
    medical_specialties: List[str] = Field(
        default=[
            "Cardiology",
            "Emergency", 
            "Pulmonology",
            "Gastroenterology",
            "Dermatology"
        ],
        description="Supported medical specialties"
    )
    
    class Config:
        env_prefix = "MODEL_"


class DataSettings(BaseSettings):
    """Data processing configuration settings."""
    
    data_path: str = Field(default="data/", description="Base data directory")
    raw_data_path: str = Field(default="data/raw/", description="Raw data directory")
    processed_data_path: str = Field(default="data/processed/", description="Processed data directory")
    features_path: str = Field(default="data/features/", description="Features directory")
    
    class Config:
        env_prefix = "DATA_"


class LoggingSettings(BaseSettings):
    """Logging configuration settings."""
    
    level: str = Field(default="INFO", description="Log level")
    format: str = Field(default="json", description="Log format (json/text)")
    file: str = Field(default="logs/app.log", description="Log file path")
    
    class Config:
        env_prefix = "LOG_"


class RedisSettings(BaseSettings):
    """Redis configuration settings."""
    
    host: str = Field(default="localhost", description="Redis host")
    port: int = Field(default=6379, description="Redis port")
    password: Optional[str] = Field(default=None, description="Redis password")
    
    @property
    def url(self) -> str:
        """Construct Redis URL."""
        if self.password:
            return f"redis://:{self.password}@{self.host}:{self.port}/0"
        return f"redis://{self.host}:{self.port}/0"
    
    class Config:
        env_prefix = "REDIS_"


class Settings(BaseSettings):
    """Main application settings."""
    
    # Environment
    environment: str = Field(default="development", description="Environment name")
    debug: bool = Field(default=True, description="Debug mode")
    testing: bool = Field(default=False, description="Testing mode")
    
    # Component settings
    database: DatabaseSettings = Field(default_factory=DatabaseSettings)
    mlflow: MLFlowSettings = Field(default_factory=MLFlowSettings)
    api: APISettings = Field(default_factory=APISettings)
    model: ModelSettings = Field(default_factory=ModelSettings)
    data: DataSettings = Field(default_factory=DataSettings)
    logging: LoggingSettings = Field(default_factory=LoggingSettings)
    redis: RedisSettings = Field(default_factory=RedisSettings)
    
    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"
        case_sensitive = False


# Global settings instance
settings = Settings()


def get_settings() -> Settings:
    """Get application settings."""
    return settings

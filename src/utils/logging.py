"""
Logging Utilities
================

Centralized logging configuration with structured logging support.
Provides JSON formatting for production and human-readable format for development.
"""

import json
import logging
import logging.config
import sys
from pathlib import Path
from typing import Any, Dict

import structlog
from rich.logging import RichHandler

from ..config import settings


class JSONFormatter(logging.Formatter):
    """Custom JSON formatter for structured logging."""
    
    def format(self, record: logging.LogRecord) -> str:
        """Format log record as JSON."""
        log_entry = {
            "timestamp": self.formatTime(record),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        
        # Add extra fields
        if hasattr(record, "user_id"):
            log_entry["user_id"] = record.user_id
        if hasattr(record, "request_id"):
            log_entry["request_id"] = record.request_id
        if hasattr(record, "model_version"):
            log_entry["model_version"] = record.model_version
            
        # Add exception info if present
        if record.exc_info:
            log_entry["exception"] = self.formatException(record.exc_info)
            
        return json.dumps(log_entry)


def setup_logging() -> None:
    """Configure application logging."""
    
    # Ensure log directory exists
    log_file_path = Path(settings.logging.file)
    log_file_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Configure based on format preference
    if settings.logging.format.lower() == "json":
        # Production-style JSON logging
        logging_config = {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "json": {
                    "()": "src.utils.logging.JSONFormatter",
                },
                "console": {
                    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                },
            },
            "handlers": {
                "file": {
                    "level": settings.logging.level,
                    "class": "logging.handlers.RotatingFileHandler",
                    "filename": settings.logging.file,
                    "maxBytes": 10485760,  # 10MB
                    "backupCount": 5,
                    "formatter": "json",
                },
                "console": {
                    "level": settings.logging.level,
                    "class": "logging.StreamHandler",
                    "stream": sys.stdout,
                    "formatter": "console",
                },
            },
            "loggers": {
                "": {  # root logger
                    "level": settings.logging.level,
                    "handlers": ["file", "console"],
                },
                "uvicorn": {
                    "level": "INFO",
                    "handlers": ["file", "console"],
                    "propagate": False,
                },
                "sqlalchemy": {
                    "level": "WARNING",
                    "handlers": ["file"],
                    "propagate": False,
                },
            },
        }
    else:
        # Development-style rich logging
        logging_config = {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {
                "rich": {
                    "format": "%(message)s",
                },
                "file": {
                    "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
                },
            },
            "handlers": {
                "rich": {
                    "level": settings.logging.level,
                    "class": "rich.logging.RichHandler",
                    "formatter": "rich",
                    "rich_tracebacks": True,
                    "markup": True,
                },
                "file": {
                    "level": settings.logging.level,
                    "class": "logging.handlers.RotatingFileHandler",
                    "filename": settings.logging.file,
                    "maxBytes": 10485760,  # 10MB
                    "backupCount": 5,
                    "formatter": "file",
                },
            },
            "loggers": {
                "": {  # root logger
                    "level": settings.logging.level,
                    "handlers": ["rich", "file"],
                },
            },
        }
    
    # Apply configuration
    logging.config.dictConfig(logging_config)
    
    # Configure structlog for structured logging
    structlog.configure(
        processors=[
            structlog.stdlib.filter_by_level,
            structlog.stdlib.add_logger_name,
            structlog.stdlib.add_log_level,
            structlog.stdlib.PositionalArgumentsFormatter(),
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            structlog.processors.JSONRenderer() if settings.logging.format.lower() == "json" 
            else structlog.processors.KeyValueRenderer(key_order=['timestamp', 'level', 'event']),
        ],
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        wrapper_class=structlog.stdlib.BoundLogger,
        cache_logger_on_first_use=True,
    )


def get_logger(name: str) -> logging.Logger:
    """Get a configured logger instance."""
    return logging.getLogger(name)


def get_structured_logger(name: str) -> Any:
    """Get a structured logger instance."""
    return structlog.get_logger(name)

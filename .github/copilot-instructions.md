<!-- Use this file to provide workspace-specific custom instructions to Copilot. For more details, visit https://code.visualstudio.com/docs/copilot/copilot-customization#_use-a-githubcopilotinstructionsmd-file -->

# Medical Document Classifier - Copilot Instructions

## Project Context
This is a professional medical document classification system built with MLOps best practices. The system classifies medical documents into 5 specialties: Cardiology, Emergency, Pulmonology, Gastroenterology, and Dermatology.

## Architecture & Technologies
- **Backend**: FastAPI for model serving and APIs
- **Frontend**: Streamlit for medical professional dashboard
- **ML Stack**: scikit-learn, MLflow, feature engineering
- **Database**: PostgreSQL for data warehouse and feature store
- **Deployment**: Docker containers with docker-compose
- **Monitoring**: Custom metrics and logging system

## Code Style & Standards
- Follow PEP 8 Python style guide
- Use type hints consistently
- Implement proper error handling and logging
- Write docstrings for all functions and classes
- Include unit tests for critical functionality
- Use dependency injection patterns for testability

## Medical Domain Considerations
- Handle medical terminology and abbreviations correctly
- Maintain HIPAA compliance principles (no real patient data)
- Use appropriate medical classification confidence thresholds
- Implement proper medical data validation
- Consider clinical workflow patterns in UI design

## MLOps Best Practices
- Version all datasets and models
- Track experiments with MLflow
- Implement proper train/validation/test splits
- Monitor for data drift and model performance
- Use feature stores for consistent ML features
- Implement proper model deployment patterns

## Key Patterns to Follow
- Repository pattern for data access
- Factory pattern for model creation
- Strategy pattern for different classification approaches
- Observer pattern for monitoring and alerts
- Config-driven development for flexibility

## Security & Compliance
- No real patient data in code or commits
- Use environment variables for sensitive config
- Implement proper API authentication
- Follow data privacy best practices
- Sanitize all user inputs

When suggesting code, prioritize production-ready, maintainable solutions that follow healthcare industry standards.

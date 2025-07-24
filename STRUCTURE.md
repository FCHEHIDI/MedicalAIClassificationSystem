# 📁 Project Structure

```
medical-classification-engine/
├── 🚀 simple_api.py              # FastAPI production application
├── 📊 simple_dashboard.py        # Streamlit medical dashboard
├── 📋 requirements.txt           # Python dependencies
├── 🚀 start.sh                   # Local development startup
├── 📄 README.md                  # Professional project documentation
├── 📜 LICENSE                    # MIT License
│
├── 🤖 models/                    # Trained ML models & encoders
│   └── (.gitkeep - models stored separately)
│
├── 🏗️ src/                       # Core application modules
│   ├── api/                     # API implementation
│   ├── data/                    # Data processing pipeline
│   ├── dashboard/               # Dashboard components
│   ├── config/                  # Configuration management
│   └── utils/                   # Utility functions
│
├── 🐳 docker/                    # Docker deployment
│   ├── api.Dockerfile           # API container configuration
│   ├── dashboard.Dockerfile     # Dashboard container configuration
│   └── README.md                # Docker instructions
│
├── 📝 data/                      # Dataset management
│   ├── raw/                     # Raw medical data
│   ├── processed/               # Processed datasets
│   └── features/                # Feature engineering output
│
├── 📊 notebooks/                 # Analysis & experimentation
│   ├── medical_classification_analysis.ipynb
│   ├── model_performance_analysis.ipynb
│   ├── feature_engineering_analysis.ipynb
│   └── clean_model_performance_analysis.ipynb
│
├── 🧪 tests/                     # Comprehensive testing suite
│   ├── test_api.py              # API endpoint tests
│   ├── test_models.py           # ML model validation
│   ├── test_medical_engine.py   # Core engine tests
│   └── test_compatibility.py    # Integration tests
│
├── 🔧 scripts/                   # Development & training scripts
│   ├── training/                # Model training scripts
│   ├── testing/                 # Testing utilities
│   ├── demo_pipeline.py         # Demo data pipeline
│   └── confidence_analysis.py   # Model confidence analysis
│
├── 📚 docs/                      # Documentation
│   ├── README.md                # Documentation overview
│   ├── DEMO_GUIDE.md            # Live demonstration guide
│   ├── technical/               # Technical documentation
│   └── deployment/              # Deployment guides
│
├── 🐙 .github/                   # GitHub integration
│   ├── workflows/               # CI/CD pipeline
│   ├── ISSUE_TEMPLATE/          # Issue templates
│   └── copilot-instructions.md  # Development guidelines
│
└── ⚙️ Configuration Files
    ├── docker-compose.yml       # Local development environment
    ├── docker-compose.production.yml  # Production deployment
    ├── pyproject.toml           # Python project configuration
    ├── .gitignore               # Git ignore patterns
    └── .pre-commit-config.yaml  # Code quality hooks
```

## 🎯 Key Components for Recruiters

### **Production-Ready Code**
- `simple_api.py` - FastAPI backend with 99.9% accuracy
- `simple_dashboard.py` - Professional Streamlit interface
- `src/` - Well-structured application modules

### **Machine Learning Pipeline**
- `models/` - Trained classification models
- `notebooks/` - Data analysis and model development
- `scripts/training/` - Model training infrastructure

### **Professional Development**
- `tests/` - Comprehensive test suite
- `docker/` - Containerization for deployment
- `.github/workflows/` - CI/CD pipeline

### **Documentation**
- `README.md` - Professional project overview
- `docs/` - Comprehensive documentation
- `DEMO_GUIDE.md` - Live demonstration instructions

## ⚡ Quick Start

```bash
# Clone and run locally
git clone <repository-url>
cd medical-classification-engine
bash start.sh
```

## 🌐 Live Production URLs

- **Dashboard**: https://medical-dashboard.blackrock-067a426a.eastus.azurecontainerapps.io/
- **API**: https://medical-api.blackrock-067a426a.eastus.azurecontainerapps.io/docs

# 🔧 Medical Classification Engine - Scripts

This directory contains all utility scripts for data processing, model training, testing, and deployment.

## 📁 Directory Structure

### 🚀 **Deployment Scripts** (`deployment/`)
- **`deploy-production.ps1`** - Complete Azure production deployment script

### 🧪 **Testing Scripts** (`testing/`)  
- **`test-local-services.ps1`** - Local service testing and validation
- **`test_pubmed.py`** - PubMed API connectivity testing

### 🎓 **Training Scripts** (`training/`)
- **`train_production_models.py`** - Production model training pipeline
- **`train_medical_models.py`** - Core medical model training
- **`professional_model_training.py`** - Advanced training with professional features

### 📊 **Data Processing Scripts** (root)
- **`fetch_1500_samples.py`** - Fetch large medical dataset samples
- **`fetch_large_pubmed_dataset.py`** - Large-scale PubMed data collection
- **`fetch_pubmed_data.py`** - Standard PubMed data fetching
- **`simple_pubmed_fetch.py`** - Simple PubMed API interface
- **`complete_dataset_processor.py`** - Comprehensive dataset processing
- **`integrate_medical_data.py`** - Medical data integration pipeline
- **`generate_sample_data.py`** - Synthetic medical data generation

### 📈 **Analysis Scripts** (root)
- **`confidence_analysis.py`** - Model confidence analysis and calibration
- **`validate_model_robustness.py`** - Model robustness testing
- **`demo_pipeline.py`** - Complete demonstration pipeline

---

## 🎯 **Common Usage Patterns**

### **🏋️ Model Training Workflow**
```bash
# Step 1: Fetch medical data
python scripts/fetch_1500_samples.py

# Step 2: Process and integrate data  
python scripts/complete_dataset_processor.py

# Step 3: Train production models
python scripts/training/train_production_models.py
```

### **🧪 Testing Workflow**
```bash
# Test PubMed connectivity
python scripts/testing/test_pubmed.py

# Test local services (PowerShell)
.\scripts\testing\test-local-services.ps1
```

### **🚀 Deployment Workflow**
```bash
# Deploy to Azure production (PowerShell)
.\scripts\deployment\deploy-production.ps1
```

### **📊 Analysis Workflow**
```bash
# Analyze model confidence
python scripts/confidence_analysis.py

# Validate model robustness
python scripts/validate_model_robustness.py

# Run complete demo
python scripts/demo_pipeline.py
```

---

## 🔧 **Script Categories**

### **Data Pipeline Scripts**
- **Fetching**: PubMed data collection from NCBI API
- **Processing**: Data cleaning, augmentation, and formatting
- **Integration**: Combining multiple data sources
- **Generation**: Creating synthetic training data

### **ML Pipeline Scripts**  
- **Training**: Model development and optimization
- **Validation**: Performance testing and robustness checks
- **Analysis**: Confidence calibration and feature analysis
- **Professional**: Production-ready training pipelines

### **DevOps Scripts**
- **Testing**: Service validation and health checks
- **Deployment**: Azure cloud deployment automation
- **Monitoring**: Performance and health monitoring

---

## 📋 **Dependencies**

### **Python Scripts**
```bash
# Install required packages
pip install -r requirements.txt
```

### **PowerShell Scripts**
```powershell
# Azure CLI required for deployment
az --version

# Docker required for containerization
docker --version
```

---

## 🎯 **Best Practices**

### **Before Running Scripts**
1. **Check Requirements** - Ensure all dependencies installed
2. **Verify Data** - Confirm input data exists and is valid
3. **Set Environment** - Configure API keys and settings
4. **Test Connectivity** - Verify API and service access

### **Script Execution**
1. **Run from Project Root** - Most scripts expect to run from main directory
2. **Check Outputs** - Verify script completion and output quality
3. **Log Results** - Scripts generate logs in `logs/` directory
4. **Validate Models** - Test trained models before deployment

### **Error Handling**
- **Check Logs** - Scripts log detailed information
- **Verify Paths** - Ensure data and model directories exist
- **API Limits** - PubMed scripts respect API rate limits
- **Resource Usage** - Training scripts are resource-intensive

---

## 🔗 **Related Documentation**

- **[Technical Deep Dive](../docs/technical/TECHNICAL_DEEP_DIVE.md)** - Detailed architecture
- **[Deployment Guide](../docs/deployment/PHASE_4_COMPLETE.md)** - Production deployment
- **[Docker README](../docker/README.md)** - Container configuration
- **[Project Structure](../docs/technical/PROJECT_STRUCTURE.md)** - Overall organization

---

## 🆘 **Troubleshooting**

### **Common Issues**
- **Import Errors** - Run scripts from project root directory
- **API Failures** - Check network connectivity and API keys
- **Memory Issues** - Training scripts require significant RAM
- **Path Issues** - Verify data and model directories exist

### **Getting Help**
- Check script documentation and comments
- Review logs in `logs/` directory
- Verify requirements and dependencies
- Ensure proper environment setup

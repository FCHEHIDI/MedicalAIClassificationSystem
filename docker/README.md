# 🐳 Docker Configuration

This directory contains all Docker-related files for the Medical Classification Engine.

## 📁 Directory Contents

### 🐋 Dockerfiles
- **`api.Dockerfile`** - FastAPI medical classification service
- **`dashboard.Dockerfile`** - Streamlit dashboard interface  
- **`mlflow.Dockerfile`** - MLflow experiment tracking (optional)
- **`streamlit.Dockerfile`** - Standard Streamlit configuration
- **`streamlit-compat.Dockerfile`** - Compatible Streamlit version
- **`Dockerfile.production`** - Production deployment container

### 🧪 Training Scripts
- **`docker_train_models.py`** - Docker-compatible model training script

## 🚀 Usage

### Quick Start with Docker Compose
```bash
# From project root directory
docker-compose up -d
```

### Individual Container Builds
```bash
# Build API service
docker build -f docker/api.Dockerfile -t medical-api .

# Build Dashboard service  
docker build -f docker/dashboard.Dockerfile -t medical-dashboard .

# Build Production container
docker build -f docker/Dockerfile.production -t medical-engine-prod .
```

### Training Models in Docker
```bash
# From project root
cd docker
python docker_train_models.py
```

## 🔧 Configuration

### Service Ports
- **API**: 8001 → 8000 (container)
- **Dashboard**: 8501 → 8501 (container)
- **MLflow**: 5000 → 5000 (container, if used)

### Volume Mounts
- `./models:/app/models:ro` - Model files (read-only)
- Health checks configured for all services

### Networks
- **medical_net** - Bridge network for service communication

## 📊 Health Checks

All services include health check endpoints:
- **API Health**: `http://localhost:8001/health`
- **Dashboard Health**: `http://localhost:8501/_stcore/health`

## 🏗️ Architecture

```
┌─────────────────┐    ┌─────────────────┐
│  Dashboard      │    │  API Service    │
│  (Streamlit)    │◄──►│  (FastAPI)      │
│  Port: 8501     │    │  Port: 8001     │
└─────────────────┘    └─────────────────┘
         │                       │
         └───────┬───────────────┘
                 │
         ┌─────────────────┐
         │  Shared Models  │
         │  Volume Mount   │
         └─────────────────┘
```

## 🔒 Security

- **Read-only model mounts** - Models mounted as read-only volumes
- **Network isolation** - Services communicate via private bridge network
- **Health monitoring** - Automatic service health checks
- **Restart policies** - Automatic restart unless stopped

## 📝 Development

To modify Docker configurations:

1. **Update Dockerfiles** in this directory
2. **Test locally** with `docker-compose up`
3. **Update references** in `../docker-compose.yml` if needed
4. **Document changes** in this README

## 🧪 Testing

```bash
# Test all services
../test-local-services.ps1

# Test individual service health
curl http://localhost:8001/health
curl http://localhost:8501/_stcore/health
```

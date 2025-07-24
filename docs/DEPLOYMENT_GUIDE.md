# 🚀 Azure Container Apps Deployment Guide

This guide contains the complete, tested deployment process for the Medical AI Classification System on Azure Container Apps.

## 📋 Overview

This deployment was successfully used to deploy the production system currently running at:
- **Dashboard**: https://medical-dashboard.blackrock-067a426a.eastus.azurecontainerapps.io/
- **API**: https://medical-api.blackrock-067a426a.eastus.azurecontainerapps.io/docs

## 🛠️ Prerequisites

### Required Software
- **Azure CLI** - [Install Guide](https://docs.microsoft.com/en-us/cli/azure/install-azure-cli)
- **Docker Desktop** - [Download](https://www.docker.com/products/docker-desktop/)
- **Git** - For repository management

### Azure Resources (One-time Setup)
```bash
# Create Resource Group
az group create --name medical-ai-rg --location eastus

# Create Container Registry
az acr create --resource-group medical-ai-rg --name medicalairegistry2025 --sku Basic

# Create Container Apps Environment
az containerapp env create \
  --name medical-ai-env \
  --resource-group medical-ai-rg \
  --location eastus
```

### Initial Container Apps Creation
```bash
# Create API Container App (one-time)
az containerapp create \
  --name medical-api \
  --resource-group medical-ai-rg \
  --environment medical-ai-env \
  --image mcr.microsoft.com/azuredocs/containerapps-helloworld:latest \
  --target-port 8000 \
  --ingress 'external' \
  --registry-server medicalairegistry2025.azurecr.io

# Create Dashboard Container App (one-time)
az containerapp create \
  --name medical-dashboard \
  --resource-group medical-ai-rg \
  --environment medical-ai-env \
  --image mcr.microsoft.com/azuredocs/containerapps-helloworld:latest \
  --target-port 8501 \
  --ingress 'external' \
  --registry-server medicalairegistry2025.azurecr.io
```

## 🚀 Deployment Process

### Option 1: Automated Deployment (Recommended)

#### For Linux/macOS:
```bash
chmod +x deploy-azure-production.sh
./deploy-azure-production.sh
```

#### For Windows:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
.\deploy-azure-production.ps1
```

### Option 2: Manual Step-by-Step Deployment

#### Step 1: Build Docker Images
```bash
# Build API image
docker build -t medical-api:v2 -f docker/api.Dockerfile .

# Build Dashboard image  
docker build -t medical-dashboard:v2 -f docker/dashboard.Dockerfile .
```

#### Step 2: Tag and Push to Registry
```bash
# Login to registry
az acr login --name medicalairegistry2025

# Tag images
docker tag medical-api:v2 medicalairegistry2025.azurecr.io/medical-api:v2
docker tag medical-dashboard:v2 medicalairegistry2025.azurecr.io/medical-dashboard:v2

# Push images
docker push medicalairegistry2025.azurecr.io/medical-api:v2
docker push medicalairegistry2025.azurecr.io/medical-dashboard:v2
```

#### Step 3: Create Container App Revisions
```bash
# Create API revision
az containerapp revision copy \
  --name medical-api \
  --resource-group medical-ai-rg \
  --revision-suffix "v2" \
  --image medicalairegistry2025.azurecr.io/medical-api:v2

# Create Dashboard revision
az containerapp revision copy \
  --name medical-dashboard \
  --resource-group medical-ai-rg \
  --revision-suffix "latest" \
  --image medicalairegistry2025.azurecr.io/medical-dashboard:v2
```

#### Step 4: Verify Deployment
```bash
# Get URLs
API_URL=$(az containerapp show --name medical-api --resource-group medical-ai-rg --query "properties.configuration.ingress.fqdn" --output tsv)
DASHBOARD_URL=$(az containerapp show --name medical-dashboard --resource-group medical-ai-rg --query "properties.configuration.ingress.fqdn" --output tsv)

# Test health endpoint
curl https://$API_URL/health
```

## 🔧 Key Features of This Deployment

### ✅ **Production-Ready Architecture**
- **Container Apps Environment**: Managed Kubernetes with auto-scaling
- **Azure Container Registry**: Private Docker image storage
- **Revision Management**: Blue-green deployments with zero downtime
- **External Ingress**: Public HTTPS endpoints with SSL certificates

### ✅ **Deployment Best Practices**
- **Versioned Images**: Using `v2` tags instead of `latest` for reliability
- **Explicit Revision Suffixes**: Avoids caching issues in Container Apps
- **Health Checks**: Automated verification of successful deployment
- **Error Handling**: Comprehensive error checking and rollback capabilities

### ✅ **Security & Compliance**
- **Private Registry**: Images stored in Azure Container Registry
- **HTTPS Only**: All endpoints secured with SSL/TLS
- **Managed Identity**: Secure authentication between services
- **Network Security**: Container Apps environment isolation

## 🚨 Common Issues & Solutions

### Issue 1: Container App Not Updating
**Problem**: New code doesn't appear after deployment
**Solution**: Use explicit revision suffixes instead of relying on 'latest' tags

```bash
# ❌ This can cause caching issues
--image registry.azurecr.io/app:latest

# ✅ This forces new revision
--revision-suffix "v2" --image registry.azurecr.io/app:v2
```

### Issue 2: Docker Build Fails
**Problem**: Docker build fails with missing files
**Solution**: Ensure you're in the project root directory and all files exist

```bash
# Check current directory
ls -la simple_api.py simple_dashboard.py requirements.txt docker/

# Verify Docker is running
docker info
```

### Issue 3: Azure Login Expired
**Problem**: Authentication errors during deployment
**Solution**: Re-authenticate with Azure

```bash
az login
az acr login --name medicalairegistry2025
```

### Issue 4: Health Check Fails
**Problem**: Health endpoint returns 404 or times out
**Solution**: Wait for container startup (30-60 seconds) and check logs

```bash
# Check container app logs
az containerapp logs show --name medical-api --resource-group medical-ai-rg --follow
```

## 📊 Performance & Scaling

### Current Configuration
- **CPU**: 0.25 vCPU per container
- **Memory**: 0.5 GB per container  
- **Scaling**: 1-10 replicas based on HTTP requests
- **Target Port**: API (8000), Dashboard (8501)

### Scaling Configuration
```bash
# Update scaling rules
az containerapp update \
  --name medical-api \
  --resource-group medical-ai-rg \
  --min-replicas 1 \
  --max-replicas 10 \
  --scale-rule-name "http-requests" \
  --scale-rule-type "http" \
  --scale-rule-metadata "concurrentRequests=10"
```

## 💡 Deployment Tips

### 1. **Version Your Images**
Always use specific version tags (v1, v2, v3) instead of 'latest' for predictable deployments.

### 2. **Test Locally First**
Test your Docker images locally before pushing to the registry:
```bash
docker run -p 8000:8000 medical-api:v2
docker run -p 8501:8501 medical-dashboard:v2
```

### 3. **Monitor Resource Usage**
Check container resource usage and adjust if needed:
```bash
az containerapp show --name medical-api --resource-group medical-ai-rg --query "properties.template.containers[0].resources"
```

### 4. **Use Staging Environment**
Consider creating a staging environment for testing before production deployment.

## 📈 Cost Optimization

### Current Monthly Cost (Estimated)
- **Container Apps**: ~$10-30/month (depending on usage)
- **Container Registry**: ~$5/month (Basic tier)
- **Data Transfer**: Minimal (within same region)

### Cost Reduction Tips
- Use scale-to-zero for development environments
- Optimize Docker image sizes
- Monitor and adjust resource allocations

## 🎯 Success Metrics

The deployment is considered successful when:
- ✅ Both applications are accessible via HTTPS
- ✅ API health endpoint returns 200 status
- ✅ Dashboard loads and can make predictions
- ✅ Auto-scaling works under load
- ✅ SSL certificates are valid

## 📞 Support

For deployment issues or questions:
- **Email**: fareschehidi7@gmail.com
- **GitHub Issues**: [Create an issue](https://github.com/FCHEHIDI/MedicalAIClassificationSystem/issues)
- **Documentation**: Check this guide and script comments

## 🔄 Updates and Maintenance

### Regular Updates
1. Update Docker images with new versions
2. Test in staging environment
3. Deploy using the same revision process
4. Monitor for any issues post-deployment

### Monitoring
- Set up Azure Monitor alerts
- Check container logs regularly
- Monitor resource usage and costs
- Verify SSL certificate renewals

---

*This deployment guide was created based on the successful production deployment of the Medical AI Classification System. Last updated: July 2025*

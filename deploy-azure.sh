#!/bin/bash
# Azure Medical AI Deployment Script

echo "🚀 Deploying Medical AI to Azure Container Instances"
echo "=================================================="

# Variables
RESOURCE_GROUP="medical-ai-rg"
REGISTRY_NAME="medicalairegistry2025"
LOCATION="eastus"

echo "📋 Deployment Configuration:"
echo "   Resource Group: $RESOURCE_GROUP"
echo "   Registry: $REGISTRY_NAME.azurecr.io"
echo "   Location: $LOCATION"
echo ""

# Deploy API
echo "🔥 Deploying Medical AI API..."
az container create \
  --resource-group $RESOURCE_GROUP \
  --name medical-api-prod \
  --image $REGISTRY_NAME.azurecr.io/medical-api:latest \
  --acr-identity $REGISTRY_NAME \
  --dns-name-label medical-api-prod-2025 \
  --ports 8000 \
  --cpu 2.0 \
  --memory 4.0 \
  --os-type Linux \
  --restart-policy Always \
  --environment-variables \
    ENV=production \
    API_HOST=0.0.0.0 \
    API_PORT=8000

echo "✅ API deployment initiated!"
echo "📡 API URL: http://medical-api-prod-2025.eastus.azurecontainer.io:8000"
echo ""

# Deploy Dashboard  
echo "🖥️ Deploying Medical AI Dashboard..."
az container create \
  --resource-group $RESOURCE_GROUP \
  --name medical-dashboard-prod \
  --image $REGISTRY_NAME.azurecr.io/medical-dashboard:latest \
  --registry-login-server $REGISTRY_NAME.azurecr.io \
  --dns-name-label medical-dashboard-prod-2025 \
  --ports 8501 \
  --cpu 1.0 \
  --memory 2.0 \
  --os-type Linux \
  --restart-policy Always \
  --environment-variables \
    STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    API_URL=http://medical-api-prod-2025.eastus.azurecontainer.io:8000

echo "✅ Dashboard deployment initiated!"
echo "🖥️ Dashboard URL: http://medical-dashboard-prod-2025.eastus.azurecontainer.io:8501"
echo ""

echo "🎉 DEPLOYMENT COMPLETE!"
echo "================================"
echo "📖 API Documentation: http://medical-api-prod-2025.eastus.azurecontainer.io:8000/docs"
echo "🖥️ Medical AI Dashboard: http://medical-dashboard-prod-2025.eastus.azurecontainer.io:8501"
echo "💰 Estimated Monthly Cost: ~$25-35"
echo ""
echo "🔍 Check deployment status:"
echo "   az container show --resource-group $RESOURCE_GROUP --name medical-api-prod --query 'provisioningState'"
echo "   az container show --resource-group $RESOURCE_GROUP --name medical-dashboard-prod --query 'provisioningState'"

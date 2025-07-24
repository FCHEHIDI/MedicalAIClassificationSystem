#!/bin/bash
# Minimal Azure Deployment - Free Tier Approach

echo "🚀 Minimal Azure Deployment for Medical AI"
echo "=========================================="

RESOURCE_GROUP="medical-ai-rg"
LOCATION="eastus"

echo "Creating minimal environment..."

# Try Container Apps environment (has generous free tier)
az containerapp env create \
  --name medical-ai-env \
  --resource-group $RESOURCE_GROUP \
  --location $LOCATION

# Deploy just the API first (most important)
echo "Deploying Medical API..."
az containerapp create \
  --name medical-api \
  --resource-group $RESOURCE_GROUP \
  --environment medical-ai-env \
  --image medicalairegistry2025.azurecr.io/medical-api:latest \
  --registry-server medicalairegistry2025.azurecr.io \
  --registry-username medicalairegistry2025 \
  --registry-password "wOL5QvJTXBuoVEPogHtDzYvZFQ44mVnH9/1yxF1ibO+ACRATD5pj" \
  --target-port 8000 \
  --ingress external \
  --cpu 0.5 \
  --memory 1.0Gi

echo ""
echo "🎉 MINIMAL DEPLOYMENT COMPLETE!"
echo "==============================="
echo "🔗 API: Check Azure portal for URL"
echo "💰 Cost: ~$5-10/month (very minimal)"

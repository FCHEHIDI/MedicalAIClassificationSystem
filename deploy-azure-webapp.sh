#!/bin/bash
# Simple Azure Deployment - Alternative Approach

echo "🚀 Alternative Azure Deployment for Medical AI"
echo "=============================================="

RESOURCE_GROUP="medical-ai-rg"
REGISTRY_NAME="medicalairegistry2025"
PLAN_NAME="medical-ai-plan"

echo "Creating App Service Plan..."
az appservice plan create \
  --name $PLAN_NAME \
  --resource-group $RESOURCE_GROUP \
  --sku B1 \
  --is-linux \
  --location eastus

echo "Creating Web App for API..."
az webapp create \
  --resource-group $RESOURCE_GROUP \
  --plan $PLAN_NAME \
  --name medical-api-webapp-2025 \
  --deployment-container-image-name $REGISTRY_NAME.azurecr.io/medical-api:latest

echo "Configuring registry credentials..."
az webapp config container set \
  --name medical-api-webapp-2025 \
  --resource-group $RESOURCE_GROUP \
  --container-image-name $REGISTRY_NAME.azurecr.io/medical-api:latest \
  --container-registry-url https://$REGISTRY_NAME.azurecr.io \
  --container-registry-user $REGISTRY_NAME \
  --container-registry-password "wOL5QvJTXBuoVEPogHtDzYvZFQ44mVnH9/1yxF1ibO+ACRATD5pj"

echo "Creating Web App for Dashboard..."
az webapp create \
  --resource-group $RESOURCE_GROUP \
  --plan $PLAN_NAME \
  --name medical-dashboard-webapp-2025 \
  --deployment-container-image-name $REGISTRY_NAME.azurecr.io/medical-dashboard:latest

echo "Configuring dashboard registry credentials..."
az webapp config container set \
  --name medical-dashboard-webapp-2025 \
  --resource-group $RESOURCE_GROUP \
  --container-image-name $REGISTRY_NAME.azurecr.io/medical-dashboard:latest \
  --container-registry-url https://$REGISTRY_NAME.azurecr.io \
  --container-registry-user $REGISTRY_NAME \
  --container-registry-password "wOL5QvJTXBuoVEPogHtDzYvZFQ44mVnH9/1yxF1ibO+ACRATD5pj"

echo ""
echo "🎉 DEPLOYMENT COMPLETE!"
echo "======================="
echo "🔗 API: https://medical-api-webapp-2025.azurewebsites.net"
echo "🔗 Dashboard: https://medical-dashboard-webapp-2025.azurewebsites.net"
echo "📖 API Docs: https://medical-api-webapp-2025.azurewebsites.net/docs"

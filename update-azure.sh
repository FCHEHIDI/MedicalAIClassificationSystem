#!/bin/bash
# Quick Azure Update Script - Updates deployed containers with latest code

echo "🚀 Updating Medical AI System on Azure Container Apps..."

# Configuration
RESOURCE_GROUP="medical-ai-rg"
REGISTRY_NAME="medicalairegistry"
API_APP="medical-api"
DASHBOARD_APP="medical-dashboard"

echo "📋 Step 1: Building updated images..."

# Build API image
docker build -f docker/api.Dockerfile -t $API_APP:latest .
if [ $? -eq 0 ]; then
    echo "✅ API image built successfully"
else
    echo "❌ API image build failed"
    exit 1
fi

# Build Dashboard image  
docker build -f docker/dashboard.Dockerfile -t $DASHBOARD_APP:latest .
if [ $? -eq 0 ]; then
    echo "✅ Dashboard image built successfully"
else
    echo "❌ Dashboard image build failed"
    exit 1
fi

echo "📤 Step 2: Pushing images to Azure Container Registry..."

# Login to ACR
az acr login --name $REGISTRY_NAME

# Tag and push API
docker tag $API_APP:latest $REGISTRY_NAME.azurecr.io/$API_APP:latest
docker push $REGISTRY_NAME.azurecr.io/$API_APP:latest

# Tag and push Dashboard
docker tag $DASHBOARD_APP:latest $REGISTRY_NAME.azurecr.io/$DASHBOARD_APP:latest
docker push $REGISTRY_NAME.azurecr.io/$DASHBOARD_APP:latest

echo "🔄 Step 3: Updating Azure Container Apps..."

# Update API container
az containerapp update \
    --name $API_APP \
    --resource-group $RESOURCE_GROUP \
    --image $REGISTRY_NAME.azurecr.io/$API_APP:latest

# Update Dashboard container
az containerapp update \
    --name $DASHBOARD_APP \
    --resource-group $RESOURCE_GROUP \
    --image $REGISTRY_NAME.azurecr.io/$DASHBOARD_APP:latest

echo "✅ Update completed! Your improved system is now live:"
echo "📊 Dashboard: https://medical-dashboard.blackrock-067a426a.eastus.azurecontainerapps.io/"
echo "🔧 API: https://medical-api.blackrock-067a426a.eastus.azurecontainerapps.io/docs"
echo ""
echo "🧪 Test the improvements:"
echo "- Try non-medical text to see 'Unknown/Low_Confidence' response"
echo "- Test medical cases to see enhanced confidence display"

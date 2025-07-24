#!/bin/bash
# Medical Classification Engine - Azure Container Apps Deployment (Final)
echo "🏥 Deploying Medical Classification Engine to Azure Container Apps"
echo "=================================================================="

# Variables
RESOURCE_GROUP="medical-ai-rg"
LOCATION="eastus"
REGISTRY_NAME="medicalairegistry2025"
ENVIRONMENT_NAME="medical-ai-env"

echo "📋 Getting registry credentials..."
REGISTRY_PASSWORD=$(az acr credential show --name $REGISTRY_NAME --query "passwords[0].value" --output tsv)

echo "🚀 Deploying Medical API..."
az containerapp create \
  --name medical-api \
  --resource-group $RESOURCE_GROUP \
  --environment $ENVIRONMENT_NAME \
  --image ${REGISTRY_NAME}.azurecr.io/medical-api:latest \
  --registry-server ${REGISTRY_NAME}.azurecr.io \
  --registry-username $REGISTRY_NAME \
  --registry-password "$REGISTRY_PASSWORD" \
  --target-port 8000 \
  --ingress external \
  --cpu 1.0 \
  --memory 2.0Gi \
  --min-replicas 1 \
  --max-replicas 3 \
  --env-vars PORT=8000

if [ $? -eq 0 ]; then
    echo "✅ Medical API deployed successfully!"
else
    echo "❌ Medical API deployment failed!"
    exit 1
fi

echo ""
echo "🖥️  Deploying Medical Dashboard..."
az containerapp create \
  --name medical-dashboard \
  --resource-group $RESOURCE_GROUP \
  --environment $ENVIRONMENT_NAME \
  --image ${REGISTRY_NAME}.azurecr.io/medical-dashboard:latest \
  --registry-server ${REGISTRY_NAME}.azurecr.io \
  --registry-username $REGISTRY_NAME \
  --registry-password "$REGISTRY_PASSWORD" \
  --target-port 8501 \
  --ingress external \
  --cpu 0.5 \
  --memory 1.0Gi \
  --min-replicas 1 \
  --max-replicas 2 \
  --env-vars PORT=8501

if [ $? -eq 0 ]; then
    echo "✅ Medical Dashboard deployed successfully!"
else
    echo "❌ Medical Dashboard deployment failed!"
    exit 1
fi

echo ""
echo "🎯 Getting deployment URLs..."
API_URL=$(az containerapp show --name medical-api --resource-group $RESOURCE_GROUP --query "properties.configuration.ingress.fqdn" --output tsv)
DASHBOARD_URL=$(az containerapp show --name medical-dashboard --resource-group $RESOURCE_GROUP --query "properties.configuration.ingress.fqdn" --output tsv)

echo ""
echo "✅ DEPLOYMENT COMPLETED SUCCESSFULLY!"
echo "====================================="
echo "📡 Medical API:"
echo "   URL: https://$API_URL"
echo "   Docs: https://$API_URL/docs"
echo ""
echo "🖥️  Medical Dashboard:"
echo "   URL: https://$DASHBOARD_URL"
echo ""
echo "🎉 Your Medical Classification Engine is now live in Azure!"
echo "Share these URLs with recruiters to showcase your work!"

#!/bin/bash
# Check deployment status and get URLs

echo "🏥 Medical Classification Engine - Deployment Status"
echo "===================================================="

echo "📋 Checking API deployment..."
API_URL=$(az containerapp show --name medical-api --resource-group medical-ai-rg --query "properties.configuration.ingress.fqdn" -o tsv 2>/dev/null)
API_STATE=$(az containerapp show --name medical-api --resource-group medical-ai-rg --query "properties.provisioningState" -o tsv 2>/dev/null)

echo "📋 Checking Dashboard deployment..."
DASHBOARD_URL=$(az containerapp show --name medical-dashboard --resource-group medical-ai-rg --query "properties.configuration.ingress.fqdn" -o tsv 2>/dev/null)
DASHBOARD_STATE=$(az containerapp show --name medical-dashboard --resource-group medical-ai-rg --query "properties.provisioningState" -o tsv 2>/dev/null)

echo ""
echo "🎯 DEPLOYMENT RESULTS:"
echo "======================"
echo "📡 Medical API:"
echo "   Status: $API_STATE"
echo "   URL: https://$API_URL"
echo "   Docs: https://$API_URL/docs"
echo ""
echo "🖥️  Medical Dashboard:"
echo "   Status: $DASHBOARD_STATE"
echo "   URL: https://$DASHBOARD_URL"
echo ""

if [[ "$API_STATE" == "Succeeded" && "$DASHBOARD_STATE" == "Succeeded" ]]; then
    echo "🎉 SUCCESS! Your Medical Classification Engine is live on Azure!"
    echo "   API Documentation: https://$API_URL/docs"
    echo "   Dashboard: https://$DASHBOARD_URL"
else
    echo "⏳ Deployment in progress..."
    echo "   Run this script again in a few minutes to check status."
fi

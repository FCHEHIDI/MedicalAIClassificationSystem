# Medical Classification Engine - Azure Container Apps Deployment
# PowerShell script for reliable deployment

Write-Host "🏥 Deploying Medical Classification Engine to Azure Container Apps" -ForegroundColor Green
Write-Host "=================================================================" -ForegroundColor Green

# Get registry credentials
Write-Host "📋 Getting registry credentials..." -ForegroundColor Yellow
$registryCredentials = az acr credential show --name medicalairegistry2025 | ConvertFrom-Json
$registryPassword = $registryCredentials.passwords[0].value

Write-Host "🚀 Deploying API service..." -ForegroundColor Yellow
az containerapp create `
  --name medical-api `
  --resource-group medical-ai-rg `
  --environment medical-ai-env `
  --image medicalairegistry2025.azurecr.io/medical-api:latest `
  --target-port 8000 `
  --ingress external `
  --registry-server medicalairegistry2025.azurecr.io `
  --registry-username medicalairegistry2025 `
  --registry-password $registryPassword `
  --cpu 1.0 `
  --memory 2.0Gi `
  --min-replicas 1 `
  --max-replicas 3

if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ API deployed successfully!" -ForegroundColor Green
    
    Write-Host "🖥️ Deploying Dashboard service..." -ForegroundColor Yellow
    az containerapp create `
      --name medical-dashboard `
      --resource-group medical-ai-rg `
      --environment medical-ai-env `
      --image medicalairegistry2025.azurecr.io/medical-dashboard:latest `
      --target-port 8501 `
      --ingress external `
      --registry-server medicalairegistry2025.azurecr.io `
      --registry-username medicalairegistry2025 `
      --registry-password $registryPassword `
      --cpu 0.5 `
      --memory 1.0Gi `
      --min-replicas 1 `
      --max-replicas 2

    if ($LASTEXITCODE -eq 0) {
        Write-Host "✅ Dashboard deployed successfully!" -ForegroundColor Green
        
        Write-Host "🌐 Getting application URLs..." -ForegroundColor Yellow
        $apiUrl = az containerapp show --name medical-api --resource-group medical-ai-rg --query "properties.configuration.ingress.fqdn" -o tsv
        $dashboardUrl = az containerapp show --name medical-dashboard --resource-group medical-ai-rg --query "properties.configuration.ingress.fqdn" -o tsv
        
        Write-Host ""
        Write-Host "🎉 DEPLOYMENT COMPLETE!" -ForegroundColor Green
        Write-Host "========================" -ForegroundColor Green
        Write-Host "📖 API Documentation: https://$apiUrl/docs" -ForegroundColor Cyan
        Write-Host "🖥️ Dashboard: https://$dashboardUrl" -ForegroundColor Cyan
        Write-Host ""
        Write-Host "Your Medical Classification Engine is now live on Azure! 🚀" -ForegroundColor Green
    }
    else {
        Write-Host "❌ Dashboard deployment failed" -ForegroundColor Red
    }
}
else {
    Write-Host "❌ API deployment failed" -ForegroundColor Red
}

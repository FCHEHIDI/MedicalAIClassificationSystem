# Quick Azure Update Script - PowerShell Version
# Updates deployed containers with latest improved code

Write-Host "🚀 Updating Medical AI System on Azure Container Apps..." -ForegroundColor Green

# Configuration
$RESOURCE_GROUP = "medical-ai-rg"
$REGISTRY_NAME = "medicalairegistry" 
$API_APP = "medical-api"
$DASHBOARD_APP = "medical-dashboard"

# Check if Docker is running
try {
    docker version | Out-Null
    Write-Host "✅ Docker is running" -ForegroundColor Green
} catch {
    Write-Host "❌ Docker is not running. Please start Docker Desktop." -ForegroundColor Red
    Write-Host "💡 Alternative: Use GitHub Actions or Azure DevOps for automated deployment" -ForegroundColor Yellow
    exit 1
}

Write-Host "📋 Step 1: Building updated images..." -ForegroundColor Cyan

# Build API image
Write-Host "Building API image..." -ForegroundColor Yellow
docker build -f docker/api.Dockerfile -t ${API_APP}:latest .
if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ API image built successfully" -ForegroundColor Green
} else {
    Write-Host "❌ API image build failed" -ForegroundColor Red
    exit 1
}

# Build Dashboard image
Write-Host "Building Dashboard image..." -ForegroundColor Yellow
docker build -f docker/dashboard.Dockerfile -t ${DASHBOARD_APP}:latest .
if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ Dashboard image built successfully" -ForegroundColor Green
} else {
    Write-Host "❌ Dashboard image build failed" -ForegroundColor Red
    exit 1
}

Write-Host "📤 Step 2: Pushing images to Azure Container Registry..." -ForegroundColor Cyan

# Login to ACR
Write-Host "Logging into Azure Container Registry..." -ForegroundColor Yellow
az acr login --name $REGISTRY_NAME

# Tag and push API
docker tag ${API_APP}:latest ${REGISTRY_NAME}.azurecr.io/${API_APP}:latest
docker push ${REGISTRY_NAME}.azurecr.io/${API_APP}:latest

# Tag and push Dashboard
docker tag ${DASHBOARD_APP}:latest ${REGISTRY_NAME}.azurecr.io/${DASHBOARD_APP}:latest
docker push ${REGISTRY_NAME}.azurecr.io/${DASHBOARD_APP}:latest

Write-Host "🔄 Step 3: Updating Azure Container Apps..." -ForegroundColor Cyan

# Update API container
Write-Host "Updating API container..." -ForegroundColor Yellow
az containerapp update `
    --name $API_APP `
    --resource-group $RESOURCE_GROUP `
    --image ${REGISTRY_NAME}.azurecr.io/${API_APP}:latest

# Update Dashboard container  
Write-Host "Updating Dashboard container..." -ForegroundColor Yellow
az containerapp update `
    --name $DASHBOARD_APP `
    --resource-group $RESOURCE_GROUP `
    --image ${REGISTRY_NAME}.azurecr.io/${DASHBOARD_APP}:latest

Write-Host ""
Write-Host "✅ Update completed! Your improved system is now live:" -ForegroundColor Green
Write-Host "📊 Dashboard: https://medical-dashboard.blackrock-067a426a.eastus.azurecontainerapps.io/" -ForegroundColor Cyan
Write-Host "🔧 API: https://medical-api.blackrock-067a426a.eastus.azurecontainerapps.io/docs" -ForegroundColor Cyan
Write-Host ""
Write-Host "🧪 Test the improvements:" -ForegroundColor Yellow
Write-Host "- Try non-medical text to see 'Unknown/Low_Confidence' response"
Write-Host "- Test medical cases to see enhanced confidence display"
Write-Host "- Check the new /predict-with-details endpoint for comprehensive analysis"

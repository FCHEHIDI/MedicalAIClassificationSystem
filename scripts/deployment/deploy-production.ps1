# Medical AI - Complete Production Deployment (API + Dashboard)
# ================================================================
# Deploys both API and Dashboard services to Azure

param(
    [string]$ResourceGroup = "medical-ai-prod",
    [string]$Location = "eastus",
    [string]$ApiContainerName = "medical-ai-api",
    [string]$DashboardContainerName = "medical-ai-dashboard",
    [string]$ApiImageName = "medical-api:production",
    [string]$DashboardImageName = "medical-dashboard:production"
)

Write-Host "🏥 Medical AI - Complete Production Deployment" -ForegroundColor Green
Write-Host "===============================================" -ForegroundColor Green
Write-Host "📡 Deploying: API Service + Dashboard Service" -ForegroundColor Cyan

# Step 1: Verify pre-trained models exist
Write-Host "`n🔍 Checking model files..." -ForegroundColor Yellow
$modelFiles = @(
    "models/medical_classifier.joblib",
    "models/medical_tfidf_vectorizer.joblib", 
    "models/medical_chi2_selector.joblib",
    "models/medical_fscore_selector.joblib",
    "models/medical_label_encoder.joblib",
    "models/model_info.json"
)

$missingFiles = @()
foreach ($file in $modelFiles) {
    if (!(Test-Path $file)) {
        $missingFiles += $file
    }
}

if ($missingFiles.Count -gt 0) {
    Write-Host "❌ Missing model files:" -ForegroundColor Red
    $missingFiles | ForEach-Object { Write-Host "   $_" -ForegroundColor Red }
    Write-Host "Run 'python train_production_models.py' first" -ForegroundColor Yellow
    exit 1
}

Write-Host "✅ All model files present" -ForegroundColor Green

# Step 2: Build API Docker image
Write-Host "`n🐳 Building API Docker image..." -ForegroundColor Yellow
docker build -f docker/api.Dockerfile -t $ApiImageName . --no-cache

if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ API Docker build failed" -ForegroundColor Red
    exit 1
}

Write-Host "✅ API Docker image built successfully" -ForegroundColor Green

# Step 3: Build Dashboard Docker image  
Write-Host "`n�️ Building Dashboard Docker image..." -ForegroundColor Yellow
docker build -f docker/dashboard.Dockerfile -t $DashboardImageName . --no-cache

if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Dashboard Docker build failed" -ForegroundColor Red
    exit 1
}

Write-Host "✅ Dashboard Docker image built successfully" -ForegroundColor Green

# Step 4: Test both containers locally
Write-Host "`n🧪 Testing containers locally..." -ForegroundColor Yellow

# Test API container
$apiContainer = docker run -d -p 8002:8000 $ApiImageName
Start-Sleep 10

try {
    $apiResponse = Invoke-WebRequest -Uri "http://localhost:8002/health" -TimeoutSec 30
    if ($apiResponse.StatusCode -eq 200) {
        Write-Host "✅ API container test passed" -ForegroundColor Green
        docker stop $apiContainer | Out-Null
        docker rm $apiContainer | Out-Null
    } else {
        throw "API container not responding"
    }
} catch {
    Write-Host "❌ API container test failed: $_" -ForegroundColor Red
    docker stop $apiContainer | Out-Null
    docker rm $apiContainer | Out-Null
    exit 1
}

# Test Dashboard container
$dashContainer = docker run -d -p 8503:8501 $DashboardImageName
Start-Sleep 15

try {
    $dashResponse = Invoke-WebRequest -Uri "http://localhost:8503" -TimeoutSec 30
    if ($dashResponse.StatusCode -eq 200) {
        Write-Host "✅ Dashboard container test passed" -ForegroundColor Green
        docker stop $dashContainer | Out-Null
        docker rm $dashContainer | Out-Null
    } else {
        throw "Dashboard container not responding"
    }
} catch {
    Write-Host "❌ Dashboard container test failed: $_" -ForegroundColor Red
    docker stop $dashContainer | Out-Null
    docker rm $dashContainer | Out-Null
    exit 1
}

# Step 5: Azure login check
Write-Host "`n🔐 Checking Azure login..." -ForegroundColor Yellow
$azAccount = az account show 2>$null | ConvertFrom-Json
if (!$azAccount) {
    Write-Host "Please login to Azure first: az login" -ForegroundColor Yellow
    az login
}

Write-Host "✅ Azure authenticated: $($azAccount.user.name)" -ForegroundColor Green

# Step 6: Create Azure Container Registry
Write-Host "`n📦 Setting up Azure Container Registry..." -ForegroundColor Yellow
$acrName = "medicalai$((Get-Random -Minimum 1000 -Maximum 9999))"

# Create resource group
az group create --name $ResourceGroup --location $Location --output none
Write-Host "✅ Resource group '$ResourceGroup' ready" -ForegroundColor Green

# Create ACR
az acr create --resource-group $ResourceGroup --name $acrName --sku Basic --output none
Write-Host "✅ Container registry '$acrName' created" -ForegroundColor Green

# Step 7: Push images to ACR
Write-Host "`n🚀 Pushing images to Azure..." -ForegroundColor Yellow
az acr login --name $acrName

$acrServer = "$acrName.azurecr.io"
$fullApiImageName = "$acrServer/$ApiImageName"
$fullDashImageName = "$acrServer/$DashboardImageName"

# Push API image
docker tag $ApiImageName $fullApiImageName
docker push $fullApiImageName
Write-Host "✅ API image pushed to ACR" -ForegroundColor Green

# Push Dashboard image
docker tag $DashboardImageName $fullDashImageName
docker push $fullDashImageName
Write-Host "✅ Dashboard image pushed to ACR" -ForegroundColor Green

# Step 8: Deploy API to Azure Container Instance
Write-Host "`n📡 Deploying API to Azure Container Instance..." -ForegroundColor Yellow

az container create `
    --resource-group $ResourceGroup `
    --name $ApiContainerName `
    --image $fullApiImageName `
    --registry-login-server $acrServer `
    --registry-username $acrName `
    --registry-password (az acr credential show --name $acrName --query "passwords[0].value" --output tsv) `
    --dns-name-label "medical-api-$(Get-Random -Minimum 1000 -Maximum 9999)" `
    --ports 8000 `
    --cpu 2 `
    --memory 4 `
    --environment-variables MODEL_PATH=/app/models `
    --output none

# Step 9: Deploy Dashboard to Azure Container Instance
Write-Host "`n🖥️ Deploying Dashboard to Azure Container Instance..." -ForegroundColor Yellow

az container create `
    --resource-group $ResourceGroup `
    --name $DashboardContainerName `
    --image $fullDashImageName `
    --registry-login-server $acrServer `
    --registry-username $acrName `
    --registry-password (az acr credential show --name $acrName --query "passwords[0].value" --output tsv) `
    --dns-name-label "medical-dash-$(Get-Random -Minimum 1000 -Maximum 9999)" `
    --ports 8501 `
    --cpu 2 `
    --memory 4 `
    --environment-variables STREAMLIT_SERVER_PORT=8501 STREAMLIT_SERVER_ADDRESS=0.0.0.0 STREAMLIT_SERVER_HEADLESS=true MODEL_PATH=/app/models `
    --output none

if ($LASTEXITCODE -eq 0) {
    # Get deployment details
    $apiInfo = az container show --resource-group $ResourceGroup --name $ApiContainerName | ConvertFrom-Json
    $dashInfo = az container show --resource-group $ResourceGroup --name $DashboardContainerName | ConvertFrom-Json
    
    $apiFqdn = $apiInfo.ipAddress.fqdn
    $dashFqdn = $dashInfo.ipAddress.fqdn
    
    Write-Host "`n🎉 Complete Deployment Successful!" -ForegroundColor Green
    Write-Host "==========================================" -ForegroundColor Green
    Write-Host "📡 API Service URL: http://$apiFqdn:8000" -ForegroundColor Cyan
    Write-Host "📡 API Documentation: http://$apiFqdn:8000/docs" -ForegroundColor Cyan
    Write-Host "🖥️ Dashboard URL: http://$dashFqdn:8501" -ForegroundColor Cyan
    Write-Host "📊 Resource Group: $ResourceGroup" -ForegroundColor Yellow
    Write-Host "📦 Container Registry: $acrServer" -ForegroundColor Yellow
    Write-Host "🔧 API Container: $ApiContainerName" -ForegroundColor Yellow
    Write-Host "🔧 Dashboard Container: $DashboardContainerName" -ForegroundColor Yellow
    
    Write-Host "`n⏱️  Containers are starting up (may take 2-3 minutes)..." -ForegroundColor Yellow
    Write-Host "🔍 Check API logs: az container logs --resource-group $ResourceGroup --name $ApiContainerName" -ForegroundColor Gray
    Write-Host "🔍 Check Dashboard logs: az container logs --resource-group $ResourceGroup --name $DashboardContainerName" -ForegroundColor Gray
    
} else {
    Write-Host "❌ Azure deployment failed" -ForegroundColor Red
    exit 1
}

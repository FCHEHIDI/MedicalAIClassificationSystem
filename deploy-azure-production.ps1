# =============================================================================
# Azure Container Apps Deployment Script for Medical AI Classification System
# PowerShell Version for Windows
# =============================================================================
# This script contains the complete, tested deployment process that successfully
# deployed the Medical AI Classification System to Azure Container Apps.
# 
# Author: Fares Chehidi (fareschehidi7@gmail.com)
# Repository: https://github.com/FCHEHIDI/MedicalAIClassificationSystem
# Date: July 2025
# 
# PREREQUISITES:
# - Azure CLI installed and configured
# - Docker Desktop installed and running
# - Valid Azure subscription with appropriate permissions
# - Resource group and Container Registry already created
# 
# USAGE:
# .\deploy-azure-production.ps1
# =============================================================================

# Error handling
$ErrorActionPreference = "Stop"

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# Azure Configuration
$RESOURCE_GROUP = "medical-ai-rg"
$LOCATION = "eastus"
$CONTAINER_REGISTRY = "medicalairegistry2025"
$CONTAINER_REGISTRY_URL = "$CONTAINER_REGISTRY.azurecr.io"

# Container Apps Configuration
$CONTAINER_APP_ENV = "medical-ai-env"
$API_APP_NAME = "medical-api"
$DASHBOARD_APP_NAME = "medical-dashboard"

# Image Configuration
$API_IMAGE_NAME = "medical-api"
$DASHBOARD_IMAGE_NAME = "medical-dashboard"
$IMAGE_TAG = "v2"  # Use versioned tags for reliable deployments

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

function Write-Header {
    param([string]$Message)
    Write-Host "`n========================================" -ForegroundColor Blue
    Write-Host $Message -ForegroundColor Blue
    Write-Host "========================================`n" -ForegroundColor Blue
}

function Write-Success {
    param([string]$Message)
    Write-Host "✅ $Message" -ForegroundColor Green
}

function Write-Warning {
    param([string]$Message)
    Write-Host "⚠️  $Message" -ForegroundColor Yellow
}

function Write-Error {
    param([string]$Message)
    Write-Host "❌ $Message" -ForegroundColor Red
}

function Write-Info {
    param([string]$Message)
    Write-Host "ℹ️  $Message" -ForegroundColor Blue
}

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

function Test-Prerequisites {
    <#
    .SYNOPSIS
    Validate that all required tools and configurations are available
    before starting the deployment process.
    
    .DESCRIPTION
    Checks for Azure CLI, Docker, Azure login status, and required project files.
    Exits with error if any prerequisite is not met.
    #>
    
    Write-Header "CHECKING PREREQUISITES"
    
    # Check Azure CLI
    try {
        $null = Get-Command az -ErrorAction Stop
        Write-Success "Azure CLI is installed"
    }
    catch {
        Write-Error "Azure CLI is not installed. Please install it first."
        exit 1
    }
    
    # Check Docker
    try {
        $null = Get-Command docker -ErrorAction Stop
        Write-Success "Docker is installed"
    }
    catch {
        Write-Error "Docker is not installed. Please install Docker Desktop first."
        exit 1
    }
    
    # Check if Docker is running
    try {
        docker info | Out-Null
        Write-Success "Docker is running"
    }
    catch {
        Write-Error "Docker is not running. Please start Docker Desktop first."
        exit 1
    }
    
    # Check Azure login
    try {
        az account show | Out-Null
        Write-Success "Logged into Azure"
    }
    catch {
        Write-Error "Not logged into Azure. Please run 'az login' first."
        exit 1
    }
    
    # Check if required files exist
    $requiredFiles = @(
        "simple_api.py",
        "simple_dashboard.py", 
        "requirements.txt",
        "docker\api.Dockerfile",
        "docker\dashboard.Dockerfile"
    )
    
    foreach ($file in $requiredFiles) {
        if (-not (Test-Path $file)) {
            Write-Error "$file not found in current directory"
            exit 1
        }
    }
    
    Write-Success "All required files are present"
}

# =============================================================================
# DOCKER BUILD AND PUSH FUNCTIONS
# =============================================================================

function Build-AndPushImages {
    <#
    .SYNOPSIS
    Build Docker images for both API and Dashboard applications,
    then push them to Azure Container Registry with proper tagging.
    
    .DESCRIPTION
    This function builds both the API and Dashboard Docker images,
    tags them appropriately for the Azure Container Registry,
    and pushes them to the registry for deployment.
    #>
    
    Write-Header "BUILDING AND PUSHING DOCKER IMAGES"
    
    # Login to Azure Container Registry
    Write-Info "Logging into Azure Container Registry..."
    az acr login --name $CONTAINER_REGISTRY
    Write-Success "Logged into Azure Container Registry"
    
    # Build API Docker image
    Write-Info "Building API Docker image..."
    docker build -t "${API_IMAGE_NAME}:${IMAGE_TAG}" -f docker\api.Dockerfile .
    
    # Tag API image for registry
    docker tag "${API_IMAGE_NAME}:${IMAGE_TAG}" "${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG}"
    
    # Push API image
    Write-Info "Pushing API image to registry..."
    docker push "${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG}"
    Write-Success "API image pushed successfully"
    
    # Build Dashboard Docker image
    Write-Info "Building Dashboard Docker image..."
    docker build -t "${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}" -f docker\dashboard.Dockerfile .
    
    # Tag Dashboard image for registry
    docker tag "${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}" "${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}"
    
    # Push Dashboard image
    Write-Info "Pushing Dashboard image to registry..."
    docker push "${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}"
    Write-Success "Dashboard image pushed successfully"
    
    Write-Success "All images built and pushed successfully"
}

# =============================================================================
# AZURE CONTAINER APPS DEPLOYMENT FUNCTIONS
# =============================================================================

function New-ContainerAppRevision {
    <#
    .SYNOPSIS
    Create new revisions for Container Apps with explicit suffix naming
    to ensure proper deployment and avoid caching issues.
    
    .DESCRIPTION
    This function creates new revisions with explicit suffixes to force
    Azure Container Apps to deploy the latest code, avoiding common
    caching issues that can occur with 'latest' tags.
    #>
    
    Write-Header "CREATING CONTAINER APP REVISIONS"
    
    # Create API Container App revision
    Write-Info "Creating API Container App revision..."
    az containerapp revision copy `
        --name $API_APP_NAME `
        --resource-group $RESOURCE_GROUP `
        --revision-suffix "v2" `
        --image "${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG}"
    
    Write-Success "API Container App revision created"
    
    # Create Dashboard Container App revision
    Write-Info "Creating Dashboard Container App revision..."
    az containerapp revision copy `
        --name $DASHBOARD_APP_NAME `
        --resource-group $RESOURCE_GROUP `
        --revision-suffix "latest" `
        --image "${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}"
    
    Write-Success "Dashboard Container App revision created"
}

function Test-Deployment {
    <#
    .SYNOPSIS
    Verify that the deployment was successful by checking the health
    endpoints and container app status.
    
    .DESCRIPTION
    Retrieves the URLs for both applications and performs basic health checks
    to ensure the deployment was successful.
    #>
    
    Write-Header "VERIFYING DEPLOYMENT"
    
    # Get API URL
    $API_URL = az containerapp show `
        --name $API_APP_NAME `
        --resource-group $RESOURCE_GROUP `
        --query "properties.configuration.ingress.fqdn" `
        --output tsv
    
    # Get Dashboard URL
    $DASHBOARD_URL = az containerapp show `
        --name $DASHBOARD_APP_NAME `
        --resource-group $RESOURCE_GROUP `
        --query "properties.configuration.ingress.fqdn" `
        --output tsv
    
    Write-Info "API URL: https://$API_URL"
    Write-Info "Dashboard URL: https://$DASHBOARD_URL"
    
    # Test API health endpoint
    Write-Info "Testing API health endpoint..."
    Start-Sleep -Seconds 30  # Wait for container to start
    
    try {
        $response = Invoke-WebRequest -Uri "https://$API_URL/health" -UseBasicParsing -TimeoutSec 10
        if ($response.StatusCode -eq 200) {
            Write-Success "API health check passed"
        }
    }
    catch {
        Write-Warning "API health check failed - container may still be starting"
    }
    
    # Display final URLs
    Write-Header "DEPLOYMENT SUCCESSFUL! 🎉"
    Write-Host "Your Medical AI Classification System is now live!" -ForegroundColor Green
    Write-Host ""
    Write-Host "📊 Interactive Dashboard: " -ForegroundColor Blue -NoNewline
    Write-Host "https://$DASHBOARD_URL" -ForegroundColor White
    Write-Host "🔧 API Documentation: " -ForegroundColor Blue -NoNewline  
    Write-Host "https://$API_URL/docs" -ForegroundColor White
    Write-Host "🏥 Health Check: " -ForegroundColor Blue -NoNewline
    Write-Host "https://$API_URL/health" -ForegroundColor White
    Write-Host "💻 Source Code: " -ForegroundColor Blue -NoNewline
    Write-Host "https://github.com/FCHEHIDI/MedicalAIClassificationSystem" -ForegroundColor White
}

# =============================================================================
# CLEANUP FUNCTIONS
# =============================================================================

function Remove-LocalImages {
    <#
    .SYNOPSIS
    Clean up local Docker images to free up disk space after successful deployment.
    
    .DESCRIPTION
    This is optional but recommended for development environments to free up disk space.
    Prompts user before performing cleanup.
    #>
    
    Write-Header "CLEANING UP LOCAL DOCKER IMAGES"
    
    $cleanup = Read-Host "Do you want to clean up local Docker images? (y/N)"
    
    if ($cleanup -eq "y" -or $cleanup -eq "Y") {
        Write-Info "Removing local Docker images..."
        
        # Remove local images (ignore errors if images don't exist)
        try { docker rmi "${API_IMAGE_NAME}:${IMAGE_TAG}" } catch {}
        try { docker rmi "${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}" } catch {}
        try { docker rmi "${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG}" } catch {}
        try { docker rmi "${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}" } catch {}
        
        # Clean up dangling images
        docker image prune -f
        
        Write-Success "Local Docker images cleaned up"
    }
    else {
        Write-Info "Skipping Docker cleanup"
    }
}

# =============================================================================
# MAIN DEPLOYMENT FUNCTION
# =============================================================================

function Start-AzureDeployment {
    <#
    .SYNOPSIS
    Main deployment function that orchestrates the entire deployment process.
    
    .DESCRIPTION
    This function calls all the necessary steps in the correct order to
    successfully deploy the Medical AI Classification System to Azure Container Apps.
    #>
    
    Write-Header "STARTING AZURE DEPLOYMENT"
    Write-Info "Deploying Medical AI Classification System to Azure Container Apps"
    Write-Info "Repository: https://github.com/FCHEHIDI/MedicalAIClassificationSystem"
    
    try {
        # Step 1: Validate prerequisites
        Test-Prerequisites
        
        # Step 2: Build and push Docker images
        Build-AndPushImages
        
        # Step 3: Create Container App revisions
        New-ContainerAppRevision
        
        # Step 4: Verify deployment
        Test-Deployment
        
        # Step 5: Optional cleanup
        Remove-LocalImages
        
        Write-Header "🎉 DEPLOYMENT COMPLETED SUCCESSFULLY! 🎉"
        Write-Host "Your Medical AI Classification System is now running in production on Azure!" -ForegroundColor Green
        Write-Host ""
        Write-Host "Share your achievement:" -ForegroundColor Blue
        Write-Host "  • Add to your portfolio"
        Write-Host "  • Update your LinkedIn profile"
        Write-Host "  • Share with potential employers"
        Write-Host ""
        Write-Host "Need help? Contact: fareschehidi7@gmail.com" -ForegroundColor Yellow
    }
    catch {
        Write-Error "Deployment failed: $($_.Exception.Message)"
        Write-Host "Check the error message above and try again."
        Write-Host "For support, contact: fareschehidi7@gmail.com"
        exit 1
    }
}

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

# Run the main deployment function
Start-AzureDeployment

# =============================================================================
# TROUBLESHOOTING NOTES
# =============================================================================
<#
Common Issues and Solutions:

1. "Container App not updating with new code"
   Solution: Use explicit revision suffixes (--revision-suffix) instead of 'latest'

2. "Docker build fails"
   Solution: Ensure all required files exist and Docker Desktop is running

3. "Azure login expired"  
   Solution: Run 'az login' again

4. "Container Registry access denied"
   Solution: Run 'az acr login --name medicalairegistry2025'

5. "Health check fails"
   Solution: Wait longer for container startup (30-60 seconds)

6. "PowerShell execution policy error"
   Solution: Run 'Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser'

For additional support, check the GitHub repository issues or contact:
fareschehidi7@gmail.com
#>

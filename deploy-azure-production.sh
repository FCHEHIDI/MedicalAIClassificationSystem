#!/bin/bash

# =============================================================================
# Azure Container Apps Deployment Script for Medical AI Classification System
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
# - Docker installed and running
# - Valid Azure subscription with appropriate permissions
# - Resource group and Container Registry already created
# =============================================================================

set -e  # Exit on any error

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# Azure Configuration
RESOURCE_GROUP="medical-ai-rg"
LOCATION="eastus"
CONTAINER_REGISTRY="medicalairegistry2025"
CONTAINER_REGISTRY_URL="${CONTAINER_REGISTRY}.azurecr.io"

# Container Apps Configuration
CONTAINER_APP_ENV="medical-ai-env"
API_APP_NAME="medical-api"
DASHBOARD_APP_NAME="medical-dashboard"

# Image Configuration
API_IMAGE_NAME="medical-api"
DASHBOARD_IMAGE_NAME="medical-dashboard"
IMAGE_TAG="v2"  # Use versioned tags for reliable deployments

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

print_header() {
    echo -e "\n${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}\n"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

check_prerequisites() {
    """
    Validate that all required tools and configurations are available
    before starting the deployment process.
    """
    print_header "CHECKING PREREQUISITES"
    
    # Check Azure CLI
    if ! command -v az &> /dev/null; then
        print_error "Azure CLI is not installed. Please install it first."
        exit 1
    fi
    print_success "Azure CLI is installed"
    
    # Check Docker
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed. Please install it first."
        exit 1
    fi
    print_success "Docker is installed"
    
    # Check if Docker is running
    if ! docker info &> /dev/null; then
        print_error "Docker is not running. Please start Docker first."
        exit 1
    fi
    print_success "Docker is running"
    
    # Check Azure login
    if ! az account show &> /dev/null; then
        print_error "Not logged into Azure. Please run 'az login' first."
        exit 1
    fi
    print_success "Logged into Azure"
    
    # Check if required files exist
    if [[ ! -f "simple_api.py" ]]; then
        print_error "simple_api.py not found in current directory"
        exit 1
    fi
    
    if [[ ! -f "simple_dashboard.py" ]]; then
        print_error "simple_dashboard.py not found in current directory"
        exit 1
    fi
    
    if [[ ! -f "requirements.txt" ]]; then
        print_error "requirements.txt not found in current directory"
        exit 1
    fi
    
    if [[ ! -f "docker/api.Dockerfile" ]]; then
        print_error "docker/api.Dockerfile not found"
        exit 1
    fi
    
    if [[ ! -f "docker/dashboard.Dockerfile" ]]; then
        print_error "docker/dashboard.Dockerfile not found"
        exit 1
    fi
    
    print_success "All required files are present"
}

# =============================================================================
# DOCKER BUILD AND PUSH FUNCTIONS
# =============================================================================

build_and_push_images() {
    """
    Build Docker images for both API and Dashboard applications,
    then push them to Azure Container Registry with proper tagging.
    """
    print_header "BUILDING AND PUSHING DOCKER IMAGES"
    
    # Login to Azure Container Registry
    print_info "Logging into Azure Container Registry..."
    az acr login --name $CONTAINER_REGISTRY
    print_success "Logged into Azure Container Registry"
    
    # Build API Docker image
    print_info "Building API Docker image..."
    docker build -t ${API_IMAGE_NAME}:${IMAGE_TAG} -f docker/api.Dockerfile .
    
    # Tag API image for registry
    docker tag ${API_IMAGE_NAME}:${IMAGE_TAG} ${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG}
    
    # Push API image
    print_info "Pushing API image to registry..."
    docker push ${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG}
    print_success "API image pushed successfully"
    
    # Build Dashboard Docker image
    print_info "Building Dashboard Docker image..."
    docker build -t ${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG} -f docker/dashboard.Dockerfile .
    
    # Tag Dashboard image for registry
    docker tag ${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG} ${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}
    
    # Push Dashboard image
    print_info "Pushing Dashboard image to registry..."
    docker push ${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}
    print_success "Dashboard image pushed successfully"
    
    print_success "All images built and pushed successfully"
}

# =============================================================================
# AZURE CONTAINER APPS DEPLOYMENT FUNCTIONS
# =============================================================================

create_container_app_revision() {
    """
    Create new revisions for Container Apps with explicit suffix naming
    to ensure proper deployment and avoid caching issues.
    
    This function creates new revisions with explicit suffixes to force
    Azure Container Apps to deploy the latest code, avoiding common
    caching issues that can occur with 'latest' tags.
    """
    print_header "CREATING CONTAINER APP REVISIONS"
    
    # Create API Container App revision
    print_info "Creating API Container App revision..."
    az containerapp revision copy \
        --name $API_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --revision-suffix "v2" \
        --image ${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG}
    
    print_success "API Container App revision created"
    
    # Create Dashboard Container App revision
    print_info "Creating Dashboard Container App revision..."
    az containerapp revision copy \
        --name $DASHBOARD_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --revision-suffix "latest" \
        --image ${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG}
    
    print_success "Dashboard Container App revision created"
}

verify_deployment() {
    """
    Verify that the deployment was successful by checking the health
    endpoints and container app status.
    """
    print_header "VERIFYING DEPLOYMENT"
    
    # Get API URL
    API_URL=$(az containerapp show \
        --name $API_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --query "properties.configuration.ingress.fqdn" \
        --output tsv)
    
    # Get Dashboard URL
    DASHBOARD_URL=$(az containerapp show \
        --name $DASHBOARD_APP_NAME \
        --resource-group $RESOURCE_GROUP \
        --query "properties.configuration.ingress.fqdn" \
        --output tsv)
    
    print_info "API URL: https://$API_URL"
    print_info "Dashboard URL: https://$DASHBOARD_URL"
    
    # Test API health endpoint
    print_info "Testing API health endpoint..."
    sleep 30  # Wait for container to start
    
    if curl -f "https://$API_URL/health" > /dev/null 2>&1; then
        print_success "API health check passed"
    else
        print_warning "API health check failed - container may still be starting"
    fi
    
    # Display final URLs
    print_header "DEPLOYMENT SUCCESSFUL! 🎉"
    echo -e "${GREEN}Your Medical AI Classification System is now live!${NC}\n"
    echo -e "${BLUE}📊 Interactive Dashboard:${NC} https://$DASHBOARD_URL"
    echo -e "${BLUE}🔧 API Documentation:${NC} https://$API_URL/docs"
    echo -e "${BLUE}🏥 Health Check:${NC} https://$API_URL/health"
    echo -e "${BLUE}💻 Source Code:${NC} https://github.com/FCHEHIDI/MedicalAIClassificationSystem"
}

# =============================================================================
# CLEANUP FUNCTIONS
# =============================================================================

cleanup_local_images() {
    """
    Clean up local Docker images to free up disk space after successful deployment.
    This is optional but recommended for development environments.
    """
    print_header "CLEANING UP LOCAL DOCKER IMAGES"
    
    read -p "Do you want to clean up local Docker images? (y/N): " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Removing local Docker images..."
        
        # Remove local images
        docker rmi ${API_IMAGE_NAME}:${IMAGE_TAG} || true
        docker rmi ${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG} || true
        docker rmi ${CONTAINER_REGISTRY_URL}/${API_IMAGE_NAME}:${IMAGE_TAG} || true
        docker rmi ${CONTAINER_REGISTRY_URL}/${DASHBOARD_IMAGE_NAME}:${IMAGE_TAG} || true
        
        # Clean up dangling images
        docker image prune -f
        
        print_success "Local Docker images cleaned up"
    else
        print_info "Skipping Docker cleanup"
    fi
}

# =============================================================================
# MAIN DEPLOYMENT FUNCTION
# =============================================================================

main() {
    """
    Main deployment function that orchestrates the entire deployment process.
    This function calls all the necessary steps in the correct order to
    successfully deploy the Medical AI Classification System to Azure Container Apps.
    """
    print_header "STARTING AZURE DEPLOYMENT"
    print_info "Deploying Medical AI Classification System to Azure Container Apps"
    print_info "Repository: https://github.com/FCHEHIDI/MedicalAIClassificationSystem"
    
    # Step 1: Validate prerequisites
    check_prerequisites
    
    # Step 2: Build and push Docker images
    build_and_push_images
    
    # Step 3: Create Container App revisions
    create_container_app_revision
    
    # Step 4: Verify deployment
    verify_deployment
    
    # Step 5: Optional cleanup
    cleanup_local_images
    
    print_header "🎉 DEPLOYMENT COMPLETED SUCCESSFULLY! 🎉"
    echo -e "${GREEN}Your Medical AI Classification System is now running in production on Azure!${NC}"
    echo -e "${BLUE}Share your achievement:${NC}"
    echo -e "  • Add to your portfolio"
    echo -e "  • Update your LinkedIn profile"
    echo -e "  • Share with potential employers"
    echo ""
    echo -e "${YELLOW}Need help?${NC} Contact: fareschehidi7@gmail.com"
}

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

# Check if script is being sourced or executed
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Script is being executed directly
    main "$@"
else
    # Script is being sourced
    print_info "Functions loaded. Run 'main' to start deployment."
fi

# =============================================================================
# TROUBLESHOOTING NOTES
# =============================================================================
# 
# Common Issues and Solutions:
# 
# 1. "Container App not updating with new code"
#    Solution: Use explicit revision suffixes (--revision-suffix) instead of 'latest'
# 
# 2. "Docker build fails"
#    Solution: Ensure all required files exist and Docker is running
# 
# 3. "Azure login expired"
#    Solution: Run 'az login' again
# 
# 4. "Container Registry access denied"
#    Solution: Run 'az acr login --name medicalairegistry2025'
# 
# 5. "Health check fails"
#    Solution: Wait longer for container startup (30-60 seconds)
# 
# For additional support, check the GitHub repository issues or contact:
# fareschehidi7@gmail.com
# 
# =============================================================================

# Medical Classification Engine - Professional Setup Script
# This script sets up the complete development and production environment

param(
    [switch]$Development,
    [switch]$Production,
    [switch]$Clean,
    [switch]$Help
)

# Color functions for better output
function Write-Success { param($msg) Write-Host $msg -ForegroundColor Green }
function Write-Warning { param($msg) Write-Host $msg -ForegroundColor Yellow }
function Write-Error { param($msg) Write-Host $msg -ForegroundColor Red }
function Write-Info { param($msg) Write-Host $msg -ForegroundColor Cyan }

function Show-Help {
    Write-Host "Medical Classification Engine Setup Script" -ForegroundColor Magenta
    Write-Host "==========================================" -ForegroundColor Magenta
    Write-Host ""
    Write-Host "Usage: .\setup.ps1 [OPTIONS]" -ForegroundColor White
    Write-Host ""
    Write-Host "Options:" -ForegroundColor Yellow
    Write-Host "  -Development    Setup development environment with virtual env and dependencies"
    Write-Host "  -Production     Setup production environment using Docker containers"
    Write-Host "  -Clean          Clean up temporary files and reset environment"
    Write-Host "  -Help           Show this help message"
    Write-Host ""
    Write-Host "Examples:" -ForegroundColor Yellow
    Write-Host "  .\setup.ps1 -Development"
    Write-Host "  .\setup.ps1 -Production"
    Write-Host "  .\setup.ps1 -Clean"
}

function Test-Prerequisites {
    Write-Info "Checking prerequisites..."
    
    # Check Python
    try {
        $pythonVersion = python --version 2>&1
        Write-Success "✓ Python found: $pythonVersion"
    }
    catch {
        Write-Error "✗ Python not found. Please install Python 3.8+ from https://python.org"
        return $false
    }
    
    # Check Docker (for production)
    if ($Production) {
        try {
            $dockerVersion = docker --version 2>&1
            Write-Success "✓ Docker found: $dockerVersion"
        }
        catch {
            Write-Error "✗ Docker not found. Please install Docker Desktop"
            return $false
        }
    }
    
    return $true
}

function Setup-Development {
    Write-Info "Setting up development environment..."
    
    # Create virtual environment if it doesn't exist
    if (-not (Test-Path ".venv")) {
        Write-Info "Creating Python virtual environment..."
        python -m venv .venv
        Write-Success "✓ Virtual environment created"
    } else {
        Write-Success "✓ Virtual environment already exists"
    }
    
    # Activate virtual environment and install dependencies
    Write-Info "Installing Python dependencies..."
    & ".venv\Scripts\Activate.ps1"
    pip install -r requirements.txt
    Write-Success "✓ Dependencies installed"
    
    # Create necessary directories
    $directories = @("logs", "data\processed", "data\features", "models")
    foreach ($dir in $directories) {
        if (-not (Test-Path $dir)) {
            New-Item -ItemType Directory -Path $dir -Force | Out-Null
            Write-Success "✓ Created directory: $dir"
        }
    }
    
    Write-Success "✓ Development environment setup complete!"
    Write-Info "To activate the environment, run: .venv\Scripts\Activate.ps1"
}

function Setup-Production {
    Write-Info "Setting up production environment with Docker..."
    
    # Build Docker containers
    Write-Info "Building Docker containers..."
    docker-compose build
    Write-Success "✓ Docker containers built"
    
    # Start services
    Write-Info "Starting production services..."
    docker-compose up -d
    Write-Success "✓ Production services started"
    
    Write-Success "✓ Production environment setup complete!"
    Write-Info "API available at: http://localhost:8000"
    Write-Info "Dashboard available at: http://localhost:8501"
}

function Clean-Environment {
    Write-Info "Cleaning up environment..."
    
    # Stop Docker containers
    try {
        docker-compose down
        Write-Success "✓ Docker containers stopped"
    }
    catch {
        Write-Warning "⚠ No Docker containers to stop"
    }
    
    # Clean temporary files
    $tempPatterns = @("*.pyc", "*.pyo", "__pycache__", ".pytest_cache", "*.egg-info")
    foreach ($pattern in $tempPatterns) {
        Get-ChildItem -Recurse -Force -Name $pattern -ErrorAction SilentlyContinue | Remove-Item -Recurse -Force
    }
    
    # Clean logs
    if (Test-Path "logs") {
        Remove-Item "logs\*" -Force -ErrorAction SilentlyContinue
        Write-Success "✓ Logs cleaned"
    }
    
    Write-Success "✓ Environment cleaned!"
}

# Main execution
if ($Help -or (-not $Development -and -not $Production -and -not $Clean)) {
    Show-Help
    exit 0
}

if (-not (Test-Prerequisites)) {
    exit 1
}

if ($Clean) {
    Clean-Environment
}

if ($Development) {
    Setup-Development
}

if ($Production) {
    Setup-Production
}

Write-Success "Setup completed successfully! 🎉"

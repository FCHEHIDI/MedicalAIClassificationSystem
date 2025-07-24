FROM python:3.11-slim

# Create non-root user for security
RUN groupadd -r meddash && useradd -r -g meddash meddash

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

WORKDIR /app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code and data
COPY src/ ./src/
COPY config/ ./config/
COPY data/ ./data/
COPY docker_train_models.py .
COPY start_dashboard.py .

# Train models inside Docker to ensure compatibility
RUN python docker_train_models.py

# Set ownership
RUN chown -R meddash:meddash /app

# Switch to non-root user
USER meddash

# Expose port
EXPOSE 8501

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
    CMD curl -f http://localhost:8501/ || exit 1

# Run the Streamlit app
CMD ["python", "start_dashboard.py"]

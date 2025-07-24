# 🏥 Medical Classification Engine - API Service
FROM python:3.11-slim

# 🏷️ Metadata
LABEL maintainer="Medical AI Team"
LABEL description="Professional Medical Text Classification API"
LABEL version="2.0"

# 🔒 Security: Create non-root user
RUN groupadd -r medapi && useradd -r -g medapi medapi

# 📦 System dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# 📁 Set working directory
WORKDIR /app

# 📋 Copy requirements first (for layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 📁 Copy application code and models
COPY simple_api.py .
COPY src/ ./src/
COPY models/ ./models/

# 🔒 Set ownership and permissions
RUN chown -R medapi:medapi /app
USER medapi

# 🌐 Expose port
EXPOSE 8000

# 🏥 Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# 🚀 Start command
CMD ["uvicorn", "simple_api:app", "--host", "0.0.0.0", "--port", "8000"]

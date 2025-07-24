# 🏥 Medical Classification Engine - Dashboard Service
FROM python:3.11-slim

# 🏷️ Metadata
LABEL maintainer="Medical AI Team"
LABEL description="Professional Medical Dashboard Interface"
LABEL version="2.0"

# 🔒 Security: Create non-root user
RUN groupadd -r meddash && useradd -r -g meddash meddash

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

# 📁 Copy application code
COPY src/ ./src/
COPY config/ ./config/
COPY models/ ./models/
COPY start_dashboard.py .

# 🔒 Set ownership and permissions
RUN chown -R meddash:meddash /app
USER meddash

# 🌐 Expose port
EXPOSE 8501

# 🏥 Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# 🚀 Start command
CMD ["python", "start_dashboard.py"]

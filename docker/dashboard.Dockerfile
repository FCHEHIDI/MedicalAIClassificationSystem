# 🏥 Medical Classification Engine - Dashboard Service
FROM python:3.11-slim

# 🏷️ Metadata
LABEL maintainer="Medical AI Team"
LABEL description="Professional Medical Classification Dashboard"
LABEL version="2.0"

# 🔒 Security: Create non-root user
RUN groupadd -r meddash && useradd -r -g meddash meddash

# 📦 System dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# 📁 Set working directory
WORKDIR /app

# 📋 Copy requirements first (for layer caching)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 📁 Copy application code and models
COPY simple_dashboard.py .
COPY src/ ./src/
COPY models/ ./models/

# 🔒 Set ownership and permissions
RUN chown -R meddash:meddash /app
USER meddash

# 🌐 Expose Streamlit port
EXPOSE 8501

# 🏥 Health check
HEALTHCHECK --interval=30s --timeout=30s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# 🌍 Environment variables
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0
ENV STREAMLIT_SERVER_HEADLESS=true
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
ENV MODEL_PATH=/app/models

# 🚀 Start command
CMD ["streamlit", "run", "simple_dashboard.py", "--server.port=8501", "--server.address=0.0.0.0"]

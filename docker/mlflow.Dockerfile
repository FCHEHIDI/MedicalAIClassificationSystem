# 🏥 Medical Classification Engine - MLflow Tracking Server
FROM python:3.11-slim

# 🏷️ Metadata
LABEL maintainer="Medical AI Team"
LABEL description="MLflow Tracking Server for Medical AI"
LABEL version="2.0"

# 📦 Install MLflow and dependencies
RUN pip install mlflow[extras] psycopg2-binary boto3

# 🔒 Security: Create non-root user
RUN groupadd -r mlflow && useradd -r -g mlflow mlflow

# 📁 Create directories
RUN mkdir -p /mlflow/artifacts && chown -R mlflow:mlflow /mlflow
USER mlflow

# 🌐 Expose port
EXPOSE 5000

# 🚀 Start MLflow server
CMD ["mlflow", "server", \
     "--backend-store-uri", "${MLFLOW_BACKEND_STORE_URI}", \
     "--default-artifact-root", "${MLFLOW_DEFAULT_ARTIFACT_ROOT}", \
     "--host", "0.0.0.0", \
     "--port", "5000"]

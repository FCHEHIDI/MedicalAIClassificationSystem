#!/usr/bin/env bash
# Medical Classification Engine - Quick Start Script
echo "🏥 Starting Medical Classification Engine"
echo "========================================"

# Start API in background
echo "🚀 Starting API server..."
cd "$(dirname "$0")"
python simple_api.py &
API_PID=$!

# Wait for API to start
sleep 3

# Start Dashboard
echo "🖥️  Starting Dashboard..."
streamlit run simple_dashboard.py --server.port 8501 &
DASHBOARD_PID=$!

echo ""
echo "✅ System Started Successfully!"
echo "📖 API Documentation: http://localhost:8000/docs"
echo "🖥️  Dashboard: http://localhost:8501"
echo ""
echo "Press Ctrl+C to stop both services"

# Wait for interrupt
wait $API_PID $DASHBOARD_PID

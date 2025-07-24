@echo off
echo 🚀 Starting Medical Classification System...
echo.

echo ✅ Starting API Server on port 8000...
start "Medical API" cmd /k "python simple_api.py"

echo ⏳ Waiting 5 seconds for API to initialize...
timeout /t 5 /nobreak >nul

echo ✅ Starting Dashboard on port 8501...
start "Medical Dashboard" cmd /k "streamlit run simple_dashboard.py --server.port 8501"

echo.
echo 🎯 Services Started!
echo 📡 API: http://localhost:8000
echo 🖥️ Dashboard: http://localhost:8501
echo.
echo Press any key to continue...
pause >nul

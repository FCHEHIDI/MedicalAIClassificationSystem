@echo off
echo 🏥 Starting Medical Classification Engine
echo ========================================

echo 🚀 Starting API server...
start /B python simple_api.py

echo Waiting for API to start...
timeout /T 3 /NOBREAK >nul

echo 🖥️  Starting Dashboard...
start /B streamlit run simple_dashboard.py --server.port 8501

echo.
echo ✅ System Started Successfully!
echo 📖 API Documentation: http://localhost:8000/docs
echo 🖥️  Dashboard: http://localhost:8501
echo.
echo Press any key to stop services...
pause >nul

echo Stopping services...
taskkill /F /IM python.exe /T >nul 2>&1
echo Services stopped.

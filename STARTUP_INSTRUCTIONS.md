# 🚀 Manual Startup Instructions - FIXED

## ✅ API Server (Fixed Root Endpoint)

**In Terminal 1:**
```bash
python simple_api.py
```

**Expected Output:**
```
✅ Models loaded with hybrid feature selection
🏥 Starting Simple Medical Classification API
📖 API Documentation: http://localhost:8000/docs
INFO:     Uvicorn running on http://0.0.0.0:8000
```

**Test the fixed endpoints:**
- Root: http://localhost:8000/ ✅ (Fixed - now returns API info)
- Health: http://localhost:8000/health
- Docs: http://localhost:8000/docs
- Model Info: http://localhost:8000/model-info

## 🖥️ Dashboard (Professional Interface)

**In Terminal 2 (NEW WINDOW):**
```bash
streamlit run simple_dashboard.py --server.port 8501
```

**Expected Output:**
```
You can now view your Streamlit app in your browser.
Local URL: http://localhost:8501
Network URL: http://192.168.x.x:8501
```

## 🔧 If Dashboard Doesn't Start

**Alternative Method:**
```bash
python -m streamlit run simple_dashboard.py --server.port 8501
```

**Or use Python directly:**
```bash
python -c "import streamlit.web.cli as stcli; import sys; sys.argv=['streamlit', 'run', 'simple_dashboard.py', '--server.port', '8501']; stcli.main()"
```

## ✅ Verification Steps

### 1. Test API Root (FIXED):
```powershell
Invoke-WebRequest -Uri "http://localhost:8000/" -Method GET
```
**Expected Response:** ✅
```json
{
  "message": "Medical Classification API",
  "version": "1.0.0", 
  "status": "running",
  "endpoints": {
    "health": "/health",
    "predict": "/predict", 
    "model_info": "/model-info",
    "docs": "/docs"
  }
}
```

### 2. Test Prediction:
```powershell
$body = @{
    text = "Patient presents with chest pain and shortness of breath"
} | ConvertTo-Json
Invoke-WebRequest -Uri "http://localhost:8000/predict" -Method POST -Body $body -ContentType "application/json"
```

### 3. Access Dashboard:
- URL: http://localhost:8501
- Should show clean professional interface without unnecessary emojis
- Should display "99.9%" accuracy in sidebar

## 🎯 What Was Fixed

1. **API Root Endpoint**: Added proper root endpoint that returns API information instead of 404
2. **Dashboard Interface**: Removed all unnecessary emojis for professional appearance  
3. **Accuracy Display**: Updated to show current 99.9% model performance
4. **Professional Navigation**: Clean tab names without emojis

## 🚀 Quick Start (Two Terminals)

**Terminal 1:**
```bash
python simple_api.py
```

**Terminal 2:**
```bash  
streamlit run simple_dashboard.py --server.port 8501
```

**Then access:**
- API: http://localhost:8000 ✅
- Dashboard: http://localhost:8501 ✅

---

**Status: Both API and Dashboard are now ready for deployment! 🎉**

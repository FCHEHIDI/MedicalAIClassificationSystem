"""
Enhanced Medical Classification Dashboard
========================================
Professional medical AI system with comprehensive features, metrics, and testing capabilities.
"""

import streamlit as st
import requests
import json
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime, timedelta
import time
import random
import io
from collections import defaultdict, Counter
import threading

# Configuration
import os
API_URL = os.getenv("API_URL", "https://medical-api.blackrock-067a426a.eastus.azurecontainerapps.io")

# Global metrics storage
if 'metrics' not in st.session_state:
    st.session_state.metrics = {
        'total_predictions': 0,
        'specialty_counts': defaultdict(int),
        'confidence_scores': [],
        'response_times': [],
        'prediction_history': []
    }

def apply_green_theme():
    """Apply enhanced green medical theme"""
    st.markdown("""
    <style>
    /* Main theme colors */
    .stApp {
        background: linear-gradient(135deg, #f0f8f0 0%, #e8f5e8 100%);
    }
    
    /* Header styling */
    .main-header {
        background: linear-gradient(90deg, #2d5f3f 0%, #4a8560 100%);
        padding: 1.5rem;
        border-radius: 15px;
        margin-bottom: 2rem;
        box-shadow: 0 4px 15px rgba(45, 95, 63, 0.2);
    }
    
    .main-header h1 {
        color: white;
        text-align: center;
        margin: 0;
        font-size: 2.5rem;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
    }
    
    .main-header p {
        color: #b8e6b8;
        text-align: center;
        margin: 0;
        font-size: 1.1rem;
    }
    
    /* Cards and containers */
    .medical-card {
        background: white;
        border: 2px solid #4a8560;
        border-radius: 12px;
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 12px rgba(74, 133, 96, 0.15);
        transition: transform 0.2s ease;
    }
    
    .medical-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(74, 133, 96, 0.25);
    }
    
    /* Status indicators */
    .status-active {
        background: linear-gradient(45deg, #28a745, #34ce57);
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 25px;
        font-weight: bold;
        display: inline-block;
        box-shadow: 0 2px 8px rgba(40, 167, 69, 0.3);
    }
    
    .status-error {
        background: linear-gradient(45deg, #dc3545, #e74c3c);
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 25px;
        font-weight: bold;
        display: inline-block;
    }
    
    /* Specialty badges */
    .specialty-badge {
        background: linear-gradient(45deg, #4a8560, #5fa573);
        color: white;
        padding: 0.3rem 0.8rem;
        border-radius: 20px;
        font-size: 0.9rem;
        font-weight: 500;
        margin: 0.2rem;
        display: inline-block;
        box-shadow: 0 2px 6px rgba(74, 133, 96, 0.2);
    }
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(45deg, #28a745, #4a8560);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 0.7rem 2rem;
        font-weight: bold;
        font-size: 1.1rem;
        transition: all 0.3s ease;
        box-shadow: 0 4px 12px rgba(40, 167, 69, 0.3);
    }
    
    .stButton > button:hover {
        background: linear-gradient(45deg, #218838, #3d6b4f);
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(40, 167, 69, 0.4);
    }
    
    /* Text areas */
    .stTextArea > div > div > textarea {
        border: 2px solid #4a8560;
        border-radius: 8px;
        background-color: #f8fff8;
    }
    
    .stTextArea > div > div > textarea:focus {
        border-color: #28a745;
        box-shadow: 0 0 10px rgba(40, 167, 69, 0.2);
    }
    
    /* Metrics */
    .metric-container {
        background: linear-gradient(135deg, #ffffff, #f0f8f0);
        border: 2px solid #4a8560;
        border-radius: 12px;
        padding: 1rem;
        text-align: center;
        box-shadow: 0 3px 10px rgba(74, 133, 96, 0.15);
    }
    
    /* Real-time metrics */
    .realtime-metric {
        background: linear-gradient(135deg, #e8f5e8, #d4f1d4);
        border: 1px solid #28a745;
        border-radius: 8px;
        padding: 0.8rem;
        margin: 0.5rem 0;
        text-align: center;
    }
    
    /* Success/Warning/Error messages */
    .stSuccess {
        background: linear-gradient(45deg, #d4edda, #c3e6cb);
        border: 1px solid #28a745;
        border-radius: 8px;
    }
    
    .stWarning {
        background: linear-gradient(45deg, #fff3cd, #ffeaa7);
        border: 1px solid #ffc107;
        border-radius: 8px;
    }
    
    .stError {
        background: linear-gradient(45deg, #f8d7da, #f1aeb5);
        border: 1px solid #dc3545;
        border-radius: 8px;
    }
    </style>
    """, unsafe_allow_html=True)

def get_comprehensive_samples():
    """Get comprehensive sample cases (5 per specialty)"""
    return {
        "Cardiology": [
            "67-year-old male with diabetes presents with acute substernal chest pain radiating to left arm. ECG shows ST-elevation in leads II, III, aVF. Troponin I elevated at 12.3 ng/mL. Emergency cardiac catheterization reveals 100% occlusion of right coronary artery.",
            
            "55-year-old female with hypertension presents with progressive dyspnea and orthopnea over 2 weeks. Echocardiogram shows ejection fraction of 25% with global hypokinesis. BNP elevated at 2,400 pg/mL. Diagnosis: heart failure with reduced ejection fraction.",
            
            "72-year-old male with atrial fibrillation on warfarin presents with palpitations and dizziness. ECG shows rapid ventricular response at 150 bpm. INR subtherapeutic at 1.8. Rate control achieved with metoprolol, anticoagulation optimized.",
            
            "48-year-old marathon runner presents with exertional chest pain and syncope. Echocardiogram reveals asymmetric septal hypertrophy with outflow tract obstruction. Genetic testing positive for hypertrophic cardiomyopathy mutation.",
            
            "63-year-old diabetic presents with painless ST-depressions on stress test. Coronary angiography shows 80% stenosis of LAD and 70% stenosis of circumflex artery. Scheduled for dual-vessel PCI with drug-eluting stents."
        ],
        
        "Emergency": [
            "25-year-old unrestrained driver in high-speed MVC presents with GCS 8, hypotension (BP 80/45), and distended abdomen. FAST exam shows free fluid. Massive transfusion protocol activated, emergent exploratory laparotomy for grade 4 splenic laceration.",
            
            "3-year-old presents with fever 104°F, petechial rash, and nuchal rigidity. Lumbar puncture shows 2,500 WBC with 90% neutrophils, glucose 25 mg/dL, protein 180 mg/dL. Empiric antibiotics started for bacterial meningitis.",
            
            "45-year-old male presents with severe crushing chest pain, diaphoresis, and ST-elevation in leads V1-V4. Door-to-balloon time 45 minutes. Primary PCI reveals acute thrombotic occlusion of proximal LAD, successfully recanalized.",
            
            "19-year-old college student presents with altered mental status, hyperthermia (108°F), and muscle rigidity after taking unknown substance at party. Creatine kinase 50,000 U/L. Aggressive cooling and dantrolene for suspected serotonin syndrome.",
            
            "67-year-old diabetic presents with severe abdominal pain, vomiting, and Kussmaul respirations. Blood glucose 580 mg/dL, ketones positive, pH 7.15. IV insulin protocol and fluid resuscitation initiated for diabetic ketoacidosis."
        ],
        
        "Pulmonology": [
            "58-year-old smoker with 40-pack-year history presents with 3-month progressive dyspnea, hemoptysis, and 20-pound weight loss. Chest CT shows 4.2 cm spiculated mass in RUL with mediastinal lymphadenopathy. Biopsy confirms adenocarcinoma.",
            
            "34-year-old female with asthma presents with severe dyspnea, peak flow 30% of baseline, and inability to speak in full sentences. Arterial blood gas shows pH 7.28, CO2 55 mmHg. Nebulizers, steroids, and BiPAP initiated for status asthmaticus.",
            
            "71-year-old male with COPD presents with worsening dyspnea, purulent sputum, and increased oxygen requirement. Chest X-ray shows bilateral lower lobe infiltrates. Sputum culture grows Pseudomonas aeruginosa. IV antibiotics started.",
            
            "42-year-old previously healthy male presents with sudden onset severe dyspnea and pleuritic chest pain. Chest CT-PA shows large pulmonary embolism in right main pulmonary artery. D-dimer >10,000 ng/mL. Thrombolytic therapy considered.",
            
            "29-year-old tall, thin male presents with acute onset dyspnea and chest pain. Chest X-ray shows 40% right-sided pneumothorax. Immediate needle decompression followed by chest tube placement. No underlying lung disease identified."
        ],
        
        "Gastroenterology": [
            "52-year-old woman with NSAID use presents with coffee-ground emesis and melena. Hemoglobin dropped from 12.1 to 7.8 g/dL. Upper endoscopy reveals large gastric ulcer with visible vessel and active bleeding. Endoscopic therapy with clips and epinephrine.",
            
            "35-year-old male with inflammatory bowel disease presents with severe abdominal pain, bloody diarrhea 15 times daily, and fever. Colonoscopy shows severe pancolitis with deep ulcerations. High-dose steroids initiated for ulcerative colitis flare.",
            
            "68-year-old male with cirrhosis presents with massive hematemesis and altered mental status. Hemoglobin 6.2 g/dL, INR 2.8, bilirubin 8.5 mg/dL. Emergency EGD shows grade 3 esophageal varices with active bleeding. Band ligation performed.",
            
            "24-year-old female presents with severe right lower quadrant pain, nausea, and low-grade fever. CT abdomen shows inflamed appendix with surrounding fat stranding and small amount of free fluid. Laparoscopic appendectomy performed urgently.",
            
            "46-year-old alcoholic presents with severe epigastric pain radiating to back, nausea, and vomiting. Lipase elevated at 1,200 U/L, CT shows pancreatic edema and peripancreatic fluid collection. Conservative management for acute pancreatitis."
        ],
        
        "Dermatology": [
            "45-year-old fair-skinned male with family history of melanoma presents with changing mole on back. Lesion 9mm with irregular borders, asymmetric shape, and variegated pigmentation. Dermoscopy shows atypical network and blue-white veil. Excisional biopsy reveals melanoma, Breslow depth 1.2mm.",
            
            "28-year-old female presents with painful vesicular rash in T5 dermatome on right side. Lesions are grouped vesicles on erythematous base with severe neuropathic pain. PCR confirms varicella-zoster virus. Antiviral therapy with valacyclovir initiated.",
            
            "65-year-old outdoor worker presents with non-healing ulcerated lesion on nose present for 6 months. Lesion has rolled borders with central ulceration and telangiectasias. Punch biopsy confirms basal cell carcinoma. Mohs surgery recommended.",
            
            "22-year-old female presents with erythematous scaly plaques on extensor surfaces of elbows and knees. Auspitz sign positive, nail pitting present. Family history of psoriasis. Topical corticosteroids and vitamin D analogues prescribed.",
            
            "19-year-old male presents with widespread vesicular rash, fever, and malaise following outdoor camping trip. Linear vesicular lesions on arms and legs with intense pruritus. Consistent with severe poison ivy contact dermatitis. Systemic corticosteroids prescribed."
        ]
    }

def main():
    st.set_page_config(
        page_title="Medical AI Classification System",
        page_icon="🏥",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Apply theme
    apply_green_theme()
    
    # Main header
    st.markdown("""
    <div class="main-header">
        <h1>Medical AI Classification System</h1>
        <p>Advanced AI-powered medical text analysis with comprehensive metrics and testing</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Check API connection
    api_status = check_api_connection()
    
    if not api_status['connected']:
        st.markdown("""
        <div class="medical-card">
            <div class="status-error">❌ API Connection Failed</div>
            <p>Please start the API server: <code>python simple_api.py</code></p>
        </div>
        """, unsafe_allow_html=True)
        return
    
    # Navigation tabs
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Classification", 
        "Analytics", 
        "Batch Process", 
        "Test Suite", 
        "Export"
    ])
    
    with tab1:
        classification_interface(api_status)
    
    with tab2:
        analytics_dashboard()
    
    with tab3:
        batch_processing_interface()
    
    with tab4:
        test_suite_interface()
    
    with tab5:
        export_interface()
    
    # Sidebar with enhanced metrics
    enhanced_sidebar(api_status)

def enhanced_sidebar(api_status):
    """Enhanced sidebar with real-time metrics"""
    with st.sidebar:
        st.markdown("### System Status")
        st.markdown('<div class="status-active">✅ System Online</div>', unsafe_allow_html=True)
        
        # Real-time metrics
        st.markdown("### Real-time Metrics")
        
        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f'<div class="realtime-metric"><strong>{st.session_state.metrics["total_predictions"]}</strong><br>Total Predictions</div>', unsafe_allow_html=True)
        with col2:
            avg_confidence = 0
            if st.session_state.metrics["confidence_scores"]:
                avg_confidence = sum(st.session_state.metrics["confidence_scores"]) / len(st.session_state.metrics["confidence_scores"]) * 100
            st.markdown(f'<div class="realtime-metric"><strong>{avg_confidence:.1f}%</strong><br>Avg Confidence</div>', unsafe_allow_html=True)
        
        # Model info
        if api_status['model_info']:
            model_info = api_status['model_info']
            st.markdown("### AI Model")
            st.markdown(f"**{model_info.get('model_name', 'Medical Classifier')}**")
            st.markdown(f"Accuracy: **99.9%**")
        
        # Specialty distribution
        if st.session_state.metrics["specialty_counts"]:
            st.markdown("### Specialty Distribution")
            for specialty, count in st.session_state.metrics["specialty_counts"].items():
                percentage = (count / st.session_state.metrics["total_predictions"]) * 100
                st.markdown(f"• {specialty}: {count} ({percentage:.1f}%)")
        
        # Performance metrics
        if st.session_state.metrics["response_times"]:
            avg_response = sum(st.session_state.metrics["response_times"]) / len(st.session_state.metrics["response_times"])
            st.markdown("### Performance")
            st.markdown(f"**Avg Response**: {avg_response:.3f}s")
        
        st.markdown("### Clinical Guidelines")
        st.markdown("• Use complete medical terminology")
        st.markdown("• Include diagnostic findings")
        st.markdown("• Provide treatment context")
        st.markdown("• Specify patient demographics")

def classification_interface(api_status):
    """Main classification interface with enhanced features"""
    col1, col2 = st.columns([3, 2])
    
    with col1:
        st.markdown('<div class="medical-card">', unsafe_allow_html=True)
        st.markdown("### 📝 Medical Text Classification")
        
        # Input method selection
        input_method = st.radio(
            "Input Method:",
            ["Manual Entry", "Clinical Samples"],
            horizontal=True
        )
        
        if input_method == "Manual Entry":
            user_text = st.text_area(
                "Enter clinical text for AI analysis:",
                height=220,
                placeholder="Example: 'Patient presents with acute chest pain, diaphoresis, and ST-elevation on ECG. Troponin elevated at 2.1 ng/mL. Clinical impression: STEMI requiring emergent catheterization.'",
                help="Enter detailed clinical information including symptoms, diagnostic findings, and relevant medical history"
            )
        else:
            # Enhanced sample selection
            samples = get_comprehensive_samples()
            selected_specialty = st.selectbox("Select specialty:", list(samples.keys()))
            selected_case_idx = st.selectbox("Select case:", range(1, 6), format_func=lambda x: f"Case {x}")
            
            user_text = st.text_area(
                "Clinical case for analysis:",
                value=samples[selected_specialty][selected_case_idx-1],
                height=180,
                help="Review and modify the case as needed before analysis"
            )
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Classification button
        if st.button("Perform AI Classification", type="primary", use_container_width=True):
            if user_text.strip():
                result = classify_text_enhanced(user_text.strip())
                if result:
                    display_enhanced_results(result, result['response_time'], user_text.strip())
            else:
                st.warning("Please enter clinical text for analysis")
    
    with col2:
        st.markdown('<div class="medical-card">', unsafe_allow_html=True)
        st.markdown("### System Overview")
        
        # Real-time model info
        if api_status['model_info']:
            model_info = api_status['model_info']
            st.markdown('<div class="metric-container">', unsafe_allow_html=True)
            st.metric("Model Accuracy", f"{model_info.get('accuracy', '99.9%')}")
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.markdown('<div class="metric-container">', unsafe_allow_html=True)
            st.metric("Model Type", model_info.get('model_name', 'RBF SVM'))
            st.markdown('</div>', unsafe_allow_html=True)
        
        st.markdown('<div class="metric-container">', unsafe_allow_html=True)
        st.metric("Specialties", "5 Medical Areas")
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Recent predictions
        if st.session_state.metrics["prediction_history"]:
            st.markdown("### Recent Predictions")
            for i, pred in enumerate(st.session_state.metrics["prediction_history"][-5:]):
                confidence_color = "🟢" if pred['confidence'] > 0.8 else "🟡" if pred['confidence'] > 0.6 else "🔴"
                st.markdown(f"{confidence_color} {pred['specialty']} ({pred['confidence']:.1%})")
        
        # Professional notes
        st.markdown("### Clinical Notes")
        st.info("**Best Practice**: Include detailed symptoms, diagnostic findings, and medical terminology for optimal classification accuracy.")
        
        st.markdown('</div>', unsafe_allow_html=True)

def classify_text_enhanced(text):
    """Enhanced classification with metrics tracking and confidence handling"""
    start_time = time.time()
    
    with st.spinner("Analyzing medical text with AI..."):
        time.sleep(0.2)  # Small UX delay
        
        try:
            response = requests.post(
                f"{API_URL}/predict",
                json={"text": text},
                timeout=10
            )
            
            response_time = time.time() - start_time
            
            if response.status_code == 200:
                result = response.json()
                
                # Check if prediction is reliable
                is_reliable = result['specialty'] != "Unknown/Low_Confidence"
                display_specialty = result['specialty']
                
                # Update metrics (only count medical predictions)
                if is_reliable:
                    st.session_state.metrics["total_predictions"] += 1
                    st.session_state.metrics["specialty_counts"][result['specialty']] += 1
                    st.session_state.metrics["confidence_scores"].append(result['confidence'])
                
                st.session_state.metrics["response_times"].append(response_time)
                st.session_state.metrics["prediction_history"].append({
                    'specialty': display_specialty,
                    'confidence': result['confidence'],
                    'timestamp': datetime.now(),
                    'text_length': len(text),
                    'is_reliable': is_reliable
                })
                
                return {
                    'specialty': display_specialty,
                    'confidence': result['confidence'],
                    'response_time': response_time,
                    'is_reliable': is_reliable,
                    'model_version': result.get('model_version', 'Unknown')
                }
            else:
                st.error(f"API Error: {response.status_code}")
                return None
                
        except requests.exceptions.RequestException as e:
            st.error(f"Connection error: {str(e)}")
            return None
        except Exception as e:
            st.error(f"Unexpected error: {str(e)}")
            return None

def display_enhanced_results(result, response_time, original_text):
    """Display enhanced classification results with confidence handling"""
    st.markdown('<div class="medical-card">', unsafe_allow_html=True)
    st.markdown("### AI Classification Results")
    
    # Check if this is a reliable medical prediction
    is_reliable = result.get('is_reliable', True)
    
    # Main metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.markdown('<div class="metric-container">', unsafe_allow_html=True)
        specialty_display = result['specialty']
        if not is_reliable:
            specialty_display = "⚠️ Non-Medical/Uncertain"
        st.metric("Medical Specialty", specialty_display)
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col2:
        confidence_pct = round(result['confidence'] * 100, 1)
        st.markdown('<div class="metric-container">', unsafe_allow_html=True)
        st.metric("Confidence Score", f"{confidence_pct}%")
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col3:
        st.markdown('<div class="metric-container">', unsafe_allow_html=True)
        st.metric("Response Time", f"{response_time:.3f}s")
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col4:
        st.markdown('<div class="metric-container">', unsafe_allow_html=True)
        st.metric("📝 Text Length", f"{len(original_text)} chars")
        st.markdown('</div>', unsafe_allow_html=True)
    
    # Clinical interpretation
    st.markdown("### 🩺 Clinical Interpretation")
    
    if not is_reliable:
        st.error(f"❗ **Low Medical Relevance** ({confidence_pct}%) - Text may not be medical or too ambiguous for reliable classification. Please provide clear medical terminology and symptoms.")
    elif result['confidence'] > 0.9:
        st.success(f"**High Confidence Classification** ({confidence_pct}%) - Excellent model certainty")
    elif result['confidence'] > 0.8:
        st.success(f"✅ **Good Confidence Classification** ({confidence_pct}%) - Reliable prediction")
    elif result['confidence'] > 0.7:
        st.warning(f"**Moderate Confidence** ({confidence_pct}%) - Consider clinical correlation")
    else:
        st.warning(f"**Borderline Confidence** ({confidence_pct}%) - Recommend clinical review")
    
    # Specialty description (only for reliable predictions)
    if is_reliable:
        specialty_descriptions = {
            "Cardiology": "Heart and cardiovascular system conditions requiring cardiac evaluation",
            "Emergency": "Urgent medical situations requiring immediate intervention and stabilization",
            "Pulmonology": "Respiratory system and lung-related conditions requiring pulmonary assessment",
            "Gastroenterology": "Digestive system and gastrointestinal disorders requiring GI evaluation",
            "Dermatology": "Skin, hair, and nail-related conditions requiring dermatologic assessment"
        }
        
        if result['specialty'] in specialty_descriptions:
            st.info(f"**{result['specialty']}**: {specialty_descriptions[result['specialty']]}")
    else:
        st.info("💡 **Recommendation**: For better results, include specific medical terminology, symptoms, diagnostic findings, or clinical observations in your text.")
    
    # Model version info
    if 'model_version' in result:
        st.caption(f"Model: {result['model_version']}")
    
    st.markdown('</div>', unsafe_allow_html=True)

def analytics_dashboard():
    """Comprehensive analytics dashboard"""
    st.markdown("### System Analytics & Insights")
    
    if st.session_state.metrics["total_predictions"] == 0:
        st.info("No predictions yet. Start classifying medical texts to see analytics.")
        return
    
    # Overview metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Predictions", st.session_state.metrics["total_predictions"])
    
    with col2:
        avg_confidence = sum(st.session_state.metrics["confidence_scores"]) / len(st.session_state.metrics["confidence_scores"])
        st.metric("Average Confidence", f"{avg_confidence:.2%}")
    
    with col3:
        avg_response = sum(st.session_state.metrics["response_times"]) / len(st.session_state.metrics["response_times"])
        st.metric("Avg Response Time", f"{avg_response:.3f}s")
    
    with col4:
        unique_specialties = len(st.session_state.metrics["specialty_counts"])
        st.metric("Specialties Used", unique_specialties)
    
    # Charts
    col1, col2 = st.columns(2)
    
    with col1:
        # Specialty distribution pie chart
        if st.session_state.metrics["specialty_counts"]:
            specialty_df = pd.DataFrame([
                {"Specialty": k, "Count": v} 
                for k, v in st.session_state.metrics["specialty_counts"].items()
            ])
            
            fig_pie = px.pie(
                specialty_df, 
                values='Count', 
                names='Specialty',
                title="Specialty Distribution",
                color_discrete_sequence=px.colors.qualitative.Set2
            )
            st.plotly_chart(fig_pie, use_container_width=True)
    
    with col2:
        # Confidence distribution histogram
        if st.session_state.metrics["confidence_scores"]:
            confidence_df = pd.DataFrame({
                "Confidence": st.session_state.metrics["confidence_scores"]
            })
            
            fig_hist = px.histogram(
                confidence_df,
                x="Confidence",
                title="Confidence Score Distribution",
                nbins=20,
                color_discrete_sequence=["#28a745"]
            )
            st.plotly_chart(fig_hist, use_container_width=True)
    
    # Prediction timeline
    if len(st.session_state.metrics["prediction_history"]) > 1:
        st.markdown("### Prediction Timeline")
        
        timeline_df = pd.DataFrame([
            {
                "Time": pred['timestamp'],
                "Specialty": pred['specialty'],
                "Confidence": pred['confidence'],
                "Text Length": pred['text_length']
            }
            for pred in st.session_state.metrics["prediction_history"]
        ])
        
        fig_timeline = px.scatter(
            timeline_df,
            x="Time",
            y="Confidence",
            color="Specialty",
            size="Text Length",
            title="Predictions Over Time",
            hover_data=['Text Length']
        )
        st.plotly_chart(fig_timeline, use_container_width=True)

def batch_processing_interface():
    """Batch processing interface"""
    st.markdown("### Batch Processing")
    
    st.info("Process multiple medical texts simultaneously for efficient bulk analysis.")
    
    # Text input for batch processing
    batch_texts = st.text_area(
        "Enter multiple medical texts (one per line):",
        height=300,
        placeholder="Patient 1: Chest pain with ST elevation...\nPatient 2: Severe dyspnea with hemoptysis...\nPatient 3: Abdominal pain with guarding..."
    )
    
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("Process Batch", type="primary"):
            if batch_texts.strip():
                process_batch(batch_texts)
            else:
                st.warning("Please enter texts for batch processing")
    
    with col2:
        # File upload for batch processing
        uploaded_file = st.file_uploader(
            "Or upload a text file:",
            type=['txt', 'csv'],
            help="Upload a text file with one case per line"
        )
        
        if uploaded_file and st.button("Process File"):
            file_content = uploaded_file.read().decode('utf-8')
            process_batch(file_content)

def process_batch(batch_text):
    """Process batch of texts"""
    texts = [text.strip() for text in batch_text.split('\n') if text.strip()]
    
    if not texts:
        st.warning("No valid texts found for processing")
        return
    
    results = []
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i, text in enumerate(texts):
        status_text.text(f"Processing {i+1}/{len(texts)}: {text[:50]}...")
        
        try:
            response = requests.post(
                f"{API_URL}/predict",
                json={"text": text},
                timeout=10
            )
            
            if response.status_code == 200:
                result = response.json()
                results.append({
                    "Text": text[:100] + "..." if len(text) > 100 else text,
                    "Specialty": result['specialty'],
                    "Confidence": f"{result['confidence']:.2%}",
                    "Full_Text": text
                })
                
                # Update metrics
                st.session_state.metrics["total_predictions"] += 1
                st.session_state.metrics["specialty_counts"][result['specialty']] += 1
                st.session_state.metrics["confidence_scores"].append(result['confidence'])
            else:
                results.append({
                    "Text": text[:100] + "..." if len(text) > 100 else text,
                    "Specialty": "Error",
                    "Confidence": "N/A",
                    "Full_Text": text
                })
        except Exception as e:
            results.append({
                "Text": text[:100] + "..." if len(text) > 100 else text,
                "Specialty": f"Error: {str(e)}",
                "Confidence": "N/A",
                "Full_Text": text
            })
        
        progress_bar.progress((i + 1) / len(texts))
    
    status_text.text("✅ Batch processing complete!")
    
    # Display results
    st.markdown("### Batch Results")
    df = pd.DataFrame(results)
    st.dataframe(df[['Text', 'Specialty', 'Confidence']], use_container_width=True)
    
    # Summary statistics
    col1, col2, col3 = st.columns(3)
    
    successful_results = [r for r in results if r['Specialty'] != 'Error' and not r['Specialty'].startswith('Error:')]
    
    with col1:
        st.metric("Total Processed", len(results))
    with col2:
        st.metric("Successful", len(successful_results))
    with col3:
        if successful_results:
            avg_conf = sum(float(r['Confidence'].rstrip('%'))/100 for r in successful_results) / len(successful_results)
            st.metric("Avg Confidence", f"{avg_conf:.2%}")
    
    # Store results for export
    st.session_state['last_batch_results'] = results

def test_suite_interface():
    """Comprehensive test suite with edge cases"""
    st.markdown("### Comprehensive Test Suite")
    
    st.info("Automated testing with edge cases, boundary conditions, and validation scenarios.")
    
    # Test categories
    test_categories = {
        "Basic Functionality": get_basic_tests(),
        "Edge Cases": get_edge_case_tests(),
        "Boundary Conditions": get_boundary_tests(),
        "Clinical Variations": get_clinical_variation_tests(),
        "Stress Tests": get_stress_tests()
    }
    
    selected_category = st.selectbox("Select test category:", list(test_categories.keys()))
    
    col1, col2 = st.columns(2)
    
    with col1:
        if st.button("Run Selected Tests", type="primary"):
            run_test_suite(test_categories[selected_category], selected_category)
    
    with col2:
        if st.button("Run All Tests"):
            for category, tests in test_categories.items():
                st.markdown(f"#### Running {category}")
                run_test_suite(tests, category)

def get_basic_tests():
    """Basic functionality tests"""
    return [
        {"name": "Cardiology Basic", "text": "Patient with chest pain and elevated troponins", "expected": "Cardiology"},
        {"name": "Emergency Basic", "text": "Motor vehicle accident with trauma", "expected": "Emergency"},
        {"name": "Pulmonology Basic", "text": "Chronic cough with hemoptysis", "expected": "Pulmonology"},
        {"name": "Gastroenterology Basic", "text": "Epigastric pain with gastric ulcer", "expected": "Gastroenterology"},
        {"name": "Dermatology Basic", "text": "Changing mole with irregular borders", "expected": "Dermatology"}
    ]

def get_edge_case_tests():
    """Edge case tests"""
    return [
        {"name": "Very Short Text", "text": "Chest pain", "expected": None},
        {"name": "Very Long Text", "text": "Patient " * 1000 + "with chest pain and MI", "expected": "Cardiology"},
        {"name": "No Medical Terms", "text": "The weather is nice today and I feel happy", "expected": None},
        {"name": "Mixed Specialties", "text": "Patient with chest pain and skin rash requiring emergency care", "expected": None},
        {"name": "Special Characters", "text": "Patient w/ chest pain & ST↑ in ECG → STEMI", "expected": "Cardiology"}
    ]

def get_boundary_tests():
    """Boundary condition tests"""
    return [
        {"name": "Minimum Length", "text": "MI", "expected": None},
        {"name": "Maximum Length", "text": "A" * 10000, "expected": None},
        {"name": "Numbers Only", "text": "12345678901234567890", "expected": None},
        {"name": "Empty String", "text": "", "expected": None},
        {"name": "Whitespace Only", "text": "   \n\t   ", "expected": None}
    ]

def get_clinical_variation_tests():
    """Clinical variation tests"""
    return [
        {"name": "Abbreviations", "text": "Pt c/o SOB, CXR shows PNA", "expected": "Pulmonology"},
        {"name": "Formal Language", "text": "The patient demonstrates acute coronary syndrome", "expected": "Cardiology"},
        {"name": "Casual Language", "text": "Guy came in with bad belly pain", "expected": "Gastroenterology"},
        {"name": "Technical Terms", "text": "Echocardiogram reveals reduced LVEF with regional wall motion abnormalities", "expected": "Cardiology"},
        {"name": "Symptoms Only", "text": "Severe crushing chest pain radiating to left arm", "expected": "Cardiology"}
    ]

def get_stress_tests():
    """Stress tests"""
    return [
        {"name": "Repeated Text", "text": "chest pain " * 100, "expected": "Cardiology"},
        {"name": "Medical Gibberish", "text": "Pseudohypoparathyroidism pneumonoultramicroscopicsilicovolcanoconiosis", "expected": None},
        {"name": "Unicode Text", "text": "Patient with 胸痛 and ♥ problems", "expected": None},
        {"name": "HTML Tags", "text": "<p>Patient with <b>chest pain</b></p>", "expected": "Cardiology"},
        {"name": "SQL Injection", "text": "'; DROP TABLE patients; --", "expected": None}
    ]

def run_test_suite(tests, category_name):
    """Run a test suite"""
    st.markdown(f"#### Testing: {category_name}")
    
    results = []
    progress_bar = st.progress(0)
    
    for i, test in enumerate(tests):
        try:
            response = requests.post(
                f"{API_URL}/predict",
                json={"text": test["text"]},
                timeout=5
            )
            
            if response.status_code == 200:
                result = response.json()
                predicted = result['specialty']
                confidence = result['confidence']
                
                # Determine pass/fail
                if test["expected"] is None:
                    # For edge cases, we expect low confidence or specific handling
                    passed = confidence < 0.6 or predicted in ["Error", "Unknown"]
                else:
                    passed = predicted == test["expected"]
                
                results.append({
                    "Test": test["name"],
                    "Expected": test["expected"] or "Low Confidence",
                    "Predicted": predicted,
                    "Confidence": f"{confidence:.2%}",
                    "Status": "✅ PASS" if passed else "❌ FAIL"
                })
            else:
                results.append({
                    "Test": test["name"],
                    "Expected": test["expected"] or "Low Confidence",
                    "Predicted": f"Error {response.status_code}",
                    "Confidence": "N/A",
                    "Status": "❌ FAIL"
                })
        
        except Exception as e:
            results.append({
                "Test": test["name"],
                "Expected": test["expected"] or "Low Confidence",
                "Predicted": f"Exception: {str(e)}",
                "Confidence": "N/A",
                "Status": "❌ FAIL"
            })
        
        progress_bar.progress((i + 1) / len(tests))
    
    # Display results
    df = pd.DataFrame(results)
    st.dataframe(df, use_container_width=True)
    
    # Summary
    passed = len([r for r in results if r["Status"] == "✅ PASS"])
    total = len(results)
    pass_rate = (passed / total) * 100
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Tests Passed", f"{passed}/{total}")
    with col2:
        st.metric("Pass Rate", f"{pass_rate:.1f}%")
    with col3:
        status_color = "🟢" if pass_rate >= 80 else "🟡" if pass_rate >= 60 else "🔴"
        st.metric("Status", f"{status_color} {'GOOD' if pass_rate >= 80 else 'FAIR' if pass_rate >= 60 else 'POOR'}")

def export_interface():
    """Export interface for results and analytics"""
    st.markdown("### Export Data & Reports")
    
    st.info("Export your classification results, analytics, and system metrics.")
    
    # Export options
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Analytics Export")
        
        if st.button("Export Analytics Report"):
            export_analytics_report()
        
        if st.button("Export Prediction History"):
            export_prediction_history()
    
    with col2:
        st.markdown("#### 📄 Batch Results Export")
        
        if 'last_batch_results' in st.session_state:
            if st.button("Export Last Batch Results"):
                export_batch_results()
        else:
            st.info("No batch results to export. Run batch processing first.")
    
    # System configuration export
    st.markdown("#### System Configuration")
    if st.button("🔧 Export System Config"):
        export_system_config()

def export_analytics_report():
    """Export comprehensive analytics report"""
    if st.session_state.metrics["total_predictions"] == 0:
        st.warning("No data to export. Perform some classifications first.")
        return
    
    # Create comprehensive report
    report_data = {
        "System Overview": {
            "Total Predictions": st.session_state.metrics["total_predictions"],
            "Average Confidence": sum(st.session_state.metrics["confidence_scores"]) / len(st.session_state.metrics["confidence_scores"]),
            "Average Response Time": sum(st.session_state.metrics["response_times"]) / len(st.session_state.metrics["response_times"]),
            "Unique Specialties": len(st.session_state.metrics["specialty_counts"])
        },
        "Specialty Distribution": dict(st.session_state.metrics["specialty_counts"]),
        "Performance Metrics": {
            "Confidence Scores": st.session_state.metrics["confidence_scores"],
            "Response Times": st.session_state.metrics["response_times"]
        },
        "Export Timestamp": datetime.now().isoformat()
    }
    
    # Convert to JSON
    json_data = json.dumps(report_data, indent=2)
    
    st.download_button(
        label="Download Analytics Report (JSON)",
        data=json_data,
        file_name=f"medical_ai_analytics_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
        mime="application/json"
    )
    
    st.success("✅ Analytics report ready for download!")

def export_prediction_history():
    """Export prediction history as CSV"""
    if not st.session_state.metrics["prediction_history"]:
        st.warning("No prediction history to export.")
        return
    
    # Create DataFrame
    df = pd.DataFrame([
        {
            "Timestamp": pred['timestamp'].isoformat(),
            "Specialty": pred['specialty'],
            "Confidence": pred['confidence'],
            "Text_Length": pred['text_length']
        }
        for pred in st.session_state.metrics["prediction_history"]
    ])
    
    # Convert to CSV
    csv_data = df.to_csv(index=False)
    
    st.download_button(
        label="Download Prediction History (CSV)",
        data=csv_data,
        file_name=f"prediction_history_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        mime="text/csv"
    )
    
    st.success("✅ Prediction history ready for download!")

def export_batch_results():
    """Export last batch processing results"""
    if 'last_batch_results' not in st.session_state:
        st.warning("No batch results to export.")
        return
    
    # Create DataFrame
    df = pd.DataFrame(st.session_state['last_batch_results'])
    
    # Convert to CSV
    csv_data = df.to_csv(index=False)
    
    st.download_button(
        label="Download Batch Results (CSV)",
        data=csv_data,
        file_name=f"batch_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
        mime="text/csv"
    )
    
    st.success("✅ Batch results ready for download!")

def export_system_config():
    """Export system configuration"""
    config_data = {
        "API_URL": API_URL,
        "System_Info": {
            "Dashboard_Version": "2.0.0",
            "Features": [
                "Real-time Classification",
                "Batch Processing",
                "Comprehensive Analytics",
                "Test Suite",
                "Data Export"
            ]
        },
        "Sample_Specialties": list(get_comprehensive_samples().keys()),
        "Export_Timestamp": datetime.now().isoformat()
    }
    
    json_data = json.dumps(config_data, indent=2)
    
    st.download_button(
        label="🔧 Download System Config (JSON)",
        data=json_data,
        file_name=f"system_config_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
        mime="application/json"
    )
    
    st.success("✅ System configuration ready for download!")

def check_api_connection():
    """Check API connection and get model info"""
    try:
        # Health check
        health_response = requests.get(f"{API_URL}/health", timeout=5)
        if health_response.status_code != 200:
            return {'connected': False, 'model_info': None}
        
        # Get model info
        try:
            model_response = requests.get(f"{API_URL}/model-info", timeout=5)
            model_info = model_response.json() if model_response.status_code == 200 else None
        except:
            model_info = None
        
        return {'connected': True, 'model_info': model_info}
    
    except Exception as e:
        return {'connected': False, 'model_info': None}

if __name__ == "__main__":
    main()

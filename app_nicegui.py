#!/usr/bin/env python3
"""BiOmics Landing Page - NiceGUI Version"""

from nicegui import ui, app
from pathlib import Path

# 导入 app_nicegui_chat 模块，注册 /platform 路由
import app_nicegui_chat

# 静态文件路径
app_path = Path(__file__).parent
app.add_static_files('/static', app_path)

# 全局样式
ui.add_head_html('''
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap" rel="stylesheet">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
<style>
    body { font-family: 'Inter', sans-serif; margin: 0; padding: 0; }
    .nicegui-content { padding: 0 !important; }
    .full-width-section { width: 100vw !important; margin-left: calc(-50vw + 50%) !important; padding-left: calc(50vw - 50%) !important; padding-right: calc(50vw - 50%) !important; }
    .hero-gradient { background: linear-gradient(135deg, #0f172a 0%, #1e3a8a 100%); }
    .glass-card { background: rgba(255, 255, 255, 0.05); backdrop-filter: blur(10px); border: 1px solid rgba(255, 255, 255, 0.1); }
    .nav-link { color: #475569; transition: color 0.2s; }
    .nav-link:hover { color: #2563eb; }
    .capability-icon { transition: all 0.3s; }
    .capability-card-cyan:hover .capability-icon { background-color: #0891b2 !important; color: white !important; }
    .capability-card-emerald:hover .capability-icon { background-color: #059669 !important; color: white !important; }
    .capability-card-amber:hover .capability-icon { background-color: #d97706 !important; color: white !important; }
</style>
''', shared=True)

def create_landing_page():
    # ==================== 导航栏 ====================
    with ui.header().classes('bg-white/80 backdrop-blur-md border-b border-slate-200 flex items-center').style('height: 80px; padding-top: 0; padding-bottom: 0;'):
        with ui.row().classes('max-w-7xl mx-auto w-full pl-2 pr-6 items-center justify-between'):
            # Logo
            with ui.row().classes('items-center gap-2').style('margin-left: -8px;'):
                ui.html('''<svg xmlns="http://www.w3.org/2000/svg" viewBox="270 410 90 25" style="height: 70px; width: auto;">
                  <style type="text/css">.st0{fill:none;}.st1{fill:#476179;font-weight:bold;}.st2{font-family:'Arial-BoldMT';}.st3{font-size:9px;}.st4{fill:#FFFFFF;stroke:#476179;stroke-width:0.5;stroke-miterlimit:10;}.st5{fill:none;stroke:#476179;stroke-width:0.5;stroke-miterlimit:10;}.st6{opacity:0.5;fill:#455D7D;}.st7{fill:#476179;stroke:#476179;stroke-width:0.5;stroke-miterlimit:10;}.st9{fill:#FFFFFF;}.st10{fill:none;stroke:#476179;stroke-width:0.1294;stroke-miterlimit:10;}.st11{fill:none;stroke:#FFFFFF;stroke-width:0.1294;stroke-miterlimit:10;}.st12{fill:#476179;stroke:#476179;stroke-width:0.0863;stroke-miterlimit:10;}.st13{fill:#FFFFFF;stroke:#FFFFFF;stroke-width:0.0863;stroke-miterlimit:10;}.st14{fill:#FFFFFF;stroke:#7A75AB;stroke-width:0.0863;stroke-miterlimit:10;}.st15{fill:#476179;stroke:#476179;stroke-width:0.0863;stroke-miterlimit:10;}.st16{fill:#FFFFFF;stroke:#476179;stroke-width:0.0789;stroke-miterlimit:10;}.st17{fill:none;stroke:#FFFFFF;stroke-width:0.0731;stroke-miterlimit:10;}.st20{fill:#FFFFFF;stroke:#476179;stroke-width:0.0863;stroke-miterlimit:10;}</style>
                  <path class="st0" d="M348.25,429.42h-70.83c-3.48,0-6.3-2.82-6.3-6.3v-3.72c0-3.48,2.82-6.3,6.3-6.3h70.83c3.48,0,6.3,2.82,6.3,6.3v3.72C354.56,426.59,351.73,429.42,348.25,429.42z"/>
                  <text transform="matrix(1 0 0 1 275.389 424.6138)" class="st1 st2 st3">Bi</text><text transform="matrix(1 0 0 1 298.3929 424.6138)" class="st1 st2 st3">mics</text>
                  <g><path class="st4" d="M288.91,417.49c0-1.32,1.03-2.41,2.34-2.48c0.05,0,0.1,0,0.15,0s0.1,0,0.14,0c2.68,0.07,4.83,2.27,4.83,4.97s-2.15,4.89-4.83,4.97c1.3-0.07,2.34-1.16,2.34-2.48c0-1.37-1.11-2.49-2.48-2.49C290.02,419.97,288.91,418.86,288.91,417.49"/><path class="st5" d="M291.23,415.01c0.05,0,0.1,0,0.14,0c0.05,0,0.1,0,0.15,0c2.68,0.07,4.83,2.27,4.83,4.97c0,2.7-2.15,4.89-4.83,4.97"/><path class="st7" d="M293.86,422.46c0,1.32-1.03,2.41-2.34,2.48c-0.05,0-0.1,0-0.15,0s-0.1,0-0.14,0c-2.68-0.07-4.83-2.27-4.83-4.97s2.15-4.89,4.83-4.97c-1.3,0.07-2.34,1.16-2.34,2.48c0,1.37,1.11,2.49,2.48,2.49C292.75,419.98,293.86,421.09,293.86,422.46"/><g><line class="st10" x1="288.93" y1="418.28" x2="290.37" y2="419.57"/><line class="st10" x1="292.23" y1="419.18" x2="290.37" y2="419.57"/><line class="st10" x1="294.13" y1="420.76" x2="292.23" y2="419.18"/><line class="st10" x1="292.09" y1="416.58" x2="292.23" y2="419.18"/><line class="st10" x1="288.93" y1="418.28" x2="292.09" y2="416.58"/><line class="st11" x1="290.37" y1="419.57" x2="289.77" y2="422.24"/><line class="st11" x1="294.13" y1="420.76" x2="289.77" y2="422.24"/><path class="st12" d="M293.87,421.08c0.17,0.14,0.43,0.12,0.57-0.06c0.14-0.17,0.12-0.43-0.06-0.57c-0.17-0.14-0.43-0.12-0.57,0.06C293.67,420.68,293.7,420.93,293.87,421.08z"/><path class="st12" d="M291.83,416.89c0.17,0.14,0.43,0.12,0.57-0.06c0.14-0.17,0.12-0.43-0.06-0.57c-0.17-0.14-0.43-0.12-0.57,0.06C291.64,416.49,291.66,416.75,291.83,416.89z"/><path class="st13" d="M288.67,418.61c0.17,0.14,0.43,0.12,0.57-0.06c0.14-0.17,0.12-0.43-0.06-0.57s-0.43-0.12-0.57,0.06C288.48,418.21,288.5,418.47,288.67,418.61z"/><path class="st14" d="M291.97,419.49c0.17,0.14,0.43,0.12,0.57-0.06c0.14-0.17,0.12-0.43-0.06-0.57s-0.43-0.12-0.57,0.06C291.77,419.09,291.8,419.35,291.97,419.49z"/><path class="st15" d="M290.11,419.88c0.17,0.14,0.43,0.12,0.57-0.06c0.14-0.17,0.12-0.43-0.06-0.57c-0.17-0.14-0.43-0.12-0.57,0.06C289.91,419.48,289.94,419.74,290.11,419.88z"/><path class="st13" d="M289.51,422.56c0.17,0.14,0.43,0.12,0.57-0.06c0.14-0.17,0.12-0.43-0.06-0.57s-0.43-0.12-0.57,0.06C289.31,422.16,289.34,422.41,289.51,422.56z"/><path class="st16" d="M288.7,418.58c0.16,0.13,0.39,0.11,0.53-0.05c0.13-0.16,0.11-0.39-0.05-0.53c-0.16-0.13-0.39-0.11-0.53,0.05C288.51,418.22,288.54,418.45,288.7,418.58z"/><path class="st16" d="M289.54,422.53c0.16,0.13,0.39,0.11,0.53-0.05c0.13-0.16,0.11-0.39-0.05-0.53c-0.16-0.13-0.39-0.11-0.53,0.05C289.35,422.16,289.38,422.4,289.54,422.53z"/><path class="st17" d="M293.91,421.03c0.15,0.12,0.37,0.1,0.49-0.05c0.12-0.15,0.1-0.37-0.05-0.49c-0.15-0.12-0.37-0.1-0.49,0.05C293.74,420.69,293.76,420.91,293.91,421.03z"/><path class="st17" d="M291.87,416.84c0.15,0.12,0.37,0.1,0.49-0.05c0.12-0.15,0.1-0.37-0.05-0.49s-0.37-0.1-0.49,0.05C291.71,416.5,291.73,416.72,291.87,416.84z"/></g></g>
                </svg>''', sanitize=False)            
            # 按钮
            ui.button('Try it now', on_click=lambda: ui.navigate.to('/platform')).props('no-caps').classes('bg-black text-white px-6 py-2.5 rounded-full font-semibold hover:bg-slate-800 shadow-md')

    # ==================== 1. Hero Section (Intro + Quick Start) ====================
    with ui.element('section').props('id="paradigm"').classes('text-white py-24 px-6 w-screen').style('margin-left: calc(-50vw + 50%); box-sizing: border-box; position: relative; overflow: hidden; background: linear-gradient(135deg, #0a1929 0%, #1a3a52 50%, #0f2942 100%);'):
        # 添加生物信息学科技风背景
        ui.html('''
            <div style="position: absolute; top: 0; left: 0; width: 100%; height: 100%; pointer-events: none; overflow: hidden;">
                <!-- 光晕效果 -->
                <div style="position: absolute; width: 500px; height: 500px; top: -10%; left: -10%; 
                     background: radial-gradient(circle, #60a5fa 0%, transparent 70%); 
                     border-radius: 50%; filter: blur(80px); opacity: 0.15;
                     animation: float-glow1 20s ease-in-out infinite;"></div>
                <div style="position: absolute; width: 400px; height: 400px; top: 30%; right: -5%; 
                     background: radial-gradient(circle, #22d3ee 0%, transparent 70%); 
                     border-radius: 50%; filter: blur(80px); opacity: 0.15;
                     animation: float-glow2 20s ease-in-out infinite; animation-delay: 3s;"></div>
                <div style="position: absolute; width: 350px; height: 350px; bottom: -5%; left: 20%; 
                     background: radial-gradient(circle, #34d399 0%, transparent 70%); 
                     border-radius: 50%; filter: blur(80px); opacity: 0.15;
                     animation: float-glow3 20s ease-in-out infinite; animation-delay: 6s;"></div>
                
                <!-- DNA 图标 -->
                <i class="fas fa-dna" style="position: absolute; top: 15%; left: 8%; font-size: 80px; 
                   color: #60a5fa; opacity: 0.12; animation: dna-rotate1 20s linear infinite;"></i>
                <i class="fas fa-dna" style="position: absolute; bottom: 20%; right: 12%; font-size: 70px; 
                   color: #34d399; opacity: 0.12; animation: dna-rotate2 20s linear infinite; animation-delay: 5s;"></i>
                
                <!-- 基因序列代码流 -->
                <div style="position: absolute; top: 20%; right: 5%; font-family: 'Courier New', monospace;
                     font-size: 12px; color: #34d399; opacity: 0.4; line-height: 1.8;
                     animation: fade-sequence 5s ease-in-out infinite;">
                    ATCG TAGC GCTA<br>
                    CGAT ATGC TACG<br>
                    TAGC CGTA GCTA<br>
                    ATCG TACG GCAT<br>
                    GCTA ATCG TAGC
                </div>
                
                <!-- 分子结构 -->
                <svg style="position: absolute; bottom: 15%; left: 12%; width: 150px; height: 150px; opacity: 0.12;">
                    <line x1="50" y1="40" x2="75" y2="75" stroke="#60a5fa" stroke-width="1" opacity="0.5"/>
                    <line x1="75" y1="75" x2="100" y2="50" stroke="#60a5fa" stroke-width="1" opacity="0.5"/>
                    <line x1="75" y1="75" x2="90" y2="110" stroke="#60a5fa" stroke-width="1" opacity="0.5"/>
                    <line x1="75" y1="75" x2="40" y2="100" stroke="#60a5fa" stroke-width="1" opacity="0.5"/>
                    <line x1="40" y1="100" x2="60" y2="130" stroke="#60a5fa" stroke-width="1" opacity="0.5"/>
                    <line x1="90" y1="110" x2="110" y2="130" stroke="#60a5fa" stroke-width="1" opacity="0.5"/>
                    <circle cx="50" cy="40" r="4" fill="#60a5fa" style="animation: pulse-atom 2s ease-in-out infinite;" />
                    <circle cx="75" cy="75" r="5" fill="#34d399" style="animation: pulse-atom 2s ease-in-out infinite; animation-delay: 0.3s" />
                    <circle cx="100" cy="50" r="4" fill="#22d3ee" style="animation: pulse-atom 2s ease-in-out infinite; animation-delay: 0.6s" />
                    <circle cx="90" cy="110" r="4" fill="#fbbf24" style="animation: pulse-atom 2s ease-in-out infinite; animation-delay: 0.9s" />
                    <circle cx="40" cy="100" r="4" fill="#60a5fa" style="animation: pulse-atom 2s ease-in-out infinite; animation-delay: 1.2s" />
                    <circle cx="60" cy="130" r="3" fill="#34d399" style="animation: pulse-atom 2s ease-in-out infinite; animation-delay: 1.5s" />
                    <circle cx="110" cy="130" r="3" fill="#22d3ee" style="animation: pulse-atom 2s ease-in-out infinite; animation-delay: 1.8s" />
                </svg>
                
                <!-- 六边形网络 -->
                <svg style="position: absolute; bottom: 20%; right: 8%; width: 180px; height: 180px; opacity: 0.08;">
                    <polygon points="90,30 120,50 120,90 90,110 60,90 60,50" stroke="#60a5fa" stroke-width="1" fill="none" style="animation: glow-hex 3s ease-in-out infinite;" />
                    <polygon points="30,30 50,40 50,60 30,70 10,60 10,40" stroke="#60a5fa" stroke-width="1" fill="none" style="animation: glow-hex 3s ease-in-out infinite; animation-delay: 0.5s" />
                    <polygon points="150,30 170,40 170,60 150,70 130,60 130,40" stroke="#60a5fa" stroke-width="1" fill="none" style="animation: glow-hex 3s ease-in-out infinite; animation-delay: 1s" />
                    <polygon points="30,110 50,120 50,140 30,150 10,140 10,120" stroke="#60a5fa" stroke-width="1" fill="none" style="animation: glow-hex 3s ease-in-out infinite; animation-delay: 1.5s" />
                    <polygon points="150,110 170,120 170,140 150,150 130,140 130,120" stroke="#60a5fa" stroke-width="1" fill="none" style="animation: glow-hex 3s ease-in-out infinite; animation-delay: 2s" />
                    <line x1="60" y1="50" x2="50" y2="60" stroke="#60a5fa" stroke-width="0.5" opacity="0.4"/>
                    <line x1="120" y1="50" x2="130" y2="60" stroke="#60a5fa" stroke-width="0.5" opacity="0.4"/>
                    <line x1="60" y1="90" x2="50" y2="120" stroke="#60a5fa" stroke-width="0.5" opacity="0.4"/>
                    <line x1="120" y1="90" x2="130" y2="120" stroke="#60a5fa" stroke-width="0.5" opacity="0.4"/>
                </svg>
                
                <!-- 数据节点网络 -->
                <svg style="position: absolute; top: 30%; left: 25%; width: 250px; height: 200px; opacity: 0.1;">
                    <line x1="60" y1="50" x2="120" y2="80" stroke="#60a5fa" stroke-width="0.5" opacity="0.6"/>
                    <line x1="120" y1="80" x2="180" y2="60" stroke="#60a5fa" stroke-width="0.5" opacity="0.6"/>
                    <line x1="120" y1="80" x2="150" y2="140" stroke="#60a5fa" stroke-width="0.5" opacity="0.6"/>
                    <line x1="60" y1="50" x2="90" y2="120" stroke="#60a5fa" stroke-width="0.5" opacity="0.6"/>
                    <line x1="90" y1="120" x2="150" y2="140" stroke="#60a5fa" stroke-width="0.5" opacity="0.6"/>
                    <circle cx="60" cy="50" r="3" fill="#34d399" style="animation: node-pulse 3s ease-in-out infinite;" />
                    <circle cx="120" cy="80" r="3.5" fill="#34d399" style="animation: node-pulse 3s ease-in-out infinite; animation-delay: 0.4s" />
                    <circle cx="180" cy="60" r="2.5" fill="#34d399" style="animation: node-pulse 3s ease-in-out infinite; animation-delay: 0.8s" />
                    <circle cx="150" cy="140" r="3" fill="#34d399" style="animation: node-pulse 3s ease-in-out infinite; animation-delay: 1.2s" />
                    <circle cx="90" cy="120" r="2.5" fill="#34d399" style="animation: node-pulse 3s ease-in-out infinite; animation-delay: 1.6s" />
                </svg>
                
                <style>
                    @keyframes float-glow1 {
                        0%, 100% { transform: translate(0, 0) scale(1); }
                        33% { transform: translate(30px, -30px) scale(1.1); }
                        66% { transform: translate(-20px, 20px) scale(0.9); }
                    }
                    @keyframes float-glow2 {
                        0%, 100% { transform: translate(0, 0) scale(1); }
                        33% { transform: translate(30px, -30px) scale(1.1); }
                        66% { transform: translate(-20px, 20px) scale(0.9); }
                    }
                    @keyframes float-glow3 {
                        0%, 100% { transform: translate(0, 0) scale(1); }
                        33% { transform: translate(30px, -30px) scale(1.1); }
                        66% { transform: translate(-20px, 20px) scale(0.9); }
                    }
                    @keyframes dna-rotate1 {
                        0%, 100% { opacity: 0.08; transform: rotate(0deg) scale(1); }
                        50% { opacity: 0.15; transform: rotate(180deg) scale(1.1); }
                    }
                    @keyframes dna-rotate2 {
                        0%, 100% { opacity: 0.08; transform: rotate(0deg) scale(1); }
                        50% { opacity: 0.15; transform: rotate(180deg) scale(1.1); }
                    }
                    @keyframes fade-sequence {
                        0%, 100% { opacity: 0.2; }
                        50% { opacity: 0.6; }
                    }
                    @keyframes pulse-atom {
                        0%, 100% { r: 4; opacity: 0.6; }
                        50% { r: 6; opacity: 1; }
                    }
                    @keyframes glow-hex {
                        0%, 100% { opacity: 0.3; }
                        50% { opacity: 0.8; stroke-width: 2; }
                    }
                    @keyframes node-pulse {
                        0%, 100% { r: 2; opacity: 0.4; }
                        50% { r: 3.5; opacity: 1; }
                    }
                </style>
            </div>
        ''', sanitize=False)
        
        with ui.column().classes('max-w-7xl mx-auto text-center items-center w-full').style('position: relative; z-index: 1;'):
            # 标签
            ui.element('div').classes('inline-block px-4 py-1.5 mb-6 rounded-full bg-blue-500/20 border border-blue-400/30 text-blue-300 text-sm font-semibold tracking-wide uppercase').props('innerHTML="Future of Bioinformatics"')
            
            # 标题
            ui.html('''
                <h1 class="text-5xl md:text-7xl font-bold mb-8 leading-tight">
                    From Disjointed to <span class="text-blue-400">Seamless</span><br>Analysis-Interpretation
                </h1>
            ''', sanitize=False)
            
            # 描述
            ui.label('Empowering researchers with a Knowledge-Aware paradigm that bridges the gap between raw omics data and biological discovery through autonomous AI interpretation.').classes('text-xl text-slate-300 max-w-3xl mx-auto mb-12 leading-relaxed') 
            # Quick Start 按钮组
            with ui.row().classes('grid lg:grid-cols-2 gap-8 mt-16 w-full max-w-4xl'):
                # 第一个按钮组：工作相关文件链接
                with ui.column().classes('glass-card p-8 rounded-2xl text-left'):
                    ui.label('Resources').classes('text-xl font-bold text-white mb-2')
                    ui.label('Access documentation, code, and datasets.').classes('text-slate-400 text-sm mb-6')
                    
                    with ui.column().classes('gap-3 w-full'):
                        ui.button('Read the Paper', icon='article', on_click=lambda: ui.navigate.to('https://www.biorxiv.org/content/10.64898/2026.01.17.699830v1', new_tab=True)).props('no-caps flat').classes('w-full justify-start text-left text-white bg-white/10 hover:bg-white/20')
                        ui.button('View Source Code', icon='code', on_click=lambda: ui.navigate.to('https://github.com/BGIResearch/BiOmics', new_tab=True)).props('no-caps flat').classes('w-full justify-start text-left text-white bg-white/10 hover:bg-white/20')
                        ui.button('Download Datasets', icon='download', on_click=lambda: ui.navigate.to('/datasets')).props('no-caps flat').classes('w-full justify-start text-left text-white bg-white/10 hover:bg-white/20')
                        ui.button('Analyze Your Own Data', icon='upload_file', on_click=lambda: ui.navigate.to('/platform')).props('no-caps flat').classes('w-full justify-start text-left text-white bg-white/10 hover:bg-white/20')
                
                # 第二个按钮组：产品可以做的工作
                with ui.column().classes('glass-card p-8 rounded-2xl text-left'):
                    ui.label('Examples').classes('text-xl font-bold text-white mb-2')
                    ui.label('Explore example analyses powered by BiOmics-Agent.').classes('text-slate-400 text-sm mb-6')
                    
                    with ui.column().classes('gap-3 w-full'):
                        ui.button('Cell type annotation', icon='label', on_click=lambda: ui.navigate.to('/platform?demo=Demonstrate%20cell%20annotation')).props('no-caps flat').classes('w-full justify-start text-left text-white bg-blue-500/30 hover:bg-blue-500/50 border border-blue-400/30')
                        ui.button('Differential gene analysis', icon='compare_arrows', on_click=lambda: ui.navigate.to('/platform?demo=Demonstrate%20differential%20gene%20analysis')).props('no-caps flat').classes('w-full justify-start text-left text-white bg-blue-500/30 hover:bg-blue-500/50 border border-blue-400/30')
                        ui.button('Proteome analysis', icon='hub', on_click=lambda: ui.navigate.to('/platform?demo=Demonstrate%20proteome%20analysis')).props('no-caps flat').classes('w-full justify-start text-left text-white bg-blue-500/30 hover:bg-blue-500/50 border border-blue-400/30')
                        ui.button('Gene regulatory network', icon='account_tree', on_click=lambda: ui.navigate.to('/platform?demo=Demonstrate%20gene%20regulatory%20network')).props('no-caps flat').classes('w-full justify-start text-left text-white bg-blue-500/30 hover:bg-blue-500/50 border border-blue-400/30')

    # ==================== 2. Demo Video Section ====================
    with ui.element('section').props('id="demo"').classes('bg-gradient-to-b from-slate-100 to-white py-24 px-6 w-screen').style('margin-left: calc(-50vw + 50%); box-sizing: border-box;'):
        with ui.column().classes('max-w-5xl mx-auto text-center items-center w-full'):
            ui.label('See BiOmics-Agent in Action').classes('text-4xl font-bold text-slate-900 mb-4')
            ui.label('Watch how BiOmics-Agent autonomously transforms raw omics data into actionable biological insights.').classes('text-slate-500 mb-12 max-w-2xl')
            
            with ui.element('div').classes('w-full rounded-2xl overflow-hidden shadow-2xl border border-slate-200'):
                ui.html('''
                    <video controls class="w-full" style="max-height: 600px;">
                        <source src="/static/Supplementary%20Video%201.mp4" type="video/mp4">
                        Your browser does not support the video tag.
                    </video>
                ''', sanitize=False)
            
            ui.label('Note: This demonstration video is presented at 2.5× playback speed to efficiently convey the complete workflow within a concise viewing time.').classes('text-slate-500 text-sm mt-4 leading-relaxed max-w-3xl')
            ui.label('The full-speed analysis showcases the autonomous nature of BiOmics-Agent from data input to biological discovery.').classes('text-slate-400 text-sm mt-2 italic')

    # ==================== 3. Capabilities Section ====================
    with ui.element('section').props('id="capabilities"').classes('py-24 px-6 max-w-7xl mx-auto'):
        with ui.row().classes('flex flex-col md:flex-row items-end justify-between mb-16 gap-4 w-full'):
            with ui.column().classes('max-w-2xl'):
                ui.label('Triadic Core Capabilities').classes('text-3xl font-bold text-slate-900 mb-4')
                ui.label('Our unified embedding space allows for precise retrieval, traceable reasoning, and novel hypothesis prediction.').classes('text-slate-600')
            ui.element('div').classes('h-px flex-grow bg-slate-200 mx-8 hidden md:block mb-4')
        
        capabilities = [
            ('fa-fingerprint', 'cyan', 'Retrieving', 'Precise biomedical knowledge', 'Long-chain relation retrieval from BiOmics-KG to identify and rank valuable biological findings.'),
            ('fa-project-diagram', 'emerald', 'Reasoning', 'Evidence-based discovery', 'Explicit reasoning space providing traceable causal discovery and autonomous scientific validation.'),
            ('fa-lightbulb', 'amber', 'Predicting', 'Novel unified embedding', 'Link prediction within a unified embedding space to forecast previously unknown biological relations.'),
        ]
        
        with ui.row().classes('grid md:grid-cols-3 gap-10 w-full'):
            for icon, color, title, subtitle, desc in capabilities:
                with ui.column().classes(f'capability-card-{color} group'):
                    ui.html(f'''
                        <div class="capability-icon w-16 h-16 bg-{color}-50 text-{color}-600 rounded-2xl flex items-center justify-center mb-6 text-2xl">
                            <i class="fas {icon}"></i>
                        </div>
                    ''', sanitize=False)
                    ui.label(title).classes('text-xl font-bold mb-3 text-slate-800')
                    ui.label(subtitle).classes('text-slate-600 mb-4 italic text-sm')
                    ui.label(desc).classes('text-slate-500 text-sm leading-relaxed')

    # ==================== 4. Framework Overview Section ====================
    with ui.element('section').props('id="framework"').classes('bg-slate-50 py-24 px-6 w-screen').style('margin-left: calc(-50vw + 50%); box-sizing: border-box;'):
        with ui.column().classes('max-w-7xl mx-auto w-full items-center'):
            # 标题
            ui.label('BiOmics Framework').classes('text-4xl font-bold text-slate-900 mb-4 text-center')
            ui.label('A tripartite architecture unifying Knowledge, Tools, and Agents for biological interpretation.').classes('text-slate-500 mb-12 text-center max-w-3xl')
            
            # 展示图片
            with ui.element('div').classes('w-full rounded-2xl overflow-hidden shadow-2xl border border-slate-200 mb-16'):
                ui.image('/static/Figure1.png').classes('w-full')
            
            # 三大组件介绍
            with ui.row().classes('grid md:grid-cols-3 gap-8 w-full items-stretch'):
                # BiOmics-KG
                with ui.card().classes('p-8 bg-white border border-slate-200 rounded-2xl shadow-lg hover:shadow-xl transition-shadow h-full'):
                    with ui.row().classes('items-center gap-3 mb-4'):
                        ui.icon('hub').classes('text-blue-600').style('font-size: 32px;')
                        ui.label('BiOmics-KG').classes('text-2xl font-bold text-slate-800')
                    ui.label('A foundational knowledge memory of 350 million daily-updated relations to ground inference and mitigate stochastic hallucinations.').classes('text-slate-600 mb-4 leading-relaxed')
                    with ui.column().classes('gap-2 text-sm text-slate-500'):
                        ui.label('• 6M+ Publications with daily PubMed updates')
                        ui.label('• 23+ Ontologies (GO, HPO, Cell Ontology)')
                        ui.label('• 89+ Public Databases (DrugBank, HMDB, ClinGen)')
                
                # BiOmics-BRICK
                with ui.card().classes('p-8 bg-white border border-slate-200 rounded-2xl shadow-lg hover:shadow-xl transition-shadow h-full'):
                    with ui.row().classes('items-center gap-3 mb-4'):
                        ui.icon('build').classes('text-amber-600').style('font-size: 32px;')
                        ui.label('BiOmics-BRICK').classes('text-2xl font-bold text-slate-800')
                    ui.label('A modular, pluggable toolchain to overcome bioinformatic interoperability bottlenecks with unified infrastructure.').classes('text-slate-600 mb-4 leading-relaxed')
                    with ui.column().classes('gap-2 text-sm text-slate-500'):
                        ui.label('• Query Graph & Rank Graph for retrieval')
                        ui.label('• Embedding & Reasoning for inference')
                        ui.label('• Preprocessing & Visualization modules')
                
                # BiOmics-Agent
                with ui.card().classes('p-8 bg-white border border-slate-200 rounded-2xl shadow-lg hover:shadow-xl transition-shadow h-full'):
                    with ui.row().classes('items-center gap-3 mb-4'):
                        ui.icon('smart_toy').classes('text-emerald-600').style('font-size: 32px;')
                        ui.label('BiOmics-Agent').classes('text-2xl font-bold text-slate-800')
                    ui.label('The logical orchestrator for high-order autonomous path planning and hypothesis generation from data to interpretation.').classes('text-slate-600 mb-4 leading-relaxed')
                    with ui.column().classes('gap-2 text-sm text-slate-500'):
                        ui.label('• Planning: Task decomposition & analysis design')
                        ui.label('• Coding & Execution: Automated preprocessing')
                        ui.label('• Interpretation: Evidence-backed reasoning')

    # ==================== Citation Section ====================
    with ui.element('section').props('id="citation"').classes('bg-white py-20 px-6 w-screen border-t border-slate-200').style('margin-left: calc(-50vw + 50%); box-sizing: border-box;'):
        with ui.column().classes('max-w-4xl mx-auto w-full items-center'):
            # 标题
            ui.label('Cite Our Work').classes('text-3xl font-bold text-slate-900 mb-4 text-center')
            ui.label('If you find BiOmics useful in your research, please consider citing our paper:').classes('text-slate-500 mb-8 text-center')
            
            # 论文信息卡片
            with ui.card().classes('w-full bg-slate-50 border border-slate-200 rounded-2xl p-6 mb-6'):
                ui.label('BiOmics: A Foundational Agent for Grounded and Autonomous Multi-omics Interpretation').classes('text-slate-800 font-semibold text-lg mb-2')
                ui.label('Cao Lei, Li Yuntian, Qin Hua, Shang Yanbang, Zhang Yilin, Jovanovic Bogdan, Djokic Lazar, Xia Tianyi, Hu Luni, Hou Haiyang, Ning Xingxing, Lin Li\'ang, Qiu Hao, Deng Ziqing, Li Yuxiang, Zhang Yong, Fang Shuangsang').classes('text-slate-500 text-sm mb-2')
                with ui.row().classes('gap-4 text-sm'):
                    ui.label('bioRxiv 2026').classes('text-blue-600')
                    ui.label('DOI: 10.64898/2026.01.17.699830').classes('text-slate-400')
            
            # BibTeX 引用框
            ui.label('BibTeX').classes('text-slate-800 font-semibold mb-3 self-start')
            bibtex_text = '''@article {Cao2026.01.17.699830,
    author = {Cao, Lei and Li, Yuntian and Qin, Hua and Shang, Yanbang and Zhang, Yilin and Jovanovic, Bogdan and Djokic, Lazar and Xia, Tianyi and Hu, Luni and Hou, Haiyang and Ning, Xingxing and Lin, Li'ang and Qiu, Hao and Deng, Ziqing and Li, Yuxiang and Zhang, Yong and Fang, Shuangsang},
    title = {BiOmics: A Foundational Agent for Grounded and Autonomous Multi-omics Interpretation},
    year = {2026},
    doi = {10.64898/2026.01.17.699830},
    publisher = {Cold Spring Harbor Laboratory},
    journal = {bioRxiv}
}'''
            with ui.element('div').classes('w-full bg-slate-100 rounded-xl p-4 border border-slate-200'):
                ui.code(bibtex_text, language='bibtex').classes('text-sm')
            
            # 复制按钮
            with ui.row().classes('mt-4 gap-4'):
                ui.button('Copy BibTeX', icon='content_copy', on_click=lambda: ui.clipboard.write(bibtex_text)).props('no-caps').classes('bg-blue-600 text-white hover:bg-blue-700')
                ui.button('View on bioRxiv', icon='open_in_new', on_click=lambda: ui.navigate.to('https://www.biorxiv.org/content/10.64898/2026.01.17.699830v1', new_tab=True)).props('no-caps outline').classes('text-slate-700')

    # ==================== Footer ====================
    with ui.element('footer').classes('bg-slate-50 py-12 border-t border-slate-200 w-screen').style('margin-left: calc(-50vw + 50%); box-sizing: border-box;'):
        with ui.row().classes('max-w-7xl mx-auto px-6 flex flex-col md:flex-row justify-between items-center text-slate-500 text-sm w-full'):
            with ui.row().classes('items-center gap-2 mb-4 md:mb-0'):
                ui.label('BiOmics').classes('font-bold text-slate-800')
                ui.label('Framework & Agent Platform')
            
            with ui.row().classes('gap-6'):
                ui.label('© 2026 Analysis-Interpretation Paradigm')
                ui.link('Contact Us', 'mailto:fangshuangsang@genomics.cn').classes('hover:text-blue-600')


# 创建页面
@ui.page('/')
def landing_page():
    create_landing_page()

@ui.page('/datasets')
def datasets_page():
    """数据下载页面"""
    # 数据文件目录
    DATA_DIR = '/home/liyuntian/Biomics_agent/data'
    
    with ui.column().classes('max-w-5xl mx-auto py-16 px-6 w-full'):
        # 返回按钮
        ui.button('Back to Home', icon='arrow_back', on_click=lambda: ui.navigate.to('/')).props('no-caps flat').classes('mb-8')
        
        # 标题
        ui.label('Download Datasets').classes('text-4xl font-bold text-slate-900 mb-4')
        ui.label('Download the datasets used in our demonstration analyses. Each dataset is paired with a specific analysis task.').classes('text-slate-500 mb-12')
        
        # 数据集列表
        datasets = [
            ('Cell Type Annotation', 'adata_new1.h5ad', 'Perform cell type annotation on this dataset', 'Single-cell RNA-seq data for cell type identification'),
            ('Cell Type Refinement', 'adata_new1.h5ad', 'Perform cell type refinement on this dataset', 'Same dataset for refining cell type labels'),
            ('Differential Gene Analysis', 'adata_new1.h5ad', 'Perform differential gene expression analysis on this dataset', 'Identify differentially expressed genes between conditions'),
            ('Drug Discovery', 'Neutrophil_adata_sub.h5ad', 'Predict therapeutic drugs for COVID-19 based on this omics data', 'Neutrophil subset data for drug target prediction'),
            ('Enrichment Analysis', 'adata_new1.h5ad', 'Perform gene enrichment analysis on this dataset', 'Gene set enrichment analysis'),
            ('GWAS Causal SNPs', 'filtered_mutation.csv', 'Identify causal SNPs associated with type 2 diabetes using this data', 'Filtered SNP mutation data for GWAS analysis'),
            ('GWAS Phenotype Prediction', 'filtered_mutation.csv', 'Predict associated phenotypes based on SNPs in this data', 'Same mutation data for phenotype prediction'),
            ('Trajectory Analysis', 'processed_wbc_m_group1.h5ad', 'Perform trajectory inference analysis on this dataset', 'White blood cell data for developmental trajectory'),
            ('Proteome Analysis', 'Phosphopeptides_glycopeptides_evidence_TiO2_TMT_HUMAN.h5ad', 'Perform proteome analysis on this dataset', 'Phosphopeptide and glycopeptide proteomics data'),
            ('Gene Regulatory Network', 'regulon_0619_modules.tsv', 'Perform gene regulatory network analysis on this dataset', 'Regulon modules for GRN inference'),
        ]
        
        with ui.column().classes('gap-4 w-full'):
            for task, filename, query, description in datasets:
                with ui.card().classes('w-full p-6'):
                    with ui.row().classes('items-start justify-between w-full gap-4'):
                        with ui.column().classes('flex-grow gap-1'):
                            ui.label(task).classes('text-lg font-bold text-slate-800')
                            ui.label(description).classes('text-slate-500 text-sm')
                            with ui.row().classes('items-center gap-2 mt-2'):
                                ui.icon('insert_drive_file').classes('text-blue-500')
                                ui.label(filename).classes('text-blue-600 font-mono text-sm')
                        # 为每个文件创建下载按钮
                        file_path = f'{DATA_DIR}/{filename}'
                        ui.button('Download', icon='download', on_click=lambda fp=file_path, fn=filename: ui.download(fp, fn)).props('no-caps outline').classes('flex-shrink-0')
        


if __name__ in {"__main__", "__mp_main__"}:
    ui.run(title='BiOmics | Seamless Omics Analysis & Interpretation', favicon='/home/liyuntian/Biomics_agent/BiomicsLOGO.svg', host='0.0.0.0', port=8080,  uvicorn_logging_level='debug', proxy_headers=True, forwarded_allow_ips="*")

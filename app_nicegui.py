from ast import Dict
from time import sleep
import asyncio
import threading
import os
import re
from datetime import datetime
from dotenv import load_dotenv
from pathlib import Path
from nicegui import app, ui
from langchain_core.messages import HumanMessage
from utils.build_config_id import build_config_id
from graph.state import BrickState
from graph.builder import build_graph_with_interaction
from utils.save_dir_name import get_save_dir
from utils.sandbox_manager import SandboxManager
from utils.plan_extracter import plan_exetract
from utils.create_notebook import create_notebook

# 演示按钮配置：按钮文本 -> (文件路径, 问题)
DEMO_BUTTON_CONFIG = {
    'Demonstrate cell annotation': ('/home/liyuntian/Biomics_agent/data/adata_new1.h5ad', 'Perform cell type annotation on this dataset'),
    'Demonstrate cell refinement': ('/home/liyuntian/Biomics_agent/data/adata_new1.h5ad', 'Perform cell type refinement on this dataset'),
    'Demonstrate differential gene analysis': ('/home/liyuntian/Biomics_agent/data/adata_new1.h5ad', 'Perform differential gene expression analysis on this dataset'),
    'Demonstrate drug discovery': ('/home/liyuntian/Biomics_agent/data/Neutrophil_adata_sub.h5ad', 'Predict therapeutic drugs for COVID-19 based on this omics data'),
    'Demonstrate enrichment analysis': ('/home/liyuntian/Biomics_agent/data/adata_new1.h5ad', 'Perform gene enrichment analysis on this dataset'),
    'Demonstrate GWAS causal SNPs analysis': ('/home/liyuntian/Biomics_agent/data/filtered_mutation.csv', 'Identify causal SNPs associated with type 2 diabetes using this data'),
    'Demonstrate GWAS phenotype prediction': ('/home/liyuntian/Biomics_agent/data/filtered_mutation.csv', 'Predict associated phenotypes based on SNPs in this data'),
    'Demonstrate trajectory analysis': ('/home/liyuntian/Biomics_agent/data/processed_wbc_m_group1.h5ad', 'Perform trajectory inference analysis on this dataset'),
}

def build_tree_nodes(data, prefix=''):
    """将嵌套字典转换为 ui.tree 所需的节点格式"""
    nodes = []
    if not isinstance(data, dict):
        return nodes
    for key, value in data.items():
        node_id = f"{prefix}_{key}" if prefix else key
        if isinstance(value, dict):
            children = build_tree_nodes(value, node_id)
            nodes.append({
                'id': node_id,
                'label': key,
                'children': children if children else None
            })
        elif isinstance(value, list):
            preview = str(value[:3]) + '...' if len(value) > 3 else str(value)
            nodes.append({
                'id': node_id,
                'label': f"{key}: {preview}"
            })
        else:
            display_val = str(value)[:50] + '...' if len(str(value)) > 50 else str(value)
            nodes.append({
                'id': node_id,
                'label': f"{key}: {display_val}"
            })
    return nodes

# 加载环境变量
config_file = Path(__file__).parent / 'graph' / 'brick_test_config.env'
load_dotenv(dotenv_path=str(config_file))
PROJECT_ROOT = os.getenv('PROJECT_ROOT', os.path.abspath(os.path.dirname(__file__)))

# 添加静态文件目录，使 logo 等资源可访问
app.add_static_files('/static', PROJECT_ROOT)

ui.query('body').style('margin: 0; padding: 0; overflow: hidden;')

# 添加 iMessage 风格的动画和自定义样式
ui.add_head_html('''
<style>
@keyframes slideInFade {
    from {
        opacity: 0;
        transform: translateY(20px) scale(0.95);
    }
    to {
        opacity: 1;
        transform: translateY(0) scale(1);
    }
}
/* 覆盖默认的绿色消息背景，改为灰色 */
.q-message-text--received {
    background-color: #e0e0e0 !important;
    color: #333 !important;
}
/* 用户发送的消息改为蓝色 */
.q-message-text--sent {
    background-color: #1976d2 !important;
    color: white !important;
}
.q-message-text--sent > div {
    color: white !important;
}
/* 隐藏消息气泡的三角形箭头 */
.q-message-text:before,
.q-message-text:after {
    display: none !important;
}
/* 缩小消息框内 Markdown 标题的字体大小 */
.q-message-text h1 {
    font-size: 1.25em !important;
    margin: 0.3em 0 !important;
}
.q-message-text h2 {
    font-size: 1.1em !important;
    margin: 0.25em 0 !important;
}
.q-message-text h3 {
    font-size: 1em !important;
    margin: 0.2em 0 !important;
}
.q-message-text h4, .q-message-text h5, .q-message-text h6 {
    font-size: 0.95em !important;
    margin: 0.15em 0 !important;
}
/* 限制 code 组件宽度，防止擑开容器 */
.nicegui-code, .nicegui-code pre, .q-card pre {
    max-width: 100% !important;
    overflow-x: auto !important;
    white-space: pre-wrap !important;
    word-break: break-word !important;
}
.q-card {
    max-width: 100% !important;
    overflow: hidden !important;
}
/* 浮动按钮栏动画 */
@keyframes fadeInFromLeft {
    from {
        opacity: 0;
        transform: translateX(-20px);
    }
    to {
        opacity: 1;
        transform: translateX(0);
    }
}
.floating-btn {
    opacity: 0;
    animation: fadeInFromLeft 0.5s ease-out forwards;
}
.floating-btn-1 { animation-delay: 3s; }
.floating-btn-2 { animation-delay: 3.3s; }
.floating-btn-3 { animation-delay: 3.6s; }
.floating-btn-4 { animation-delay: 3.9s; }
.floating-btn-5 { animation-delay: 4.2s; }
.floating-btn-6 { animation-delay: 4.5s; }
.floating-btn-7 { animation-delay: 4.8s; }
.floating-btn-8 { animation-delay: 5.1s; }
.floating-btn-9 { animation-delay: 5.4s; }

</style>
''')

# === 首次打开弹窗 ===
welcome_dialog = ui.dialog()

with welcome_dialog, ui.card().style('min-width: 400px; padding: 24px;'):
    with ui.row().style('width: 100%; justify-content: space-between; align-items: center; margin-bottom: 12px;'):
        ui.label('Welcome to BiOmics Agent - Usage Tips').style('font-size: 18px; font-weight: bold;')
        ui.button(icon='close', on_click=welcome_dialog.close).props('flat round dense color=red')
    ui.label('• Sessions are temporary and not persisted. Please do not close the session midway and download result files promptly after completion.').style('font-size: 14px; color: #555; margin-bottom: 8px;')
    ui.label('• Our computing resources are limited. Please avoid processing computationally complex tasks.').style('font-size: 14px; color: #555;')
    with ui.row().style('width: 100%; justify-content: flex-end;'):
        ui.label("   ")

# 每个客户端加载页面时弹出一次
ui.timer(0.1, lambda: welcome_dialog.open(), once=True)

# 主容器:占满整个视口
with ui.column().style('width: 100vw; height: 100vh; margin: 0; padding: 0;'):
    
    # === 加载提示悬浮框（小巧透明圆角，初始隐藏） ===
    loading_banner = ui.element('div').style(
        'display: none; position: fixed; bottom: 120px; left: 50%; '
        'transform: translateX(-50%); '
        'background: rgba(25, 118, 210, 0.9); color: white; '
        'padding: 16px 24px; border-radius: 24px; '
        'box-shadow: 0 4px 12px rgba(0,0,0,0.15); '
        'backdrop-filter: blur(10px); z-index: 9999;'
    )
    with loading_banner:
        with ui.row().style('align-items: center; gap: 12px;'):
            ui.spinner(size='sm', color='white')
            ui.label('Agent running, please wait...').style('color: white; font-weight: 500; font-size: 14px;')
    
    # === 顶栏 ===
    with ui.row().style(
        'height: 10%; min-height: 60px; width: 100%; '
        'align-items: center; padding: 0 24px; '
        'border-bottom: 1px solid #ddd; background: #f5f5f5; '
        'justify-content: space-between;'
    ):
        # 左侧：Logo 和标题
        with ui.row().style('align-items: center;'):
            ui.html('''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 412.27 159.57" style="height: 40px; width: auto;">
              <defs><style>.b{fill:#231815;font-family:Arial-BoldMT,Arial;font-size:86.59px;font-weight:700}.c{fill:#dcdddd}.c,.d,.e{stroke:#231815;stroke-miterlimit:10}.c,.e{stroke-width:.75px}.d{fill:#c9caca;stroke-linecap:round;stroke-width:.25px}.e{fill:none}</style></defs>
              <text class="b" transform="translate(0 74.28)"><tspan x="0" y="0">Bi</tspan></text>
              <g><path class="e" d="M210.19,72.43c0-29.58-23.11-53.68-51.96-54.49,14.05.82,25.19,12.7,25.19,27.22s-11.97,27.27-26.77,27.27-26.75,12.22-26.75,27.27,11.14,26.41,25.19,27.22c.52.05,1.04.05,1.56.05s1.06,0,1.58-.05c28.85-.82,51.96-24.91,51.96-54.49Z"/><path class="c" d="M183.42,45.16c0-14.53-11.14-26.41-25.19-27.22-.52-.02-1.06-.05-1.58-.05s-1.04.02-1.56.05c-28.83.82-51.96,24.91-51.96,54.49s23.14,53.68,51.96,54.49c-14.05-.82-25.19-12.7-25.19-27.22s11.97-27.27,26.75-27.27,26.77-12.22,26.77-27.27Z"/><ellipse class="d" cx="125.18" cy="81.51" rx="3.57" ry="3.64"/><ellipse class="d" cx="145.94" cy="63.41" rx="3.57" ry="3.64"/><ellipse class="d" cx="166.39" cy="79.09" rx="3.57" ry="3.64"/><ellipse class="d" cx="188.77" cy="63.41" rx="3.57" ry="3.64"/><ellipse class="d" cx="145.94" cy="41.7" rx="3.57" ry="3.64"/><ellipse class="d" cx="166.39" cy="101.88" rx="3.57" ry="3.64"/><line class="d" x1="145.94" y1="63.41" x2="166.39" y2="79.09"/><line class="d" x1="166.39" y1="79.09" x2="166.39" y2="101.88"/><line class="d" x1="166.39" y1="79.09" x2="188.77" y2="63.41"/><line class="d" x1="125.18" y1="81.51" x2="145.94" y2="63.41"/><line class="d" x1="145.94" y1="63.41" x2="145.94" y2="41.7"/></g>
              <text class="b" transform="translate(214.91 126.97)"><tspan x="0" y="0">mics</tspan></text>
            </svg>''', sanitize=False)
            ui.label('Agent').style('font-size: 24px; font-weight: bold; margin-left: 8px;')
        
        # 右侧：按钮组
        with ui.row().style('gap: 20px; margin-right: 40px;'):
            with ui.button(icon='', on_click=lambda: ui.navigate.to('https://github.com/BGIResearch/BiOmics', new_tab=True)).props('outline').style('border-color: #333; background-color: #fff !important; padding-left: 1px;').classes('text-black'):
                ui.html('<svg height="20" width="20" viewBox="0 0 16 16" style="fill: #333;"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>', sanitize=False).style('margin: -4px;')
                ui.label('Github').style('margin-left: 15px; color: #333; text-transform: none;')
            ui.button('Contact Us', icon='mail', on_click=lambda: ui.navigate.to('mailto:fangshuangsang@genomics.cn')).props('outline no-caps').style('border-color: #333; background-color: #fff !important;').classes('text-black')

    # === 中间区域===
    with ui.element('div').style(
        'height: 75%; width: 100%; margin: 0 auto; display: flex; flex-direction: row;'
    ):
        
        # 左侧:对话区
        with ui.element('div').style(
            'width: 50%; height: 100%; padding: 16px; '
            'border-right: 1px solid #ddd; display: flex; flex-direction: column;'
        ):
            ui.label('💬 BiOmics Chat').style('font-size: 18px; font-weight: 600; margin-bottom: 8px;')
            biomics_chat = ui.scroll_area().style('width: 100%; flex: 1;')


        
        # 右侧:代码区
        with ui.element('div').style(
            'width: 50%; height: 100%; padding: 16px; display: flex; flex-direction: column;'
        ):
            ui.label('✨ BiOmics Co-pilot').style('font-size: 18px; font-weight: 600; margin-bottom: 8px;')
            biomics_co_pilot = ui.scroll_area().style('width: 100%; flex: 1; ')


    # === 浮动按钮栏（透明，位于对话栏上方左侧） ===
    floating_btn_bar = ui.element('div').style(
        'position: fixed; bottom: 100px; left: 100px; '
        'background: transparent; z-index: 100; '
        'display: flex; gap: 12px;'
    )
    with floating_btn_bar:
        ui.button('Demonstrate cell annotation', on_click=lambda: handle_demo_button_click('Demonstrate cell annotation')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-1')
        ui.button('Demonstrate cell refinement', on_click=lambda: handle_demo_button_click('Demonstrate cell refinement')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-2')
        ui.button('Demonstrate differential gene analysis', on_click=lambda: handle_demo_button_click('Demonstrate differential gene analysis')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-3')
        ui.button('Demonstrate drug discovery', on_click=lambda: handle_demo_button_click('Demonstrate drug discovery')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-4')
        ui.button('Demonstrate enrichment analysis', on_click=lambda: handle_demo_button_click('Demonstrate enrichment analysis')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-5')
        ui.button('Demonstrate GWAS causal SNPs analysis', on_click=lambda: handle_demo_button_click('Demonstrate GWAS causal SNPs analysis')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-6')
        ui.button('Demonstrate GWAS phenotype prediction', on_click=lambda: handle_demo_button_click('Demonstrate GWAS phenotype prediction')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-7')
        ui.button('Demonstrate trajectory analysis', on_click=lambda: handle_demo_button_click('Demonstrate trajectory analysis')).props('outline no-caps').style('border-radius: 8px; font-size: 10px; white-space: nowrap;').classes('floating-btn floating-btn-8')

    # === 底部对话栏 ===
    with ui.row().style(
        'height: 10%; min-height: 60px; width: 100%; '
        'align-items: center; padding: 0 24px; gap: 12px;'
    ):
        # 左侧:上传图标 + 文件名
        with ui.row().style('align-items: center; gap: 4px;'):
            # 隐藏的上传控件
            file_upload = ui.upload(
                auto_upload=True,
            ).style('display: none;')  # 完全隐藏
            
            # 显示的图标按钮，点击时触发上传
            upload_button = ui.button(icon='file_upload', on_click=lambda: file_upload.run_method('pickFiles')).props('flat round dense')
            
            # 文件名标签
            upload_name_label = ui.label('').style(
                'font-size: 12px; color: #666; max-width: 200px; '
                'overflow: hidden; text-overflow: ellipsis; white-space: nowrap;'
            )
        
        # 中间:输入框
        user_input = ui.input(placeholder='Enter bioinformatics analysis task...').style('flex: 1;')
        
        # 右侧:重置按钮
        reset_button = ui.button('Reset', icon='restart_alt').props('outlined')


UPLOAD_DIR = os.path.join(PROJECT_ROOT, 'data', 'uploaded_file')

async def handle_demo_button_click(button_name: str) -> None:
    """处理演示按钮点击：隐藏按钮栏、设置文件路径、启动任务"""
    if button_name not in DEMO_BUTTON_CONFIG:
        ui.notify(f'Unknown demo: {button_name}', type='warning')
        return
    
    file_path, question = DEMO_BUTTON_CONFIG[button_name]
    
    # 隐藏浮动按钮栏
    floating_btn_bar.style('display: none;')
    
    # 设置文件路径到 storage
    app.storage.client['uploaded_file_path'] = file_path
    upload_name_label.text = os.path.basename(file_path)
    upload_button.props('icon=check_circle color=positive')
    
    # 设置问题并启动任务（handle_user_input 会显示 chat_message）
    user_input.value = question
    await handle_user_input()

def set_graph_running(is_running: bool) -> None:
    """统一控制图是否在运行，以及相关控件的可用状态"""
    app.storage.client['graph_running'] = is_running
    
    if is_running:
        # 禁用输入与上传
        user_input.props('disable')
        upload_button.props('disable')
        
        # 显示顶部加载横幅
        loading_banner.style('display: flex;')
    else:
        # 恢复输入与上传
        user_input.props(remove='disable')
        upload_button.props(remove='disable')
        
        # 隐藏顶部加载横幅
        loading_banner.style('display: none;')
def save_uploaded_file(e) -> str:
    """保存上传文件到固定目录，并返回保存路径"""
    os.makedirs(UPLOAD_DIR, exist_ok=True)
    # NiceGUI 的 UploadEventArguments: 文件名在 e.name，内容在 e.content
    # 如果 e.name 不存在，尝试从 content 获取
    try:
        file_name = e.name
    except AttributeError:
        # 如果没有 name 属性，尝试使用默认名
        file_name = getattr(e, 'filename', 'uploaded_file.h5ad')
    
    save_path = os.path.join(UPLOAD_DIR, file_name)
    with open(save_path, 'wb') as f:
        f.write(e.content.read())
    return save_path
async def handle_file_upload(e) -> None:
    """上传事件回调：在时间戳文件夹中保存文件并更新图标旁文字"""
    
    # 如果图正在运行，禁止上传
    if app.storage.client.get('graph_running', False):
        ui.notify('Current task is running, cannot upload file now', type='warning')
        return
    
    os.makedirs(UPLOAD_DIR, exist_ok=True)
    
    # NiceGUI 的 SmallFileUpload 对象，read() 是异步方法
    file_name = e.file.name
    
    # 检查文件大小（500MB限制）
    file_content = await e.file.read()  # 必须 await
    file_size = len(file_content)
    max_size = 500 * 1024 * 1024
    if file_size > max_size:
        ui.notify(f'File size exceeds 500MB limit (current: {file_size / (1024*1024):.1f}MB)', type='negative')
        return
    
    # 创建以时间戳命名的子文件夹
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    upload_subfolder = os.path.join(UPLOAD_DIR, f'upload_{timestamp}')
    os.makedirs(upload_subfolder, exist_ok=True)
    
    # 保存文件到时间戳文件夹
    save_path = os.path.join(upload_subfolder, file_name)
    with open(save_path, 'wb') as f:
        f.write(file_content)
    
    # 更新右侧的小文字为文件名
    upload_name_label.text = os.path.basename(save_path)
    # 把路径存到 client storage，后续 start_graph 时可以使用
    app.storage.client['uploaded_file_path'] = save_path
    # 图标变成绿色打勾，提示上传成功
    upload_button.props('icon=check_circle color=positive')
    ui.notify(f'File uploaded: {os.path.basename(save_path)}', type='positive')
def agent_update_chat(event) -> None:
    agent_name = event.get('agent')   
    agent_thought = ""
    agent_output = ""
    
    # 添加动画样式：模仿 iMessage
    message_animation = '''
        opacity: 0;
        animation: slideInFade 0.4s ease-out forwards;
    '''
    
    if agent_name=="supervisor":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=[agent_thought, agent_output], name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)  # 滚动到底部
    elif agent_name=="env_checker":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
            with ui.chat_message(name=agent_name):
                ui.markdown(agent_output)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="data_analyzer":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
            
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="analyze_planner":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="planner":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="plan_executor":
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_output, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="coder":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="code_runner":
        pf = event.get('process_flag')
        if pf==1:
            ui.chat_message(text="Code execution completed.Ready to process the next step.", name=agent_name).style(message_animation)
        else:
            ui.chat_message(text="Code execution failed.Ready to debug.", name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="code_debugger":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="responder":
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_output, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="notebook_searcher":
        pass
    elif agent_name=="general_responder":
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_output, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    else:
        print("未获取到", agent_name)
def agent_update_copilot(event) -> None:
    agent_name = event.get('agent')   
    if agent_name=="supervisor":
        pass
    elif agent_name=="env_checker":
        di = event.get('data_info')
        if event.get("status")=="AWAITING_CONFIRMATION":
            with biomics_co_pilot:
                with ui.card().style('width: 100%;'):
                    ui.label('Env Checker waiting for confirmation:')
                    ui.label('Please input your comfirmation in the box below')
        elif event.get("status")=="VALIDATED":
            with biomics_co_pilot:
                with ui.card().style('width: 100%;'):
                    ui.label('Env Checker checked the data.')
                    if di:
                        if isinstance(di, dict):
                            tree_nodes = build_tree_nodes(di.get('data_info', di))
                            ui.tree(tree_nodes, label_key='label', children_key='children').props('default-expand-all')
                        else:
                            ui.markdown(di)
        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="data_analyzer":
        data_report = event.get('output')
        with biomics_co_pilot:
            with ui.card().style('width: 100%;'):
                ui.label('Data Analyzer generated a report:')
                ui.markdown(data_report)

        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="analyze_planner":
        a_plan = event.get('output')
        with biomics_co_pilot:
            with ui.card().style('width: 100%;'):
                ui.label('Analyze Planner generated a plan:')
                ui.markdown(a_plan)

        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="planner":
        plan = event.get('output')
        plan = plan_exetract(plan)
        with biomics_co_pilot:
            with ui.card().style('width: 100%;'):
                ui.label('Plan Check List').style('font-size: 16px; font-weight: bold; margin-bottom: 12px;')
                for idx, step in enumerate(plan, 1):
                    with ui.row().style('width: 100%; align-items: center; padding: 8px 0; border-bottom: 1px solid #eee;'):
                        ui.label(str(idx)).style('width: 24px; height: 24px; background: #1976d2; color: white; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 12px; font-weight: bold;')
                        ui.label(step).style('flex: 1; margin: 0 12px;')
                        ui.icon('radio_button_unchecked').style('color: #bbb; font-size: 20px;')

        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="plan_executor":
        step_num = event.get("step_num")
        current_step = event.get("current_step")
        current_plan = event.get("current_plan")
        is_end = current_step >= step_num
        with biomics_co_pilot:
            if not is_end:
                with ui.card().style('width: 100%;'):
                    ui.label(f'Executing Step {current_step}').style('font-size: 16px; font-weight: bold; margin-bottom: 12px;')
                    for idx, step in enumerate(current_plan, 1):
                        step_name = step.get('type', str(step)) if isinstance(step, dict) else str(step)
                        with ui.row().style('width: 100%; align-items: center; padding: 8px 0; border-bottom: 1px solid #eee;'):
                            ui.label(str(idx)).style('width: 24px; height: 24px; background: #1976d2; color: white; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 12px; font-weight: bold;')
                            if idx < current_step:
                                ui.label(step_name).style('flex: 1; margin: 0 12px; color: #333;')
                            elif idx == current_step:
                                ui.label(step_name).style('flex: 1; margin: 0 12px; color: #1976d2; font-weight: bold;')
                            else:
                                ui.label(step_name).style('flex: 1; margin: 0 12px; color: #bbb;')
                            if idx < current_step:
                                ui.icon('check_circle').style('color: #4caf50; font-size: 20px;')
                            else:
                                ui.icon('radio_button_unchecked').style('color: #bbb; font-size: 20px;')
            else:
                with ui.card().style('width: 100%;'):
                    ui.label('Plan Execution Completed').style('font-size: 16px; font-weight: bold; margin-bottom: 12px;')
                    for idx, step in enumerate(current_plan, 1):
                        step_name = step.get('type', str(step)) if isinstance(step, dict) else str(step)
                        with ui.row().style('width: 100%; align-items: center; padding: 8px 0; border-bottom: 1px solid #eee;'):
                            ui.label(str(idx)).style('width: 24px; height: 24px; background: #1976d2; color: white; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 12px; font-weight: bold;')
                            ui.label(step_name).style('flex: 1; margin: 0 12px; color: #333;')
                            ui.icon('check_circle').style('color: #4caf50; font-size: 20px;')
        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="coder":
        code = event.get('output')
        with biomics_co_pilot:
            with ui.card().style('width: 100%;'):
                ui.label('Coder generated code:')
                ui.code(code)

        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="code_runner":
        res = event.get('complete_output')
        save_dir = event.get('save_dir')
        if res:
            with biomics_co_pilot:
                with ui.card().style('width: 100%; max-width: 100%; overflow: hidden;'):
                    ui.label('Code Runner Output:').style('font-weight: bold;')
                    # stdout
                    stdout = res.get('stdout', [])
                    if stdout:
                        stdout_text = ''.join(stdout)
                        ui.code(stdout_text, language='text').style('background: #f5f5f5; width: 100%; max-width: 100%; overflow-x: auto; white-space: pre-wrap; word-break: break-all;')
                    # stderr
                    stderr = res.get('stderr', [])
                    if stderr:
                        stderr_text = ''.join(stderr)
                        ui.code(stderr_text, language='text').style('background: #fff3cd; color: #856404; width: 100%; max-width: 100%; overflow-x: auto; white-space: pre-wrap; word-break: break-all;')
                    # result (expression return value)
                    result_val = res.get('result')
                    if result_val:
                        ui.label('Out:').style('color: #d63384; font-weight: bold;')
                        ui.code(str(result_val), language='python').style('width: 100%; max-width: 100%; overflow-x: auto; white-space: pre-wrap; word-break: break-all;')
                    # images
                    images = res.get('images', [])
                    import base64
                    for idx, img in enumerate(images):
                        img_type = img.get('type', 'png')
                        img_data = img.get('data', '')
                        if img_data:
                            ui.image(f'data:image/{img_type};base64,{img_data}').style('width: 400px; height: auto;')
                            # 保存图片到 save_dir
                            if save_dir:
                                os.makedirs(save_dir, exist_ok=True)
                                img_path = os.path.join(save_dir, f'output_{idx}.{img_type}')
                                with open(img_path, 'wb') as f:
                                    f.write(base64.b64decode(img_data))
                                print(f'[INFO] 图片已保存: {img_path}')
                    # error
                    error = res.get('error')
                    if error:
                        ui.label('Error:').style('color: #dc3545; font-weight: bold;')
                        # 清除 ANSI 转义码
                        clean_error = re.sub(r'\x1b\[[0-9;]*[a-zA-Z]|\[\d+(?:;\d+)*m', '', error)
                        ui.code(clean_error, language='text').style('background: #f8d7da; color: #721c24; width: 100%; max-width: 100%; overflow-x: auto; white-space: pre-wrap; word-break: break-all;')
        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="code_debugger":
        code = event.get('output')
        with biomics_co_pilot:
            with ui.card().style('width: 100%;'):
                ui.label('Code Debugger generated code:')
                ui.code(code)

        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="responder":
        sid = event.get('sandbox_id')
        sd = event.get('save_dir')
        notebook_cells = event.get('notebook_cells')
        rt = event.get('relation_frame')
        print("rt", rt)
        if sid:
            sandbox_manager = SandboxManager()
            sandbox_manager.close_sandbox(sid)
            print(f"[INFO] 沙箱已关闭: {sid}")
        
        create_notebook(notebook_cells, os.path.join(sd, 'analysis.ipynb'))
        print(f"[INFO] Notebook 已保存: {os.path.join(sd, 'analysis.ipynb')}")
        

        if rt:
            import pandas as pd
            # 从字典转回 DataFrame
            df = pd.DataFrame(rt) if isinstance(rt, list) else rt
            with biomics_co_pilot:
                with ui.card().style('width: 100%;'):
                    ui.label('Relation Frame:').style('font-weight: bold;')
                    ui.aggrid.from_pandas(df).classes('w-full')
        # 压缩保存目录并提供下载按钮
        if sd and os.path.isdir(sd):
            import shutil
            zip_path = shutil.make_archive(sd, 'zip', sd)
            print(f"[INFO] 已压缩: {zip_path}")
            with biomics_co_pilot:
                with ui.card().style('width: 100%;'):
                    ui.label('Task completed, results saved').style('font-weight: bold;')
                    ui.button('Download Results', icon='download', on_click=lambda: ui.download(zip_path)).props('color=primary')
        biomics_co_pilot.scroll_to(percent=1.0)
    elif agent_name=="notebook_searcher":
        pass
    elif agent_name=="general_responder":

        biomics_co_pilot.scroll_to(percent=1.0)
    else:
        print("未获取到", agent_name)
def start_graph(question: str, file_path: str, config: dict, save_dir: str = None):
    state_data = {
        "question": question,
        "messages": [HumanMessage(content=question)],
    }
    if file_path:
        state_data["data_path"] = file_path
    if save_dir:
        state_data["save_dir"] = save_dir
    initial_state = BrickState(**state_data)
    graph = build_graph_with_interaction(interrupt_after=["env_checker","data_analyzer"],interrupt_before=[])
    initial_state_dict = initial_state.model_dump()
    events = graph.stream(initial_state_dict, config=config, stream_mode="values")
    return events, graph
def process_events(graph, events, config, cancel_event: threading.Event = None):
    """循环消费 events，直到遇到中断或 FINISHED
    返回值: (waiting_kind, graph, sandbox_id) 元组
    - waiting_kind: None / 'update_data_info' / 'update_data_repo'
    - graph: 更新后的 graph 实例
    - sandbox_id: 当前会话的沙箱ID
    """
    sandbox_id = None
    for event in events:
        # 检查是否被取消
        if cancel_event and cancel_event.is_set():
            print("[DEBUG] process_events 被取消")
            return (None, graph, sandbox_id)
        
        status = event.get("status")
        # 提取 sandbox_id 并立即保存到 storage
        if event.get('sandbox_id'):
            sandbox_id = event.get('sandbox_id')
            app.storage.client['sandbox_id'] = sandbox_id
            print(f"[DEBUG] sandbox_id 已保存: {sandbox_id}")
        agent_update_chat(event)  # 统一更新界面
        agent_update_copilot(event)
        if status == "AWAITING_CONFIRMATION":
            # env_checker 需要 update_data_info：在 UI 上提示，并记录"等待类型"
            print("[DEBUG] process_events 返回 waiting_kind = 'update_data_info'")
            return ('update_data_info', graph, sandbox_id)
        
        elif status == "Revise":
            # data_analyzer 需要 update_data_repo
            print("[DEBUG] process_events 返回 waiting_kind = 'update_data_repo'")
            return ('update_data_repo', graph, sandbox_id)
        
        elif status == "VALIDATED":
            # 根据项目记忆：VALIDATED 也需要手动继续 stream
            events = graph.stream(None, config=config, stream_mode="values")
            first_event = next(events, None)
            return process_events(graph, events, config, cancel_event)

        elif status == "NOT_FINISHED":
            continue

        elif status == "FINISHED":
            user_question = event.get("question")
            if event.get("agent")=='general_responder':
                with biomics_chat:
                    with ui.card().style('width: 100%;'):
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.label(f"Your question '{user_question}' is outside our capability scope. Try tasks like cell annotation or gene enrichment?")
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.button('Ask Again', on_click=lambda: reset_button.run_method('click')).props('flat round dense')
            elif event.get('agent')=='responder':
                with biomics_chat:
                    with ui.card():
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.label(f"Your task '{user_question}' is completed. This session will end.")
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.button('Ask Again', on_click=lambda: reset_button.run_method('click')).props('flat round dense')
            return (None, graph, sandbox_id)
    
    # 如果 events 消费完没有任何特殊状态，返回 None
    return (None, graph, sandbox_id)
async def handle_user_input():
    """统一处理用户输入：
    - 没有等待状态时，作为新问题启动一条图
    - 等待状态为 AWAITING_CONFIRMATION / Revise 时，作为反馈或修改继续图
    """

    text = (user_input.value or '').strip()
    if not text:
        ui.notify('Please enter content', type='warning')
        return
    with biomics_chat:
        ui.chat_message(text=text, name='user', sent=True).style('margin-left: auto; max-width: 80%;')
    biomics_chat.scroll_to(percent=1.0)  # 用户消息也滚动到底部
    # 清空输入框（防止重复发送）
    user_input.value = ''

    # 如果有后台任务还在跑，先等它完成（确保 waiting_kind 已经被正确设置）
    current_task = app.storage.client.get('background_task')
    if current_task and not current_task.done():
        try:
            await current_task
        except:
            pass  # 忽略后台任务异常
    
    # 读取当前会话状态：是否在等待某类中断输入
    waiting_kind = app.storage.client.get('waiting_kind')
    is_running = app.storage.client.get('graph_running', False)
    print(f"[DEBUG] handle_user_input 读取状态: waiting_kind={waiting_kind}, is_running={is_running}")

    # ============= 情况一：当前没有等待中断输入，视为"新问题" =============
    if not waiting_kind:
        if is_running:
            ui.notify('A task is currently running, please wait before asking', type='warning')
            return
        # 为本次会话生成 config（使用独立 thread_id）
        thread_id = build_config_id()
        config = {"configurable": {"thread_id": thread_id}, "recursion_limit": 200}
    
        try:
            # 隐藏浮动按钮栏
            floating_btn_bar.style('display: none;')
                
            set_graph_running(True)

            # 启动一次图（如果有上传文件，则使用上传的文件路径）
            uploaded_file = app.storage.client.get('uploaded_file_path', '')
            events, graph = start_graph(
                question=text,
                file_path=uploaded_file,
                config=config,
                save_dir=thread_id,
            )

            # 存储 graph 和 config，方便中断后继续使用
            app.storage.client['graph'] = graph
            app.storage.client['config'] = config

            # 创建取消事件，用于中断后台任务
            cancel_event = threading.Event()
            app.storage.client['cancel_event'] = cancel_event

            # 交给统一的事件处理逻辑（后台线程执行，避免阻塞UI）
            task = asyncio.create_task(asyncio.to_thread(process_events, graph, events, config, cancel_event))
            app.storage.client['background_task'] = task
            result = await task
            
            # 后台线程返回的状态，由主线程设置到 app.storage.client
            if result:
                waiting_kind, graph, sandbox_id = result
                app.storage.client['graph'] = graph
                app.storage.client['config'] = config
                app.storage.client['sandbox_id'] = sandbox_id
                if waiting_kind:
                    app.storage.client['waiting_kind'] = waiting_kind
                    set_graph_running(False)
                    print(f"[DEBUG] 主线程设置 waiting_kind = {waiting_kind}")
                else:
                    set_graph_running(False)

        except Exception as e:
            error_msg = f"❌ Failed to start task: {e}"
            with biomics_chat:
                ui.chat_message(text=error_msg, name='System')
            ui.notify(error_msg, type='negative')
        finally:
            # 是否置回 False，要看 process_events 是否进入等待状态
            # 这里仅在没有等待标记时重置
            if not app.storage.client.get('waiting_kind'):
                set_graph_running(False)

        return

    # ============= 情况二：正在等待中断节点输入（AWAITING_CONFIRMATION / Revise） =============
    # 等待输入时，不允许新的问题，text 被视为对方案的“确认/修改”反馈
    try:
        graph = app.storage.client.get('graph')
        config = app.storage.client.get('config')
        if graph is None or config is None:
            ui.notify('Internal state lost, please restart session', type='negative')
            return


        if waiting_kind == 'update_data_info':
            print("[DEBUG] 进入 update_data_info 分支")
            graph.update_state(config=config, values={'update_data_info': text})
            app.storage.client['waiting_kind'] = None
            set_graph_running(True)
            events = graph.stream(None, config=config, stream_mode='values')
            first_event = next(events, None)
            cancel_event = threading.Event()
            app.storage.client['cancel_event'] = cancel_event
            task = asyncio.create_task(asyncio.to_thread(process_events, graph, events, config, cancel_event))
            app.storage.client['background_task'] = task
            result = await task
            
            # 后台线程返回的状态，由主线程设置
            if result:
                waiting_kind, graph, sandbox_id = result
                app.storage.client['graph'] = graph
                app.storage.client['config'] = config
                app.storage.client['sandbox_id'] = sandbox_id
                if waiting_kind:
                    app.storage.client['waiting_kind'] = waiting_kind
                    set_graph_running(False)
                else:
                    set_graph_running(False)

        elif waiting_kind == 'update_data_repo':
            print("[DEBUG] 进入 update_data_repo 分支")
            graph.update_state(config=config, values={'update_data_repo': text})
            app.storage.client['waiting_kind'] = None
            set_graph_running(True)
            events = graph.stream(None, config=config, stream_mode='values')
            first_event = next(events, None)
            cancel_event = threading.Event()
            app.storage.client['cancel_event'] = cancel_event
            task = asyncio.create_task(asyncio.to_thread(process_events, graph, events, config, cancel_event))
            app.storage.client['background_task'] = task
            result = await task
            
            # 后台线程返回的状态，由主线程设置
            if result:
                waiting_kind, graph, sandbox_id = result
                app.storage.client['graph'] = graph
                app.storage.client['config'] = config
                app.storage.client['sandbox_id'] = sandbox_id
                if waiting_kind:
                    app.storage.client['waiting_kind'] = waiting_kind
                    set_graph_running(False)
                else:
                    set_graph_running(False)
        else:
            ui.notify(f'Unknown waiting state: {waiting_kind}', type='negative')
            return

    except Exception as e:
        error_msg = f"❌ Failed to continue task: {e}"
        with biomics_chat:
            ui.chat_message(text=error_msg, name='系统')
        ui.notify(error_msg, type='negative')
    finally:
        # 与上面一致，仅当没有新的等待状态时，认为任务结束
        if not app.storage.client.get('waiting_kind'):
            set_graph_running(False)
def reset_agent():
    """清空界面并重置会话相关的 session 状态"""

    biomics_chat.clear()
    biomics_co_pilot.clear()

    # 1. 通过 cancel_event 中断后台任务
    cancel_event = app.storage.client.get('cancel_event')
    if cancel_event:
        cancel_event.set()
        print("[DEBUG] 已设置 cancel_event，通知后台任务停止")

    # 2. 关闭会话对应的 sandbox
    sandbox_id = app.storage.client.get('sandbox_id')
    if sandbox_id:
        try:
            sandbox_manager = SandboxManager()
            sandbox_manager.close_sandbox(sandbox_id)
            print(f"[INFO] 沙箱已关闭: {sandbox_id}")
        except Exception as e:
            print(f"[WARN] 关闭沙箱失败: {e}")
    else:
        print("[INFO] 无沙箱 ID，无需关闭沙箱")
        

    # 3. 重置上传按钮和文件名标签
    upload_button.props('icon=file_upload')
    upload_button.props(remove='color')
    upload_name_label.text = ''
    app.storage.client['uploaded_file_path'] = ''
    file_upload.reset()  # 重置上传控件，允许重新上传

    # 4. 清理会话状态
    set_graph_running(False)
    app.storage.client['waiting_kind'] = None
    app.storage.client['thread_id'] = ''
    app.storage.client['graph'] = None
    app.storage.client['config'] = None
    app.storage.client['background_task'] = None
    app.storage.client['cancel_event'] = None
    app.storage.client['sandbox_id'] = None
    
    # 5. 重新显示浮动按钮栏
    floating_btn_bar.style('display: flex;')


if __name__ in {"__main__", "__mp_main__"}:
    file_upload.on_upload(handle_file_upload)
    user_input.on('keydown.enter', handle_user_input)
    reset_button.on_click(reset_agent)
    ui.run(title='Biomics Agent', favicon='/home/liyuntian/Biomics_agent/BiomicsLOGO.svg')

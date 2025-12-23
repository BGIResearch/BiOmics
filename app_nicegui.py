from ast import Dict
from time import sleep
import asyncio
import os
from dotenv import load_dotenv
from pathlib import Path
from nicegui import app, ui
from langchain_core.messages import HumanMessage
from utils.build_config_id import build_config_id
from graph.state import BrickState
from graph.builder import build_graph_with_interaction
from utils.save_dir_name import get_save_dir

# 加载环境变量
config_file = Path(__file__).parent / 'graph' / 'brick_test_config.env'
load_dotenv(dotenv_path=str(config_file))
PROJECT_ROOT = os.getenv('PROJECT_ROOT', os.path.abspath(os.path.dirname(__file__)))
"""
组件：
上方的菜单栏
左侧对话框
右侧操作台
下方输入文本框
文本框右侧重置按钮
"""


ui.query('body').style('margin: 0; padding: 0; overflow: hidden;')

# 添加 iMessage 风格的动画
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
</style>
''')

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
            ui.label('智能体运行中，请稍候...').style('color: white; font-weight: 500; font-size: 14px;')
    
    # === 顶栏 ===
    with ui.row().style(
        'height: 10%; min-height: 60px; width: 100%; '
        'align-items: center; padding: 0 24px; '
        'border-bottom: 1px solid #ddd; background: #f5f5f5;'
    ):
        ui.icon('science', size='md')
        ui.label('BRICK Agent').style('font-size: 24px; font-weight: bold; margin-left: 12px;')
        ui.space()
        ui.label('Demo Layout').style('font-size: 14px; color: #888;')

    # === 中间区域===
    with ui.element('div').style(
        'height: 75%; width: 100%; margin: 0 auto; display: flex; flex-direction: row;'
    ):
        
        # 左侧:对话区
        with ui.element('div').style(
            'width: 50%; height: 100%; padding: 16px; '
            'border-right: 1px solid #ddd; display: flex; flex-direction: column;'
        ):
            ui.label('💬 对话').style('font-size: 18px; font-weight: 600; margin-bottom: 8px;')
            biomics_chat = ui.scroll_area().style('width: 100%; flex: 1;')


        
        # 右侧:代码区
        with ui.element('div').style(
            'width: 50%; height: 100%; padding: 16px; display: flex; flex-direction: column;'
        ):
            ui.label('🧾 代码').style('font-size: 18px; font-weight: 600; margin-bottom: 8px;')
            biomics_co_pilot = ui.scroll_area().style('width: 100%; flex: 1;')


    # === 底部对话栏 ===
    with ui.row().style(
        'height: 10%; min-height: 60px; width: 100%; '
        'align-items: center; padding: 0 24px; gap: 12px; '
        'border-top: 1px solid #ddd;'
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
        user_input = ui.input(placeholder='请输入消息或指令...').style('flex: 1;')
        
        # 右侧:重置按钮
        reset_button = ui.button('重置', icon='restart_alt').props('outlined')

"""
函数：
1. 对话区，根据event更新向对话区添加信息
2. 代码区，根据event更新右侧的copilot
3. 重置页面按钮对应函数，清除页面内容
4. 

"""
UPLOAD_DIR = os.path.join(PROJECT_ROOT, 'data')

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
    """上传事件回调：保存文件并更新图标旁文字"""
    
    # 如果图正在运行，禁止上传
    if app.storage.client.get('graph_running', False):
        ui.notify('当前任务正在执行，暂时不能上传文件', type='warning')
        return
    
    os.makedirs(UPLOAD_DIR, exist_ok=True)
    
    # NiceGUI 的 SmallFileUpload 对象，read() 是异步方法
    file_name = e.file.name
    file_content = await e.file.read()  # 必须 await
    
    # 保存文件
    save_path = os.path.join(UPLOAD_DIR, file_name)
    with open(save_path, 'wb') as f:
        f.write(file_content)
    
    # 更新右侧的小文字为文件名
    upload_name_label.text = os.path.basename(save_path)
    # 把路径存到 client storage，后续 start_graph 时可以使用
    app.storage.client['uploaded_file_path'] = save_path
    # 图标变成绿色打勾，提示上传成功
    upload_button.props('icon=check_circle color=positive')
    ui.notify(f'文件已上传: {os.path.basename(save_path)}', type='positive')
def agent_update_chat(event) -> None:
    """
    supervisor
    env_checker
    general_responder
    data_analyzer
    analyze_planner
    planner
    plan_executor
    coder
    code_runner
    code_debugger
    responder 
    notebook_searcher
    """

    
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
            ui.chat_message(text=[agent_thought, agent_output], name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="data_analyzer":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=[agent_thought, agent_output], name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="analyze_planner":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=[agent_thought, agent_output], name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="planner":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="plan_executor":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="coder":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=agent_thought, name=agent_name).style(message_animation)
            with ui.chat_message(name=agent_name):
                ui.code(agent_output)
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
            with ui.chat_message(name=agent_name):
                ui.code(agent_output)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="responder":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=[agent_thought, agent_output], name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
    elif agent_name=="notebook_searcher":
        pass
    elif agent_name=="general_responder":
        agent_thought = event.get('thought')
        agent_output = event.get('output')
        with biomics_chat:
            ui.chat_message(text=[agent_thought, agent_output], name=agent_name).style(message_animation)
        biomics_chat.scroll_to(percent=1.0)
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
def process_events(graph, events, config):
    """循环消费 events，直到遇到中断或 FINISHED
    返回值: (waiting_kind, graph) 元组
    - waiting_kind: None / 'update_data_info' / 'update_data_repo'
    - graph: 更新后的 graph 实例
    """
    for event in events:
        status = event.get("status")
        agent_update_chat(event)  # 统一更新界面

        if status == "AWAITING_CONFIRMATION":
            # env_checker 需要 update_data_info：在 UI 上提示，并记录“等待类型”
            print("[DEBUG] process_events 返回 waiting_kind = 'update_data_info'")
            return ('update_data_info', graph)

        elif status == "Revise":
            # data_analyzer 需要 update_data_repo
            print("[DEBUG] process_events 返回 waiting_kind = 'update_data_repo'")
            return ('update_data_repo', graph)

        elif status == "VALIDATED":
            # 根据项目记忆：VALIDATED 也需要手动继续 stream
            events = graph.stream(None, config=config, stream_mode="values")
            first_event = next(events, None)
            return process_events(graph, events, config)

        elif status == "NOT_FINISHED":
            continue

        elif status == "FINISHED":
            user_question = event.get("question")
            if event.get("agent")=='general_responder':
                with biomics_chat:
                    with ui.card():
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.label(f"您的问题“{user_question}”不在我们能力范围内，试试让我们做细胞注释、基因富集等任务？")
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.button('重新提问', on_click=lambda: reset_button.run_method('click')).props('flat round dense')
            elif event.get('agent')=='responder':
                with biomics_chat:
                    with ui.card():
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.label(f"您的任务需求“{user_question}”已完成，此会话将已结束")
                        with ui.row().props('no-wrap').style('align-items: center;'):
                            ui.button('重新提问', on_click=lambda: reset_button.run_method('click')).props('flat round dense')
            return (None, graph)
    
    # 如果 events 消费完没有任何特殊状态，返回 None
    return (None, graph)
async def handle_user_input():
    """统一处理用户输入：
    - 没有等待状态时，作为新问题启动一条图
    - 等待状态为 AWAITING_CONFIRMATION / Revise 时，作为反馈或修改继续图
    """

    text = (user_input.value or '').strip()
    if not text:
        ui.notify('请输入内容', type='warning')
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

    # ============= 情况一：当前没有等待中断输入，视为“新问题” =============
    if not waiting_kind:
        if is_running:
            ui.notify('当前已有任务在执行，请稍候再提问', type='warning')
            return
        # 为本次会话生成 config（使用独立 thread_id）
        thread_id = build_config_id()
        config = {"configurable": {"thread_id": thread_id}}

        try:
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

            # 交给统一的事件处理逻辑（后台线程执行，避免阻塞UI）
            task = asyncio.create_task(asyncio.to_thread(process_events, graph, events, config))
            app.storage.client['background_task'] = task
            result = await task
            
            # 后台线程返回的状态，由主线程设置到 app.storage.client
            if result:
                waiting_kind, graph = result
                app.storage.client['graph'] = graph
                app.storage.client['config'] = config
                if waiting_kind:
                    app.storage.client['waiting_kind'] = waiting_kind
                    set_graph_running(False)
                    print(f"[DEBUG] 主线程设置 waiting_kind = {waiting_kind}")
                else:
                    set_graph_running(False)

        except Exception as e:
            error_msg = f"❌ 启动任务失败：{e}"
            with biomics_chat:
                ui.chat_message(text=error_msg, name='系统')
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
            ui.notify('内部状态丢失，请重新开始会话', type='negative')
            return


        if waiting_kind == 'update_data_info':
            print("[DEBUG] 进入 update_data_info 分支")
            graph.update_state(config=config, values={'update_data_info': text})
            app.storage.client['waiting_kind'] = None
            set_graph_running(True)
            events = graph.stream(None, config=config, stream_mode='values')
            first_event = next(events, None)
            task = asyncio.create_task(asyncio.to_thread(process_events, graph, events, config))
            app.storage.client['background_task'] = task
            result = await task
            
            # 后台线程返回的状态，由主线程设置
            if result:
                waiting_kind, graph = result
                app.storage.client['graph'] = graph
                app.storage.client['config'] = config
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
            task = asyncio.create_task(asyncio.to_thread(process_events, graph, events, config))
            app.storage.client['background_task'] = task
            result = await task
            
            # 后台线程返回的状态，由主线程设置
            if result:
                waiting_kind, graph = result
                app.storage.client['graph'] = graph
                app.storage.client['config'] = config
                if waiting_kind:
                    app.storage.client['waiting_kind'] = waiting_kind
                    set_graph_running(False)
                else:
                    set_graph_running(False)
        else:
            ui.notify(f'未知等待状态：{waiting_kind}', type='negative')
            return

    except Exception as e:
        error_msg = f"❌ 继续执行任务失败：{e}"
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

    set_graph_running(False)
    app.storage.client['waiting_kind'] = None

    # 为新的会话生成一个 thread_id（5 位数字字符串）
    app.storage.client['thread_id'] = ''


if __name__ in {"__main__", "__mp_main__"}:
    file_upload.on_upload(handle_file_upload)
    user_input.on('keydown.enter', handle_user_input)
    reset_button.on_click(reset_agent)
    ui.run()

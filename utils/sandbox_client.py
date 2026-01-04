import requests
import websocket
import json
import uuid
import os
import base64
import time
from pathlib import Path
from typing import List, Dict, Any, Optional

class DockerSandbox:

    def __init__(self, server_url="127.0.0.1:8888", workspace="/workspace"):
        self.http_url = f"http://{server_url}"
        self.ws_url = f"ws://{server_url}"
        self.kernel_id = None
        self.ws = None
        self.workspace = workspace
        self.session_id = str(uuid.uuid4())

    def start(self) -> bool:
        print(f"正在连接沙盒服务 {self.http_url}...")
        try:
            response = requests.post(f"{self.http_url}/api/kernels")
            if response.status_code == 201:
                self.kernel_id = json.loads(response.text)["id"]
                print(f"✅ 沙箱创建成功! ID: {self.kernel_id}")
                self._connect_ws()
                # 设置工作目录
                self._set_working_directory()
                return True
            else:
                raise Exception(f"创建失败: {response.text}")
        except Exception as e:
            print(f"❌ 连接失败: {e}")
            return False

    def _connect_ws(self) -> bool:
        """连接到内核的 WebSocket 通道，返回是否连接成功"""
        url = f"{self.ws_url}/api/kernels/{self.kernel_id}/channels"
        try:
            # 如已有连接先关掉，避免残留
            if getattr(self, "ws", None) is not None:
                try:
                    self.ws.close()
                except Exception:
                    pass

            self.ws = websocket.create_connection(url, timeout=60)
            return True
        except Exception as e:
            print(f"⚠️ 连接内核 WebSocket 失败: {e}")
            self.ws = None
            return False
    
    def _set_working_directory(self):
        """设置内核的工作目录"""
        code = f"import os; os.chdir('{self.workspace}')"
        self._run_silent(code)

    def _run_silent(self, code: str) -> bool:
        """静默执行代码（不打印调试信息）
        
        Args:
            code: 要执行的 Python 代码
            
        Returns:
            bool: 执行成功返回 True，失败返回 False
        """
        if not self.ws:
            return False

        msg_id = str(uuid.uuid4())
        message = {
            "header": {
                "msg_id": msg_id,
                "username": "agent",
                "session": self.session_id,
                "msg_type": "execute_request",
                "version": "5.3"
            },
            "content": {"code": code, "silent": True},
            "parent_header": {},
            "metadata": {},
            "channel": "shell"
        }

        self.ws.send(json.dumps(message))
        
        # 等待执行完成
        while True:
            try:
                response = json.loads(self.ws.recv())
                msg_type = response["msg_type"]
                content = response["content"]
                
                if msg_type == "status" and content["execution_state"] == "idle":
                    return True
                elif msg_type == "error":
                    return False
            except Exception:
                return False
    
    def run(self, code: str, verbose: bool = True, max_wait_time: int = 300) -> Dict[str, Any]:
        """运行代码并返回结果
        
        Args:
            code: 要执行的 Python 代码
            verbose: 是否打印详细信息
            max_wait_time: 最大等待时间（秒），默认 300秒
            
        Returns:
            Dict: 包含 stdout, stderr, result, images 等字段的字典
        """
        if not self.ws:
            return {"error": "沙箱未启动"}

        msg_id = str(uuid.uuid4())
        message = {
            "header": {
                "msg_id": msg_id,
                "username": "agent",
                "session": self.session_id,
                "msg_type": "execute_request",
                "version": "5.3"
            },
            "content": {"code": code, "silent": False},
            "parent_header": {},
            "metadata": {},
            "channel": "shell"
        }

        if verbose:
            print(f"\n>>> 正在执行代码 ({len(code)} chars)...")
        self.ws.send(json.dumps(message))

        # 收集结果
        result = {
            "stdout": [],
            "stderr": [],
            "result": None,
            "images": [],
            "error": None
        }
        
        # 计时器：用于总超时和 idle 判断
        start_time = time.time()
        last_idle_time = None  # 记录最后一次收到 idle 的时间
        idle_threshold = 10  # idle 后等待 N 秒仍无 execute_reply 就结束
        
        # 设置较短的 recv 超时，以便定期检查 idle 超时
        original_timeout = self.ws.gettimeout()
        self.ws.settimeout(2.0)  # 2秒超时，定期检查 idle 状态

        while True:
            try:
                response = json.loads(self.ws.recv())
                msg_type = response["msg_type"]
                content = response.get("content", {})
                parent_header = response.get("parent_header", {})
                
                # 检查总超时
                elapsed = time.time() - start_time
                if elapsed > max_wait_time:
                    result["error"] = f"执行超时(>{max_wait_time}秒)"
                    if verbose:
                        print(f"\n⚠️ 执行超时！已经过去 {elapsed:.1f} 秒")
                    break
                
                # 只处理与当前请求相关的消息
                if parent_header.get("msg_id") != msg_id and msg_type != "status":
                    continue
                
                if verbose:
                    print(f"  [Debug] 收到消息类型: {msg_type}")
                
                # 1. 标准输出
                if msg_type == "stream":
                    text = content.get('text', '')
                    if content.get('name') == 'stdout':
                        result["stdout"].append(text)
                        if verbose:
                            print(f"[输出]: {text}", end="")
                    elif content.get('name') == 'stderr':
                        result["stderr"].append(text)
                        if verbose:
                            print(f"[错误]: {text}", end="")
                
                # 2. 错误信息
                elif msg_type == "error":
                    trace_str = '\n'.join(content.get('traceback', []))
                    result["error"] = trace_str
                    if verbose:
                        print(f"\n❌ 运行时错误:\n{trace_str}")
                
                # 3. 执行结果
                elif msg_type == "execute_result":
                    result["result"] = content.get('data', {}).get('text/plain', '')
                    if verbose:
                        print(f"\n[结果]: {result['result']}")
                    
                    # 检查是否有图像数据
                    data = content.get('data', {})
                    if 'image/png' in data:
                        result["images"].append({
                            'type': 'png',
                            'data': data['image/png']
                        })

                # 4. 富文本/图表
                elif msg_type == "display_data":
                    data = content.get('data', {})
                    if verbose:
                        print(f"\n[显示数据]: {list(data.keys())}")
                    
                    # 提取图像
                    if 'image/png' in data:
                        result["images"].append({
                            'type': 'png',
                            'data': data['image/png']
                        })
                    if 'image/jpeg' in data:
                        result["images"].append({
                            'type': 'jpeg',
                            'data': data['image/jpeg']
                        })
                
                # 5. 执行回复（真正的结束标志）
                elif msg_type == "execute_reply":
                    status_str = content.get('status')
                    if verbose:
                        print(f"\n  [Debug] execute_reply 状态: {status_str}")
                        print("\n>>> 执行结束.")
                    # 直接结束，不再等待 idle
                    break

                # 6. 状态变化（用于 idle + 超时判断）
                elif msg_type == "status":
                    state = content.get("execution_state")
                    if verbose:
                        print(f"  [Debug] 内核状态: {state}")
                    
                    # 如果收到 idle，记录时间
                    if state == "idle":
                        if last_idle_time is None:
                            last_idle_time = time.time()
                            if verbose:
                                print(f"  [Debug] 内核进入 idle 状态，开始等待 execute_reply...")
                        else:
                            # 已经 idle 一段时间了，检查是否超过阈值
                            idle_duration = time.time() - last_idle_time
                            if idle_duration > idle_threshold:
                                if verbose:
                                    print(f"\n  [Debug] idle {idle_duration:.1f}秒仍无 execute_reply，判定为执行结束")
                                    print("\n>>> 执行结束（通过 idle 超时）.")
                                break
                    elif state == "busy":
                        # 如果又变 busy 了，重置 idle 计时器
                        last_idle_time = None
            
            except websocket.WebSocketTimeoutException:
                # recv 超时，检查是否已经在 idle 状态超过阈值
                if last_idle_time is not None:
                    idle_duration = time.time() - last_idle_time
                    if idle_duration > idle_threshold:
                        if verbose:
                            print(f"\n  [Debug] idle {idle_duration:.1f}秒仍无 execute_reply，判定为执行结束")
                            print("\n>>> 执行结束（通过 idle 超时）.")
                        break
                # 检查总超时
                elapsed = time.time() - start_time
                if elapsed > max_wait_time:
                    result["error"] = f"执行超时(>{max_wait_time}秒)"
                    if verbose:
                        print(f"\n⚠️ 执行超时！已经过去 {elapsed:.1f} 秒")
                    break
                # 否则继续等待
                continue
            except Exception as e:
                result["error"] = str(e)
                if verbose:
                    print(f"\n❌ 通信异常: {e}")
                break
        
        # 恢复原始超时设置
        try:
            self.ws.settimeout(original_timeout)
        except Exception:
            pass
        
        return result

    def upload_file(self, local_path: str, remote_path: Optional[str] = None) -> bool:
        """上传文件到沙箱
        
        Args:
            local_path: 本地文件路径
            remote_path: 沙箱中的目标路径（相对于 workspace），默认使用文件名
            
        Returns:
            bool: 上传成功返回 True
        """
        local_file = Path(local_path)
        if not local_file.exists():
            print(f"❌ 本地文件不存在: {local_path}")
            return False
        
        if remote_path is None:
            remote_path = local_file.name
        
        # 确保是相对路径
        remote_path = remote_path.lstrip('/')
        full_remote_path = f"{self.workspace}/{remote_path}"
        
        print(f"📤 上传文件: {local_path} -> {full_remote_path}")
        
        try:
            # 读取文件内容并 base64 编码
            with open(local_file, 'rb') as f:
                content = f.read()
            
            content_b64 = base64.b64encode(content).decode('utf-8')
            
            # 使用 Python 代码写入文件
            code = f"""
import base64
import os
from pathlib import Path

remote_path = '{full_remote_path}'
Path(os.path.dirname(remote_path)).mkdir(parents=True, exist_ok=True)

content = base64.b64decode('{content_b64}')
with open(remote_path, 'wb') as f:
    f.write(content)

print(f'文件已写入: {{remote_path}}')
"""
            result = self.run(code, verbose=False)
            
            if result.get("error"):
                print(f"❌ 上传失败: {result['error']}")
                return False
            
            print(f"✅ 文件上传成功: {remote_path}")
            return True
            
        except Exception as e:
            print(f"❌ 上传失败: {e}")
            return False
    
    def list_files(self, directory: str = ".") -> List[Dict[str, Any]]:
        """列出沙箱目录中的文件
        
        Args:
            directory: 目录路径（相对于 workspace），默认为当前目录
            
        Returns:
            List[Dict]: 文件列表，每个文件包含 name, type, size 等信息
        """
        directory = directory.lstrip('/')
        full_path = f"{self.workspace}/{directory}" if directory != "." else self.workspace
        
        code = f"""
import os
import json
from pathlib import Path

target_dir = Path('{full_path}')
if not target_dir.exists():
    print(json.dumps({{"error": "目录不存在"}}))
else:
    files = []
    for item in sorted(target_dir.iterdir()):
        files.append({{
            'name': item.name,
            'type': 'dir' if item.is_dir() else 'file',
            'size': item.stat().st_size if item.is_file() else 0,
            'path': str(item.relative_to('{self.workspace}'))
        }})
    print(json.dumps(files))
"""
        result = self.run(code, verbose=False)
        
        if result.get("error"):
            print(f"❌ 列出文件失败: {result['error']}")
            return []
        
        try:
            # 从 stdout 中提取 JSON
            output = ''.join(result.get("stdout", []))
            files = json.loads(output)
            
            if isinstance(files, dict) and "error" in files:
                print(f"❌ {files['error']}")
                return []
            
            return files
        except json.JSONDecodeError:
            print(f"❌ 解析文件列表失败")
            return []
    
    def download_file(self, remote_path: str, local_path: Optional[str] = None) -> bool:
        """从沙箱下载文件
        
        Args:
            remote_path: 沙箱中的文件路径（相对于 workspace）
            local_path: 本地保存路径，默认使用文件名
            
        Returns:
            bool: 下载成功返回 True
        """
        remote_path = remote_path.lstrip('/')
        full_remote_path = f"{self.workspace}/{remote_path}"
        
        if local_path is None:
            local_path = Path(remote_path).name
        
        print(f"📥 下载文件: {full_remote_path} -> {local_path}")
        
        code = f"""
import base64
import os

remote_path = '{full_remote_path}'
if not os.path.exists(remote_path):
    print('ERROR:FILE_NOT_FOUND')
else:
    with open(remote_path, 'rb') as f:
        content = f.read()
    print(base64.b64encode(content).decode('utf-8'))
"""
        result = self.run(code, verbose=False)
        
        if result.get("error"):
            print(f"❌ 下载失败: {result['error']}")
            return False
        
        output = ''.join(result.get("stdout", []))
        
        if output.startswith('ERROR:FILE_NOT_FOUND'):
            print(f"❌ 文件不存在: {remote_path}")
            return False
        
        try:
            # 解码 base64 内容
            content = base64.b64decode(output)
            
            # 写入本地文件
            with open(local_path, 'wb') as f:
                f.write(content)
            
            print(f"✅ 文件下载成功: {local_path}")
            return True
            
        except Exception as e:
            print(f"❌ 下载失败: {e}")
            return False
    
    def save_image(self, image_data: str, image_type: str, save_path: str) -> bool:
        """保存图像数据到本地文件
        
        Args:
            image_data: Base64 编码的图像数据
            image_type: 图像类型（png, jpeg）
            save_path: 保存路径
            
        Returns:
            bool: 保存成功返回 True
        """
        try:
            content = base64.b64decode(image_data)
            with open(save_path, 'wb') as f:
                f.write(content)
            print(f"✅ 图像已保存: {save_path}")
            return True
        except Exception as e:
            print(f"❌ 保存图像失败: {e}")
            return False
    
    def close(self):
        """关闭沙箱（删除内核）"""
        if self.kernel_id:
            try:
                requests.delete(f"{self.http_url}/api/kernels/{self.kernel_id}")
                if self.ws:
                    self.ws.close()
                print(f"✅ 沙箱已关闭 (ID: {self.kernel_id})")
            except Exception as e:
                print(f"⚠️ 关闭沙箱时出错: {e}")
    
    def get_id(self) -> Optional[str]:
        """获取沙箱 ID"""
        return self.kernel_id

if __name__ == "__main__":
    # 创建沙箱
    box = DockerSandbox(server_url="127.0.0.1:8888", workspace="/workspace")
    
    if not box.start():
        print("沙箱启动失败")
        exit(1)
    
    print(f"\n沙箱 ID: {box.get_id()}")
    try:
        # ========== 测试 1: 运行代码 ==========
        print("\n" + "="*50)
        print("测试 1: 运行简单代码")
        print("="*50)
        result = box.run("""

""")
        print(result)

    finally:
        # 关闭沙箱
        print("\n" + "="*50)
        box.close()
        print("="*50)


# docker stop my_code_sandbox
# docker run -d \
#   --name my_code_sandbox \
#   --network="host" \
#   -e PYTHONUNBUFFERED=1 \
#   -v /home/liyuntian/Biomics_agent/data:/workspace/data \
#   biomics_agent:v6 \
#   jupyter kernelgateway \
#   --KernelGatewayApp.ip=0.0.0.0 \
#   --KernelGatewayApp.port=8888 \
#   --KernelGatewayApp.auth_token="" \
#   --JupyterWebsocketPersonality.list_kernels=True \
#   --KernelManager.ip=0.0.0.0
# docker rm -f my_code_sandbox
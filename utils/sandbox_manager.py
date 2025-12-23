"""
Jupyter 沙箱管理工具
提供查看、删除 Jupyter 沙箱（内核）的功能
"""

import requests
import sys
from typing import List, Dict, Optional


class SandboxManager:
    """Jupyter 沙箱管理器"""
    
    def __init__(self, server_url: str = "http://127.0.0.1:8888"):
        """
        初始化 Jupyter 沙箱管理器
        
        Args:
            server_url: Jupyter Kernel Gateway 服务地址
        """
        self.server_url = server_url.rstrip('/')
        self.base_url = f"{self.server_url}/api/kernels"
    
    def list_sandboxes(self) -> List[Dict]:
        """
        列出所有正在运行的沙箱（Jupyter 内核）
        
        Returns:
            沙箱信息列表，每个元素包含 id, name, last_activity, execution_state, connections
        """
        try:
            response = requests.get(self.base_url, timeout=5)
            response.raise_for_status()
            
            kernels = response.json()
            return kernels
            
        except requests.exceptions.RequestException as e:
            print(f"❌ 获取沙箱列表失败: {str(e)}")
            print(f"提示: 请确保 Jupyter Kernel Gateway 正在运行于 {self.server_url}")
            return []
        except Exception as e:
            print(f"❌ 发生错误: {str(e)}")
            return []
    
    def close_sandbox(self, kernel_id: str) -> bool:
        """
        关闭（删除）指定的沙箱
        
        Args:
            kernel_id: 沙箱 ID（内核 ID）
        
        Returns:
            删除是否成功
        """
        try:
            url = f"{self.base_url}/{kernel_id}"
            response = requests.delete(url, timeout=5)
            response.raise_for_status()
            
            print(f"✅ 成功关闭沙箱: {kernel_id}")
            return True
            
        except requests.exceptions.RequestException as e:
            print(f"❌ 关闭沙箱失败: {kernel_id}")
            print(f"错误信息: {str(e)}")
            return False
        except Exception as e:
            print(f"❌ 发生错误: {str(e)}")
            return False
    
    def close_all_sandboxes(self, confirm: bool = True) -> int:
        """
        关闭所有正在运行的沙箱
        
        Args:
            confirm: 是否需要用户确认
        
        Returns:
            成功关闭的沙箱数量
        """
        sandboxes = self.list_sandboxes()
        
        if not sandboxes:
            print("📭 没有找到任何运行中的沙箱")
            return 0
        
        print(f"🔍 找到 {len(sandboxes)} 个运行中的沙箱:")
        for sb in sandboxes:
            print(f"  - ID: {sb['id']} | 名称: {sb.get('name', 'N/A')} | 状态: {sb.get('execution_state', 'N/A')}")
        
        if confirm:
            response = input(f"\n⚠️  确定要关闭这 {len(sandboxes)} 个沙箱吗? [y/N]: ")
            if response.lower() not in ['y', 'yes']:
                print("❌ 操作已取消")
                return 0
        
        success_count = 0
        for sb in sandboxes:
            if self.close_sandbox(sb['id']):
                success_count += 1
        
        print(f"\n✅ 成功关闭 {success_count}/{len(sandboxes)} 个沙箱")
        return success_count
    
    def show_sandboxes(self):
        """
        以友好的格式显示所有沙箱信息
        """
        sandboxes = self.list_sandboxes()
        
        if not sandboxes:
            print("📭 没有找到任何运行中的沙箱")
            return
        
        print(f"\n🧪 Jupyter 沙箱列表 (共 {len(sandboxes)} 个):")
        print("=" * 100)
        print(f"{'KERNEL ID':<40} {'NAME':<20} {'STATE':<15} {'CONNECTIONS':<12} {'LAST ACTIVITY':<20}")
        print("-" * 100)
        
        for sb in sandboxes:
            kernel_id = sb.get('id', 'N/A')
            name = sb.get('name', 'N/A')
            state = sb.get('execution_state', 'N/A')
            connections = str(sb.get('connections', 0))
            last_activity = sb.get('last_activity', 'N/A')
            
            # 截取最后活动时间（只显示前20个字符）
            if len(last_activity) > 20:
                last_activity = last_activity[:20]
            
            print(f"{kernel_id:<40} {name:<20} {state:<15} {connections:<12} {last_activity:<20}")
        
        print("=" * 100)


def main():
    """命令行入口"""
    # 从命令行参数获取服务器地址，默认为本地
    server_url = "http://127.0.0.1:8888"
    if "--server" in sys.argv:
        idx = sys.argv.index("--server")
        if idx + 1 < len(sys.argv):
            server_url = sys.argv[idx + 1]
    
    manager = SandboxManager(server_url=server_url)
    
    if len(sys.argv) < 2:
        print("用法:")
        print("  python sandbox_manager.py list [--server <url>]     # 列出所有沙箱")
        print("  python sandbox_manager.py close <kernel_id>         # 关闭指定沙箱")
        print("  python sandbox_manager.py close-all                 # 关闭所有沙箱")
        print("\n示例:")
        print("  python sandbox_manager.py list")
        print("  python sandbox_manager.py list --server http://127.0.0.1:8888")
        print("  python sandbox_manager.py close abc123-def456-789")
        print("  python sandbox_manager.py close-all")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "list":
        manager.show_sandboxes()
    
    elif command == "close":
        if len(sys.argv) < 3:
            print("❌ 错误: 请指定要关闭的沙箱 ID")
            print("用法: python sandbox_manager.py close <kernel_id>")
            sys.exit(1)
        
        kernel_id = sys.argv[2]
        manager.close_sandbox(kernel_id)
    
    elif command == "close-all":
        manager.close_all_sandboxes(confirm=True)
    
    else:
        print(f"❌ 未知命令: {command}")
        print("可用命令: list, close, close-all")
        sys.exit(1)


if __name__ == "__main__":
    main()

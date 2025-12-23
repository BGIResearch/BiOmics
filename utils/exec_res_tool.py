"""
Jupyter 代码执行结果解析工具

提供一组函数用于解析和格式化 Jupyter 内核执行代码后返回的结果字典。
结果字典格式：
{
    "stdout": [...],
    "stderr": [...],
    "result": ...,
    "images": [{"type": "png", "data": "base64..."}, ...],
    "error": "..."
}
"""

import os
import re
import base64
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional
from dotenv import load_dotenv

# 加载环境变量
config_file = Path(__file__).parent.parent / 'graph' / 'brick_test_config.env'
load_dotenv(dotenv_path=str(config_file))
PROJECT_ROOT = os.getenv('PROJECT_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def format_stdout(result: Dict[str, Any]) -> List[str]:
    """格式化标准输出
    
    Args:
        result: Jupyter 执行结果字典
        
    Returns:
        List[str]: 标准输出的列表，每个元素是一行输出
    """
    stdout = result.get("stdout", [])
    if not stdout:
        return []
    
    # 如果已经是列表，直接返回
    if isinstance(stdout, list):
        return stdout
    
    # 如果是字符串，按行分割
    if isinstance(stdout, str):
        return [stdout]
    
    return []


def format_stderr(result: Dict[str, Any]) -> List[str]:
    """格式化标准错误输出
    
    Args:
        result: Jupyter 执行结果字典
        
    Returns:
        List[str]: 标准错误输出的列表，每个元素是一行输出
    """
    stderr = result.get("stderr", [])
    if not stderr:
        return []
    
    # 如果已经是列表，直接返回
    if isinstance(stderr, list):
        return stderr
    
    # 如果是字符串，按行分割
    if isinstance(stderr, str):
        return [stderr]
    
    return []


def format_result(result: Dict[str, Any]) -> str:
    """格式化执行结果（最后一行表达式的返回值）
    
    Args:
        result: Jupyter 执行结果字典
        
    Returns:
        str: 执行结果字符串，如果没有结果返回空字符串
    """
    res = result.get("result")
    
    if res is None:
        return ""
    
    # 转换为字符串
    if isinstance(res, str):
        return res
    
    return str(res)


def format_images(result: Dict[str, Any], 
                  save_dir: str = None) -> List[str]:
    """格式化图像数据并保存到本地
    
    Args:
        result: Jupyter 执行结果字典
        save_dir: 图像保存目录，默认为 {PROJECT_ROOT}/data/images
        
    Returns:
        List[str]: 保存的图像文件相对路径列表（相对于项目根目录）
    """
    if save_dir is None:
        save_dir = os.path.join(PROJECT_ROOT, "data/images")
    
    images = result.get("images", [])
    if not images:
        return []
    
    # 确保保存目录存在
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    saved_paths = []
    
    for i, img in enumerate(images):
        try:
            img_type = img.get("type", "png")
            img_data = img.get("data", "")
            
            if not img_data:
                continue
            
            # 生成时间戳文件名
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]  # 毫秒级
            filename = f"image_{timestamp}_{i}.{img_type}"
            file_path = save_path / filename
            
            # 解码 base64 并保存
            content = base64.b64decode(img_data)
            with open(file_path, 'wb') as f:
                f.write(content)
            
            # 计算相对路径（相对于项目根目录）
            workspace_root = Path(PROJECT_ROOT)
            relative_path = file_path.relative_to(workspace_root)
            saved_paths.append(str(relative_path))
            
        except Exception as e:
            print(f"保存图像 {i} 时出错: {e}")
            continue
    
    return saved_paths


def format_error(result: Dict[str, Any]) -> str:
    """格式化错误信息，移除 ANSI 转义序列
    
    Args:
        result: Jupyter 执行结果字典
        
    Returns:
        str: 纯文本错误信息，如果没有错误返回空字符串
    """
    error = result.get("error")
    
    if not error:
        return ""
    
    # ANSI 转义序列的正则表达式
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    
    # 移除所有 ANSI 转义序列
    clean_error = ansi_escape.sub('', error)
    
    return clean_error.strip()


def has_error(result: Dict[str, Any]) -> bool:
    """判断代码块执行是否出现错误
    
    Args:
        result: Jupyter 执行结果字典
        
    Returns:
        bool: 如果有错误返回 True，否则返回 False
    """
    error = result.get("error")
    
    # 检查 error 字段是否有内容
    if error and isinstance(error, str) and error.strip():
        return True
    
    return False


# 示例用法
if __name__ == "__main__":
    # 测试数据
    test_result = {
        "stdout": ["Line 1\n", "Line 2\n"],
        "stderr": ["Warning: something\n"],
        "result": "42",
        "images": [
            {
                "type": "png",
                "data": "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=="
            }
        ],
        "error": "\x1b[0;31mNameError\x1b[0m: name 'x' is not defined"
    }
    
    print("=== 测试格式化函数 ===\n")
    
    print("1. stdout:", format_stdout(test_result))
    print("\n2. stderr:", format_stderr(test_result))
    print("\n3. result:", format_result(test_result))
    print("\n4. images:", format_images(test_result))
    print("\n5. error:", format_error(test_result))
    print("\n6. has_error:", has_error(test_result))
    
    # 测试无错误的情况
    test_result_no_error = {
        "stdout": ["Success\n"],
        "stderr": [],
        "result": None,
        "images": [],
        "error": None
    }
    
    print("\n=== 测试无错误情况 ===\n")
    print("has_error:", has_error(test_result_no_error))

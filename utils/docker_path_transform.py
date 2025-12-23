"""
路径转换工具
将本地路径转换为 Docker 容器内路径
"""
import os
from pathlib import Path
from dotenv import load_dotenv

# 加载环境变量
config_file = Path(__file__).parent.parent / 'graph' / 'brick_test_config.env'
load_dotenv(dotenv_path=str(config_file))
PROJECT_ROOT = os.getenv('PROJECT_ROOT', os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def transform_path(path: str) -> str:
    """
    将本地路径转换为 Docker 容器内路径
    
    将字符串中的 {PROJECT_ROOT}/data 替换为 /workspace/data
    
    Args:
        path: 本地路径字符串，例如 "{PROJECT_ROOT}/data/pp.txt"
    
    Returns:
        Docker 容器内路径字符串，例如 "/workspace/data/pp.txt"
    

    """
    # 定义替换规则
    local_prefix = os.path.join(PROJECT_ROOT, "data")
    docker_prefix = "/workspace/data"
    
    # 如果路径包含本地前缀，则替换为 Docker 前缀
    if path.startswith(local_prefix):
        return path.replace(local_prefix, docker_prefix, 1)
    
    return path


def transform_path_reverse(path: str) -> str:
    """
    将 Docker 容器内路径转换为本地路径（反向转换）
    
    将字符串中的 /workspace/data 替换为 {PROJECT_ROOT}/data
    
    Args:
        path: Docker 容器内路径字符串，例如 "/workspace/data/pp.txt"
    
    Returns:
        本地路径字符串，例如 "{PROJECT_ROOT}/data/pp.txt"
    
    Examples:
        >>> transform_path_reverse("/workspace/data/pp.txt")
        '{PROJECT_ROOT}/data/pp.txt'
        
        >>> transform_path_reverse("/workspace/data/images/test.png")
        '{PROJECT_ROOT}/data/images/test.png'
    """
    # 定义替换规则
    docker_prefix = "/workspace/data"
    local_prefix = os.path.join(PROJECT_ROOT, "data")
    
    # 如果路径包含 Docker 前缀，则替换为本地前缀
    if path.startswith(docker_prefix):
        return path.replace(docker_prefix, local_prefix, 1)
    
    return path


if __name__ == "__main__":
    # 测试代码
    print(f"PROJECT_ROOT: {PROJECT_ROOT}")
    print("\n测试 transform_path:")
    test_cases = [
        os.path.join(PROJECT_ROOT, "data/pp.txt"),
        os.path.join(PROJECT_ROOT, "data/images/test.png"),
        os.path.join(PROJECT_ROOT, "data/adata_new1.h5ad"),
        "/other/path/file.txt",
        "/workspace/data/file.txt"
    ]
    
    for test_path in test_cases:
        result = transform_path(test_path)
        print(f"  {test_path}")
        print(f"  -> {result}")
        print()
    
    print("\n测试 transform_path_reverse:")
    reverse_test_cases = [
        "/workspace/data/pp.txt",
        "/workspace/data/images/test.png",
        "/workspace/data/adata_new1.h5ad",
        "/other/path/file.txt"
    ]
    
    for test_path in reverse_test_cases:
        result = transform_path_reverse(test_path)
        print(f"  {test_path}")
        print(f"  -> {result}")
        print()

import os
import yaml
from typing import Dict, Any
from dotenv import load_dotenv
import json


def replace_env_vars(value: str) -> str:
    """Replace environment variables in string values."""
    if not isinstance(value, str):
        return value
    if value.startswith("$"):
        env_var = value[1:]
        return os.getenv(env_var, value)
    return value


def process_dict(config: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively process dictionary to replace environment variables."""
    result = {}
    for key, value in config.items():
        if isinstance(value, dict):
            result[key] = process_dict(value)
        elif isinstance(value, str):
            result[key] = replace_env_vars(value)
        else:
            result[key] = value
    return result


_config_cache: Dict[str, Dict[str, Any]] = {}


def load_yaml_config(file_path: str) -> Dict[str, Any]:
    """Load and process YAML configuration file."""
    # 如果文件不存在，返回{}
    if not os.path.exists(file_path):
        return {}

    # 检查缓存中是否已存在配置
    if file_path in _config_cache:
        return _config_cache[file_path]

    # 如果缓存中不存在，则加载并处理配置
    with open(file_path, "r") as f:
        config = yaml.safe_load(f)
    processed_config = process_dict(config)

    # 将处理后的配置存入缓存
    _config_cache[file_path] = processed_config
    return processed_config

def load_whole_env_config(env_file_path: str) -> Dict[str, Any]:
    """Load and process environment variables from a file."""
    # 如果文件不存在，返回{}
    if not os.path.exists(env_file_path):
        return {}
        
    # 检查缓存中是否已存在配置
    if env_file_path in _config_cache:
        return _config_cache[env_file_path]
        
    # 如果缓存中不存在，则加载并处理配置
    load_dotenv(env_file_path)
    
    # 获取所有环境变量并处理JSON类型
    config = {}
    for key, value in os.environ.items():
        if value.startswith('{') and value.endswith('}'):
            try:
                config[key] = json.loads(value)
            except json.JSONDecodeError:
                config[key] = value
        else:
            config[key] = value
    
    # 将处理后的配置存入缓存
    _config_cache[env_file_path] = config
    return config

def load_env_config(env_file_path: str) -> Dict[str, Any]:
    if not os.path.exists(env_file_path):
        return {}
        
    if env_file_path in _config_cache:
        return _config_cache[env_file_path]
    
    # 只读取.env文件中的变量
    with open(env_file_path) as f:
        env_vars = {}
        for line in f:
            if '=' in line and not line.startswith('#'):
                parts = line.split('=', 1)
                key = parts[0]
                value = parts[1].strip()
                if '://' in value:
                    env_vars[key] = value
                elif value.startswith('{') and value.endswith('}'):
                    try:
                        env_vars[key] = json.loads(value)
                    except json.JSONDecodeError:
                        env_vars[key] = value
                else:
                    env_vars[key] = value
    
    _config_cache[env_file_path] = env_vars
    return env_vars
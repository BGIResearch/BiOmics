"""
统一的配置加载模块
从 brick_test_config.env 读取配置参数
"""
import os
import json
from pathlib import Path
from dotenv import dotenv_values

# 获取项目根目录
_current_dir = Path(__file__).parent
_project_root = _current_dir.parent
_config_path = _project_root / "graph" / "brick_test_config.env"

# 加载配置
_config = dotenv_values(_config_path)

def get_kg_config():
    """获取知识图谱配置"""
    url = _config.get("KG_URL")
    if not url:
        raise ValueError("KG_URL not found in brick_test_config.env")
    
    kg_auth = _config.get("KG_AUTH")
    kg_pass = _config.get("KG_PASS")
    if not kg_auth or not kg_pass:
        raise ValueError("KG_AUTH or KG_PASS not found in brick_test_config.env")
    
    auth = (kg_auth, kg_pass)
    return url, auth

def get_llm_config(model_key="BASIC_MODEL"):
    """获取 LLM 配置
    
    Args:
        model_key: 配置键名，可选值：
            - BASIC_MODEL
            - CODE_MODEL
            - REASONING_MODEL
            - CLAUDE_MODEL
            等
    
    Returns:
        dict: 包含 api_key, base_url, model_name, temperature 等配置
    """
    model_str = _config.get(model_key)
    if not model_str:
        raise ValueError(f"{model_key} not found in brick_test_config.env")
    
    try:
        return json.loads(model_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Failed to parse {model_key} from config: {e}")

def configure_brick(model_key="BASIC_MODEL"):
    """一键配置 BRICK（知识图谱 + LLM）
    
    Args:
        model_key: LLM 配置键名，默认使用 BASIC_MODEL
    """
    import BRICK
    
    # 配置知识图谱
    url, auth = get_kg_config()
    BRICK.config(url=url, auth=auth)
    
    # 配置 LLM
    llm_config = get_llm_config(model_key)
    BRICK.config_llm(
        modeltype='ChatOpenAI',
        api_key=llm_config.get("openai_api_key"),
        base_url=llm_config.get("base_url"),
        llm_params={
            'model_name': llm_config.get("model_name"),
            'temperature': llm_config.get("temperature", 0.3)
        }
    )

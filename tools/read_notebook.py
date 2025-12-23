from typing import Optional, Dict, Any, List
from pathlib import Path
from pydantic import BaseModel
from langchain_core.tools import tool
import json

class Input(BaseModel):
    agentname: str
    notebook_path: str

class Output(BaseModel):
    status: str
    cells_text: Optional[str] = None

@tool(args_schema=Input)
def read_notebook_tool(agentname: str, notebook_path: str) -> Output:
    """
    Use this tool to read and parse Jupyter notebook (.ipynb) files
    Only extracts markdown and code cell content, skipping outputs to reduce token usage.
    
    Args:
        agentname: The name of the current agent, like "analyzer"
        notebook_path: The path to the .ipynb file to be read,must be 'str' type
        
    Return:
        output: A dict with status and cells_text containing only markdown and code content.
    """
    try:
        # 检查文件路径是否存在
        notebook_file = Path(notebook_path)
        if not notebook_file.exists():
            return Output(status="error", cells_text=f"Notebook file not found: {notebook_path}")
        
        if not notebook_file.suffix.lower() == '.ipynb':
            return Output(status="error", cells_text=f"File is not a Jupyter notebook: {notebook_path}")
        
        # 读取并解析notebook文件
        with open(notebook_file, 'r', encoding='utf-8') as f:
            notebook_content = json.load(f)
        
        # 提取markdown、代码cell的内容以及非图片输出
        cells_text_parts = []
        cells_count = 0
        
        if 'cells' in notebook_content:
            for cell in notebook_content['cells']:
                cell_type = cell.get('cell_type', 'unknown')
                
                # 处理markdown和code类型的cell
                if cell_type in ['markdown', 'code']:
                    cells_count += 1
                    
                    # 添加cell类型标识
                    cells_text_parts.append(f"\n=== Cell {cells_count} ({cell_type}) ===\n")
                    
                    # 提取cell源代码/文本内容
                    if 'source' in cell:
                        source = cell['source']
                        if isinstance(source, list):
                            source_text = ''.join(source)
                        else:
                            source_text = str(source)
                        cells_text_parts.append(source_text)
                        
                        # 添加分隔符
                        cells_text_parts.append("\n")
                    
                    # 处理code cell的输出（排除图片）
                    if cell_type == 'code' and 'outputs' in cell:
                        for output in cell['outputs']:
                            output_type = output.get('output_type', '')
                            
                            # 处理stream类型输出（如stdout, stderr）
                            if output_type == 'stream':
                                stream_name = output.get('name', '')
                                if 'text' in output:
                                    text_content = output['text']
                                    if isinstance(text_content, list):
                                        text_content = ''.join(text_content)
                                    cells_text_parts.append(f"[{stream_name} output]:\n{text_content}\n")
                            
                            # 处理execute_result和display_data类型输出（排除图片）
                            elif output_type in ['execute_result', 'display_data']:
                                if 'data' in output:
                                    data = output['data']
                                    # 优先处理text/plain输出
                                    if 'text/plain' in data:
                                        plain_text = data['text/plain']
                                        if isinstance(plain_text, list):
                                            plain_text = ''.join(plain_text)
                                        cells_text_parts.append(f"[output]:\n{plain_text}\n")
                                    # 处理text/html输出（如果没有text/plain）
                                    elif 'text/html' in data and not any(key.startswith('image/') for key in data.keys()):
                                        html_text = data['text/html']
                                        if isinstance(html_text, list):
                                            html_text = ''.join(html_text)
                                        cells_text_parts.append(f"[html output]:\n{html_text}\n")
                            
                            # 处理error类型输出
                            elif output_type == 'error':
                                if 'traceback' in output:
                                    traceback_text = '\n'.join(output['traceback'])
                                    cells_text_parts.append(f"[error]:\n{traceback_text}\n")
        
        # 合并所有文本内容
        combined_text = '\n'.join(cells_text_parts)
        
        return Output(
            status="success",
            cells_text=combined_text
        )
        
    except json.JSONDecodeError as e:
        return Output(status="error", cells_text=f"Invalid JSON format in notebook file: {str(e)}")
    except Exception as e:
        return Output(status="error", cells_text=f"Error reading notebook: {str(e)}")


# 内部函数，直接处理notebook读取逻辑
def _read_notebook_internal(notebook_path: str) -> Output:
    """
    内部函数：直接处理notebook读取逻辑，避免LangChain工具调用问题
    只提取markdown和代码cell内容，跳过输出以减少token使用
    """
    try:
        # 检查文件路径是否存在
        notebook_file = Path(notebook_path)
        if not notebook_file.exists():
            return Output(status="error", cells_text=f"Notebook file not found: {notebook_path}")
        
        if not notebook_file.suffix.lower() == '.ipynb':
            return Output(status="error", cells_text=f"File is not a Jupyter notebook: {notebook_path}")
        
        # 读取并解析notebook文件
        with open(notebook_file, 'r', encoding='utf-8') as f:
            notebook_content = json.load(f)
        
        # 只提取markdown和代码cell的内容
        cells_text_parts = []
        cells_count = 0
        
        if 'cells' in notebook_content:
            for cell in notebook_content['cells']:
                cell_type = cell.get('cell_type', 'unknown')
                
                # 只处理markdown和code类型的cell
                if cell_type in ['markdown', 'code']:
                    cells_count += 1
                    
                    # 添加cell类型标识
                    cells_text_parts.append(f"\n=== Cell {cells_count} ({cell_type}) ===\n")
                    
                    # 提取cell源代码/文本内容
                    if 'source' in cell:
                        source = cell['source']
                        if isinstance(source, list):
                            source_text = ''.join(source)
                        else:
                            source_text = str(source)
                        cells_text_parts.append(source_text)
                        
                        # 添加分隔符
                        cells_text_parts.append("\n")
        
        # 合并所有文本内容
        combined_text = '\n'.join(cells_text_parts)
        
        return Output(
            status="success",
            cells_text=combined_text
        )
        
    except json.JSONDecodeError as e:
        return Output(status="error", cells_text=f"Invalid JSON format in notebook file: {str(e)}")
    except Exception as e:
        return Output(status="error", cells_text=f"Error reading notebook: {str(e)}")


# 便捷函数，用于直接获取notebook的文本内容
def extract_notebook_text(notebook_path: str) -> str:
    """
    便捷函数：直接提取notebook的所有文本内容
    
    Args:
        notebook_path: notebook文件路径
        
    Returns:
        提取的文本内容，如果出错则返回空字符串
    """
    result = _read_notebook_internal(notebook_path)
    if result.status == "success" and result.cells_text:
        return result.cells_text
    return ""


# 便捷函数，用于获取notebook的基本信息
def get_notebook_info(notebook_path: str) -> Dict[str, Any]:
    """
    便捷函数：获取notebook的基本信息
    
    Args:
        notebook_path: notebook文件路径
        
    Returns:
        包含notebook基本信息的字典
    """
    result = _read_notebook_internal(notebook_path)
    if result.status == "success":
        # 计算markdown和代码cell的数量
        cell_count = 0
        if result.cells_text:
            cell_count = result.cells_text.count("=== Cell")
        
        return {
            "status": "success",
            "cells_count": cell_count,
            "file_path": notebook_path,
            "has_content": bool(result.cells_text)
        }
    else:
        return {
            "status": "error",
            "message": result.cells_text,  # 错误信息现在存储在cells_text中
            "file_path": notebook_path
        }
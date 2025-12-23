from langchain_core.tools import tool
import re

@tool
def extract_python_blocks_tool(text: str) -> str:
    """
    Use this tool to extract all Python code blocks from the given text.

    Args:
        text: The input string containing potential Python code blocks.

    Returns:
        A string of cleaned Python code. Returns the original text if no markdown
        code blocks are found (assuming it's already pure Python code).
    """
    separator = "\n"
    pattern = r"```python\s*([\s\S]*?)\s*```"
    code_blocks = re.findall(pattern, text, re.DOTALL)
    
    if code_blocks:
        # 如果找到markdown代码块，提取并返回
        cleaned_blocks = [
            block.strip()
            for block in code_blocks
            if block.strip()
        ]
        return separator.join(cleaned_blocks) if cleaned_blocks else ""
    else:
        # 如果没有找到markdown代码块，假设整个文本就是Python代码
        return text.strip()
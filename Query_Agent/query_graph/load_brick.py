import json
import re
from langchain_core.tools import tool


@tool
def load_BRICK_info_tool(brick_path):
    """
    Use this tool to load the BRICK information from the specified file path. This file contains all the useful information required for using BRICK.
    Return the BRICK information in list format.
    """
    with open(brick_path, 'r', encoding='utf-8') as f:
        text = f.read()

    pattern = r"# Function:\s*'''(.*?)'''"
    matches = re.findall(pattern, text, re.DOTALL)

    info_list = []

    for match in matches:
        match = match.strip()
        try:
            data = json.loads(match)
            info_list.append(data)
        except json.JSONDecodeError as e:
            print(f"解析错误: {e}")
            continue
    return info_list

if __name__ == "__main__":
    result = load_BRICK_info_tool.invoke({"brick_path": "/home/lyt/checker_finallap/BRICK说明_v2.txt"})
    print(result)

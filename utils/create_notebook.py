"""
生成 Jupyter Notebook 文件的工具函数
"""
import nbformat
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell, new_output
import os


def create_notebook(cells_data: list, save_path: str) -> str:
    """
    根据 cells_data 生成 ipynb 文件
    
    Args:
        cells_data: 包含代码执行记录的列表，每个元素结构为:
            {
                "mdata": str,      # 步骤标题，如 "Step1:Data Loading."
                "thought": str,    # coder 的思考过程
                "code": str,       # 执行的代码
                "outputs": dict,   # 执行结果，包含 stdout, stderr, result, images, error
            }
        save_path: ipynb 文件保存路径
    
    Returns:
        保存的文件路径
    """
    nb = new_notebook()
    nb.metadata['kernelspec'] = {
        'display_name': 'Python 3',
        'language': 'python',
        'name': 'python3'
    }
    
    for cell in cells_data:
        mdata = cell.get('mdata', '')
        thought = cell.get('thought', '')
        code = cell.get('code', '')
        outputs = cell.get('outputs', {})
        
        # 添加 markdown cell（标题 + 思考）
        if mdata or thought:
            md_content = ""
            if mdata:
                md_content += f"## {mdata}\n\n"
            if thought:
                md_content += f"{thought}\n"
            nb.cells.append(new_markdown_cell(source=md_content.strip()))
        
        # 添加 code cell
        if code:
            code_cell = new_code_cell(source=code)
            code_cell.outputs = []
            
            if outputs:
                # stdout
                stdout = outputs.get('stdout', [])
                if stdout:
                    stdout_text = ''.join(stdout) if isinstance(stdout, list) else str(stdout)
                    code_cell.outputs.append(nbformat.v4.new_output(
                        output_type='stream',
                        name='stdout',
                        text=stdout_text
                    ))
                
                # stderr
                stderr = outputs.get('stderr', [])
                if stderr:
                    stderr_text = ''.join(stderr) if isinstance(stderr, list) else str(stderr)
                    code_cell.outputs.append(nbformat.v4.new_output(
                        output_type='stream',
                        name='stderr',
                        text=stderr_text
                    ))
                
                # result (execute_result)
                result_val = outputs.get('result')
                if result_val:
                    code_cell.outputs.append(nbformat.v4.new_output(
                        output_type='execute_result',
                        data={'text/plain': str(result_val)},
                        execution_count=1
                    ))
                
                # images
                images = outputs.get('images', [])
                for img in images:
                    img_type = img.get('type', 'png')
                    img_data = img.get('data', '')
                    if img_data:
                        mime_type = f'image/{img_type}'
                        code_cell.outputs.append(nbformat.v4.new_output(
                            output_type='display_data',
                            data={mime_type: img_data}
                        ))
                
                # error
                error = outputs.get('error')
                if error:
                    code_cell.outputs.append(nbformat.v4.new_output(
                        output_type='error',
                        ename='Error',
                        evalue=str(error),
                        traceback=[str(error)]
                    ))
            
            nb.cells.append(code_cell)
    
    # 确保目录存在
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    
    # 保存 notebook
    with open(save_path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
    
    return save_path

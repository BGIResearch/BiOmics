import os
import subprocess
import threading
from datetime import datetime
from typing import Union


def run_subprocess(work_dir, *exe_args):
    """运行子进程并捕获输出"""
    def read_stream(stream, output_list):
        for line in iter(stream.readline, b""):
            try:
                decoded_line = line.decode("utf-8").strip()
            except UnicodeDecodeError:
                # 尝试GBK编码
                decoded_line = line.decode("gbk", errors="replace").strip()
            output_list.append(decoded_line)

    process = subprocess.Popen(
        args=[*exe_args],
        cwd=work_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
        bufsize=-1,
        universal_newlines=False,
    )

    stdout_lines = []
    stderr_lines = []
    
    # 启动两个线程分别读取 stdout 和 stderr
    threading.Thread(
        target=read_stream, args=(process.stdout, stdout_lines), daemon=True
    ).start()
    threading.Thread(
        target=read_stream, args=(process.stderr, stderr_lines), daemon=True
    ).start()

    process.wait()  # 等待子进程结束

    return {
        "returncode": process.returncode,
        "stdout": "\n".join(stdout_lines),
        "stderr": "\n".join(stderr_lines),
    }


def write_file(code, env, work_dir):
    """将代码写入文件"""
    # 确保工作目录存在
    os.makedirs(work_dir, exist_ok=True)
    
    # 生成文件名
    formatted_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    suffix = ".r" if env == "r" else ".py"
    fname = formatted_time + suffix
    file_path = os.path.join(work_dir, fname)
    
    # 写入文件
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(code)
    
    return fname


def debug(step, code, work_dir='output', env='python', python_path='python', r_path='R', rag_result=None):
    """
    简化版的debug函数，只执行代码并返回结果，不进行错误处理
    
    参数:
        step: 当前执行步骤的描述
        code: 要执行的代码字符串
        work_dir: 工作目录，默认为'output'
        env: 执行环境，支持'python'和'r'
        python_path: Python解释器路径
        r_path: R解释器路径
        rag_result: RAG检索结果（可选，此版本中未使用）
    
    返回:
        tuple: (success, final_code, output)
            - success: 布尔值，表示执行是否成功
            - final_code: 最终的代码版本（与输入相同）
            - output: 执行输出或错误信息
    """
    
    # 如果work_dir为None，创建默认目录
    if work_dir is None:
        work_dir = os.path.join("output", datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    
    try:
        # 1. 写入文件
        print(f"正在将代码写入文件...")
        file_name = write_file(code, env, work_dir)
        file_path = os.path.join(work_dir, file_name)
        print(f"代码已写入: {file_path}")
        
        # 2. 运行代码
        print(f"正在执行代码...")
        if env == "r":
            result = run_subprocess(work_dir, r_path, "<", file_name)
        elif env == "python":
            result = run_subprocess(work_dir, python_path, file_name)
        else:
            raise ValueError(f"不支持的环境: {env}")
        
        # 3. 处理结果
        success = result["returncode"] == 0
        
        if success:
            print("代码执行成功!")
            output = result["stdout"]
        else:
            print("代码执行失败!")
            print("错误信息:")
            print(result["stderr"])
            output = result["stderr"]
        
        return success, code, output
        
    except Exception as e:
        print(f"执行过程中发生异常: {str(e)}")
        return False, code, str(e)


if __name__ == "__main__":
    # 测试代码
    test_code = """
import os
print("Hello, World!")
print("当前工作目录:", os.getcwd())
print("测试完成")
"""
    
    success, final_code, output = debug("测试步骤", test_code)
    print("\n" + "="*50)
    print("测试结果:")
    print(f"成功: {success}")
    print(f"输出: {output}")
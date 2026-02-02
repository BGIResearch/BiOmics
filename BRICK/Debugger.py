from langchain_openai import ChatOpenAI
from typing_extensions import TypedDict
from IPython.display import Image, display
from langgraph.graph import StateGraph, START, END
from datetime import datetime
import os
import subprocess
import threading
from typing import Any, Dict, List, Literal
import subprocess
import logging
import inspect
import time

def timed_node(name: str):
    def decorator(fn):
        def wrapper(*args, **kwargs):
            start = time.time()
            result = fn(*args, **kwargs)
            duration = time.time() - start
            print(f"[time] Node '{name}' took {duration:.3f} seconds.")
            return result
        return wrapper
    return decorator

class State(TypedDict):
    code: list[str]
    run_msg: list[object]
    error_msg: list[str]
    solutions: list[str]
    success: bool

    work_dir: str
    file_path: str
    env: str
    python_path: str
    r_path: str


def get_func_name():
    frame = inspect.currentframe()
    # 获取当前帧的上一级帧（即调用当前函数的帧）
    caller_frame = frame.f_back
    res = ""
    # 获取类名
    if "self" in caller_frame.f_locals:
        res = caller_frame.f_locals["self"].__class__.__name__ + "."
    # 获取当前函数的名字
    function_name = caller_frame.f_code.co_name
    return res + function_name


def run_subprocess(*exe_args):

    def read_stream(stream, output_list):
        for line in iter(stream.readline, b""):  # 注意b''表示字节模式
            try:
                decoded_line = line.decode("utf-8").strip()
            except UnicodeDecodeError:
                # 尝试GBK编码
                decoded_line = line.decode("gbk", errors="replace").strip()
            # print(f"{decoded_line}")
            logger.info(decoded_line)
            output_list.append(decoded_line)

    process = subprocess.Popen(
        [*exe_args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,  # 单独捕获 stderr
        text=False,
        bufsize=-1,
        universal_newlines=False,  # 禁用自动换行转换
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


class WriteFile:

    def __call__(self, state: State) -> State:
        logger.info("=" * 40 + get_func_name() + "=" * 40)
        self.work_dir = state["work_dir"]
        code = state.get("code", [])[-1]
        env = state["env"]
        file_path = self.write_file(code, env)
        logger.info("写入文件：" + file_path)
        return {"file_path": file_path}

    def write_file(self, code, env):
        formatted_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        suffix = ".r" if env == "r" else ".py"
        new_file = os.path.join(self.work_dir, formatted_time + suffix)
        with open(new_file, "w", encoding="utf-8") as f:
            f.write(code)
        return new_file


class RunCode:
    def __call__(self, state: State) -> State:
        logger.info("=" * 40 + get_func_name() + "=" * 40)
        env = state["env"]
        file_path = state["file_path"]
        if env == "r":
            res = run_subprocess(state["r_path"], '<', file_path)
        elif env == "python":
            res = run_subprocess(state["python_path"], file_path)
        return {"run_msg": state.get("run_msg", []) + [res], "env": env}


def run_route(state: State) -> Literal["ErrorSolutionGen", "Success"]:
    logger.info("=" * 40 + get_func_name() + "=" * 40)
    res = state.get("run_msg", [])[-1]
    # 判断成功运行
    if res["returncode"]:
        logger.info("run error!")
        return "ErrorSolutionGen"
    else:  # returncode 0表示进程成功完成没有错误
        logger.info("run success!")
        return "Success"


def success(state: State):
    return {"success": True}

@timed_node("package_install")
def package_install(state: State):
    logger.info("=" * 40 + get_func_name() + "=" * 40)
    solution = state["solutions"][-1]
    env = state["env"]
    cmd = model.invoke(
        f"""### 任务说明：当前环境是{env}，请根据解决方案，如果是r环境请给出安装语句，如果是python环境只输出包的名字。
### 解决方案：{solution}
### 案例：
R语言：
install.packages(c("pkg1", "pkg2", "pkg3"), repos = "https://cloud.r-project.org/", quiet = TRUE)
install.packages(c("pkg1", "pkg2", "pkg3"), repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/", quiet = TRUE)
install.packages("BiocManager", repos = "https://cloud.r-project.org/", quiet = TRUE); options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor"); BiocManager::install(c("pkg1", "pkg2", "pkg3"))
install.packages("devtools", repos = "https://cloud.r-project.org/", quiet = TRUE); devtools::install_github(c("hadley/ggplot2", "tidyverse/dplyr"))
install.packages("remotes", repos = "https://cloud.r-project.org/", quiet = TRUE); remotes::install_github(c("hadley/ggplot2", "tidyverse/dplyr"))
---
Python语言：
pkg1, pkg2
### 输出要求：请以普通字符串的形式仅输出命令语句，不要输出任何额外的符号如```等，不要输出解释信息。"""
    ).content
    if env == "python":
        res = run_subprocess(state["python_path"], "-m", "pip", "install", cmd)
    elif env == "r":
        res = run_subprocess(state["r_path"], "<", cmd)
    logger.info(res)
    # if res.returncode:
    #     print('install error!')
    # else: # returncode 0表示进程成功完成没有错误
    #     print('install success!')

@timed_node("input_recollect")
def input_recollect(state: State):
    logger.info("=" * 40 + get_func_name() + "=" * 40)
    str = input("请根据以下意见给出修改：" + state["solutions"][-1])
    code = state["code"][-1]
    code = model.invoke(
        f"""<p>任务说明：请根据用户的输入信息对代码进行必要的修改：</p>
<p>用户输入信息：{str}</p>
<p>原代码：<code>{code}</code></p>
<p>输出要求：请以普通字符串的形式只输入修改后的代码本身，不要输出任何额外的符号如```等，不要输出任何格式信息。</p>"""
    ).content
    return {"code": state["code"] + [code]}

@timed_node("error_revisor")
def error_revisor(state: State):
    logger.info("=" * 40 + get_func_name() + "=" * 40)
    code = model.invoke(
        f"""#任务说明：请根据给出的修改方案，对以下代码进行修改，并返回修改后的代码：
------
#解决方案：
{state.get("solutions")[-1]}
------   
#待修改的代码：
{state.get("code")[-1]}
------
#输出要求：请确保输出内容只包含修改后的代码本身，不要输出额外的符号如```和别的格式信息。"""
    ).content
    return {"code": state.get("code", []) + [code]}

@timed_node("_check_if_install")
def _check_if_install(solution):
    res = model.invoke(
        f"""请判断解决方案中是否明确说明需要安装包，只输出yes/no这几个字符，不输出思考过程\n\n{solution}"""
    ).content
    logger.info(res)
    return "yes" in res.lower()

@timed_node("_check_if_input_miss")
def _check_if_input_miss(solution):
    res = model.invoke(
        f"""请根据以下信息判断是否可以只通过修改代码来修正错误，如果可以修正则输出no；
除非很大概率可以确定用户输入的文件或提供的信息导致了错误，需要通过用户重新输入信息来修正错误，则输出yes。回答只输出yes或no这个字符串，不输出别的任何内容\n\n{solution}"""
    ).content
    logger.info(res)
    return "yes" in res.lower()

@timed_node("error_solution_gen")
def error_solution_gen(state: State):
    logger.info("=" * 40 + get_func_name() + "=" * 40)
    run_msg = state["run_msg"][-1]
    error_msg = model.invoke(
        f"""请从如下输出信息中，抽取最关键的错误内容、最先报错的用户代码段、建议解决方案，请保证抽取的内容为原文：\n\n{run_msg['stderr']} """
    ).content
    logger.error("-" * 40 + "报错总结：" + "-" * 40)
    logger.error(error_msg)
    # search_str = model.invoke(f"请从以下信息中抽取最关键的信息以进行搜索引擎搜索，包含最最关键的错误信息和用户代码。只返回该问句。").content
    # google_res = search_google(search_str, num_result_pages=2, search_language="zh")

    solution = model.invoke(
        "请根据以下报错信息给出最可能的**一个**修改意见或者解决方案，要求简洁明了：\n"
        + error_msg
    ).content
    logger.error("-" * 40 + "解决方案：" + "-" * 40)
    logger.error(solution)
    return {
        "solutions": state.get("solutions", []) + [solution],
        "error_msg": state.get("error_msg", []) + [error_msg],
    }


def solution_route(
    state: State,
) -> Literal["PackageInstall", "ErrorRevisor", "InputRecollect", "__end__"]:
    logger.info("=" * 40 + get_func_name() + "=" * 40)
    if len(state["solutions"]) > 20:
        logger.info("出错20次了，不干了")
        return "__end__"

    solution = state["solutions"][-1]

    if _check_if_install(solution):
        return "PackageInstall"
    if _check_if_input_miss(solution):
        return "InputRecollect"
    else:
        return "ErrorRevisor"


def get_logger(work_dir: str):
    formatted_time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    logger = logging.getLogger("executor_logger")
    logger.setLevel(logging.INFO)  # 设置日志级别

    # 创建一个文件处理器，用于写入日志文件
    fname = os.path.join(work_dir, formatted_time + ".log")
    dirname = os.path.dirname(fname)
    if dirname and not os.path.exists(dirname):  # 检查目录是否存在
        os.makedirs(dirname, exist_ok=True)
    file_handler = logging.FileHandler(fname, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)  # 设置文件处理器的日志级别
    file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_formatter)

    # 创建一个控制台处理器，用于输出到控制台
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)  # 设置控制台处理器的日志级别
    console_formatter = logging.Formatter("%(levelname)s - %(message)s")
    console_handler.setFormatter(console_formatter)

    # 将处理器添加到 Logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

def model_init():
    # TODO 
    deepseek = ChatOpenAI(
        model='deepseek-v3-hs', # 30s
        # model='deepseek-r1-hs', # 14.9s
        api_key="sk-kpsteSkpDGl1xBmDEcC7D51b968e43499092826f17286b55",
        base_url='http://10.224.28.80:3000/v1',
    )
    return deepseek

@timed_node("debug")
def debug(code, work_dir='output', env='python', python_path='python', r_path='R'):
    global model, logger
    model = model_init()
    logger = get_logger(work_dir)

    state = {
        "code": [code],
        "work_dir": work_dir,
        "env": env,
        "python_path": python_path,
        "r_path": r_path,
        "success": False
    }

    builder = StateGraph(State)
    builder.add_node("WriteFile", WriteFile())
    builder.add_node("RunCode", RunCode())
    builder.add_node("ErrorSolutionGen", error_solution_gen)
    builder.add_node("PackageInstall", package_install)
    builder.add_node("ErrorRevisor", error_revisor)
    builder.add_node("InputRecollect", input_recollect)
    builder.add_node("Success", success)

    # Logic
    builder.add_edge(START, "WriteFile")
    builder.add_edge("WriteFile", "RunCode")
    builder.add_conditional_edges("RunCode", run_route)
    builder.add_conditional_edges("ErrorSolutionGen", solution_route)
    builder.add_edge("ErrorRevisor", "WriteFile")
    builder.add_edge("InputRecollect", "WriteFile")
    builder.add_edge("PackageInstall", "RunCode")
    builder.add_edge("Success", END)

    # Add
    graph = builder.compile()

    # View
    # display(Image(graph.get_graph().draw_mermaid_png()))
    state = graph.invoke(state, {"recursion_limit": 100})
    return state["success"], state["code"][-1], state["run_msg"][-1]["stdout"]

if __name__ == '__main__':
    # with open(r'D:\1Grad\Code\GraAgent\brick\test\complete.py', 'r', encoding='utf-8') as f:
    #     code = f.read()
    with open(r'1.py', 'r', encoding='utf-8') as f:
        code = f.read()
    state, code, output = debug(code)
    print('state', state)
    print('code', code)
    print('output', output)

import os
import ast
import glob
from typing import Literal, Union, List, ClassVar,Dict, Any, Optional
from langgraph.graph import StateGraph, START, END
from pydantic import BaseModel
import json
from IPython.display import Image, display
from langchain_community.vectorstores import FAISS
from pathlib import Path
from rapidfuzz import fuzz
from langchain.docstore.document import Document

from config import model_v3, output_parser, str_output_parser, model_q, embedding_model
#from utils import timed_node
import logging
#from logger import setup_logger
#logger = setup_logger('rag', 'logs/rag.log', level=logging.INFO)

#@timed_node("notebook_creator_cell")
def notebook_creator_cell(file_path, include_outputs=False, max_output_length=20, remove_newline=True, db_name="BRICK_notebook", embedding_model=None):
    if os.path.isdir(file_path):
        file_path = os.path.join(file_path, "*.ipynb")

    files = glob.glob(file_path)
    if not files:
        raise FileNotFoundError(f"No .ipynb files found in the specified path: {file_path}")

    documents = []

    for file in files:
        with open(file, "r", encoding="utf-8") as f:
            notebook = json.load(f)

        for idx, cell in enumerate(notebook.get("cells", [])):
            cell_type = cell.get("cell_type", "")
            source = ''.join(cell.get("source", []))  # 合并多行 source
            
            if remove_newline:
                source = source.replace('\n', ' ').strip()

            # Optionally include output
            output_text = ""
            if include_outputs and cell_type == "code":
                outputs = cell.get("outputs", [])
                for output in outputs:
                    if "text" in output:
                        output_text += ''.join(output["text"])
                    elif "data" in output and "text/plain" in output["data"]:
                        output_text += ''.join(output["data"]["text/plain"])

                if max_output_length:
                    output_text = output_text[:max_output_length]
                if remove_newline:
                    output_text = output_text.replace('\n', ' ').strip()
                if output_text:
                    source += f"\n[Output]: {output_text}"

            metadata = {
                "source_file": os.path.basename(file),
                "cell_index": idx,
                "cell_type": cell_type
            }

            documents.append(Document(page_content=source, metadata=metadata))

    db = FAISS.from_documents(documents, embedding_model)
    db.save_local(db_name)
    return db

#@timed_node("extract_functions_name_from_py")
def extract_functions_name_from_py(file_path):
    """
    Extracts all functions and their docstrings from a Python (.py) file.

    Args:
        file_path (str): The path to the Python file.

    Returns:
        list: A list of dictionaries, each representing a function with the following keys:
            - ``name`` (str): The function name.
            - ``doc`` (str): The function's docstring. If no docstring is available, returns "No docstring available."
            - ``file`` (str): The path of the source file.
    """
    with open(file_path, "r", encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=file_path)

    functions = []
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):  # 仅提取函数
            func_name = node.name
            docstring = ast.get_docstring(node) or "No docstring available."
            functions.append({"name": func_name, "doc": docstring, "file": file_path})
    return functions

#@timed_node("extract_functions_from_py")
def extract_functions_from_py(file_path):
    """
    Extracts all functions and their docstrings and source code from a Python (.py) file.

    Args:
        file_path (str): The path to the Python file.

    Returns:
        list: A list of dictionaries, each representing a function with:
            - ``name`` (str): The function name.
            - ``doc`` (str): The function's docstring.
            - ``file`` (str): The path of the source file.
            - ``source`` (str): The full source code of the function.
    """
    with open(file_path, "r", encoding="utf-8") as f:
        source_code = f.read()
        tree = ast.parse(source_code, filename=file_path)

    functions = []
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            func_name = node.name
            docstring = ast.get_docstring(node) or "No docstring available."
            try:
                func_source = ast.get_source_segment(source_code, node)
            except Exception:
                # fallback if get_source_segment fails
                start_line = node.lineno - 1
                end_line = getattr(node, 'end_lineno', start_line + 1)
                func_source = "\n".join(source_code.splitlines()[start_line:end_line])
            functions.append({
                "name": func_name,
                "doc": docstring,
                "file": file_path,
                "source": func_source
            })

    return functions

#@timed_node("extract_single_function")
def extract_single_function(file_path):
    """
    Extract functions (yield them out) one by one from a Python file. 
    Args:
        file_path (str): Python file path. 
    Yields:
        dict: Information of a single function.
    """
    with open(file_path, "r", encoding="utf-8") as f:
        source_code = f.read()
        tree = ast.parse(source_code, filename=file_path)

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            func_name = node.name
            docstring = ast.get_docstring(node) or "No docstring available."
            try:
                func_source = ast.get_source_segment(source_code, node)
            except Exception:
                # fallback
                start_line = node.lineno - 1
                end_line = getattr(node, 'end_lineno', start_line + 1)
                func_source = "\n".join(source_code.splitlines()[start_line:end_line])
            print()
            yield {
                "name": func_name,
                "doc": docstring,
                "file": file_path,
                "source": func_source
            }

#@timed_node("extract_from_folder")
def extract_from_folder(folder_path):
    """
    Recursively traverses the specified folder and extracts all functions and their docstrings from Python (.py) files.

    Args:
        folder_path (str): The path of the target folder. The function will recursively traverse all subfolders within it.

    Returns:
        list: A list of dictionaries, each representing a function with the following keys:
            - ``name`` (str): The function name.
            - ``doc`` (str): The function's docstring. If no docstring is available, returns "No docstring available."
            - ``file`` (str): The path of the Python file containing the function.
    """
    all_functions = []
    for root, _, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".py"):  # 仅处理 Python 文件
                file_path = os.path.join(root, file)
                all_functions.extend(extract_single_function(file_path))
    return all_functions

#@timed_node("pycode_creator")
def pycode_creator(data_type, data_path, db_name, embedding_model=None):
    if data_type not in {"file", "folder"}:
        raise ValueError(f"Invalid data_type: {data_type}，data_type must be 'file' or 'folder'.")

    if not os.path.exists(data_path):
        raise FileNotFoundError(f"The specified path does not exist: {data_path}")

    if data_type == "file":
        code_functions = extract_functions_from_py(data_path)
    elif data_type == "folder":
        code_functions = extract_from_folder(data_path)

    documents = [
        Document(
            page_content=f"""Function: {func['name']}
Docstring: {func['doc']}
Full Code:
{func['source']}""",
            metadata={"name": func["name"], "file": func["file"]}
        )
        for func in code_functions
    ]
    #print(documents)
    db = FAISS.from_documents(documents, embedding_model)
    db.save_local(db_name)
    
    return db

def ask_vectorization_choice() -> bool:
    """询问用户是否需要向量化操作"""
    while True:
        choice = input("是否需要对新文件进行向量化操作？(yes/no): ").lower()
        if choice in ["yes", "y"]:
            return True
        elif choice in ["no", "n"]:
            return False
        else:
            print("请输入yes或no")

def get_file_or_folder_path() -> str:
    """获取用户输入的文件或文件夹路径"""
    while True:
        path = input("请输入要向量化的文件或文件夹路径: ")
        path = path.strip('"\'') 
        if os.path.exists(path):
            return path
        else:
            print("路径不存在，请重新输入")

def get_vectorstore_path() -> str:
    """获取用户输入的向量库路径"""
    while True:
        path = input("请输入要加载的向量库路径: ")
        path = path.strip('"\'')
        if os.path.exists(path):
            return path
        else:
            print("路径不存在，请重新输入")

class RAGState(BaseModel):
    query: str = ""
    vectorstore: Optional[Union[list[str],str,dict]] = None
    vect_name: Union[list[str],str,dict] = {}
    code_output: Union[list[str],str,dict] = ""
    notebook_output: Union[list[str],str,dict] = ""
    status: str = "START"
    search_k: int = 3
    score_threshold: float = 0.6
    run_msg: list = ["START"]
    final_result: str = ""
    next: Literal["rag_supervisor","load_vectorstore","search_code","search_notebook","generate_final_answer"] = "rag_supervisor"

#@timed_node("rag_supervisor")
def rag_supervisor(state: RAGState):
    template = f""" 
    #Role# 
    You are a supervisor of RAG (retrieval-augmented generation) system.

    #Task# 
    Base on the question intention and provided information, choose appropriate tool to do RAG.

    #Question# 
    query: {state.query}

    #Content#
    message: {state.run_msg[-1]}

    #Tool Description#
    "load_vectorstore": it can load vectorstore.
    "search_code": it can search python code information.
    "search_notebook": it can search code execution examples.
    "generate_final_answer": it can generate final answer.

    #Instruction#
    1. Check message to know your are in which step: the whole plan is "load_vectorstore" -> "search_code" -> "search_notebook" -> "generate_final_answer".
    2. If one step in the plan is skipped, go to next step.
    3. Do not repeat the steps.

    #Format#
    Return your answer with a valid JSON object in the following format: {{"thought": "your thinking process", "next_agent": "the selected tool name"}}. 
    Do not include any other text in your response, only the JSON object.  
    """
    chain = model_v3 | output_parser
    result = chain.invoke(template)
    print("rag_supervisor:", result["thought"])
    state.run_msg.append(result["thought"])
    return {
        "run_msg": state.run_msg,
        "next": result["next_agent"]
    }

#@timed_node("rag_supervisor")
def rag_supervisor_func(state: RAGState):
    template = f""" 
    #Role# 
    You are a supervisor of RAG (retrieval-augmented generation) system.

    #Task# 
    Base on the question intention and provided information, choose appropriate tool to do RAG.

    #Question# 
    query: {state.query}

    #Content#
    message: {state.run_msg[-1]}

    #Tool Description#
    "load_vectorstore": it can load vectorstore.
    "search_code": it can search python code information.
    "generate_final_answer": it can generate final answer.

    #Instruction#
    1. Check message to know your are in which step: the whole plan is "load_vectorstore" -> "search_code"  -> "generate_final_answer".
    2. If one step in the plan is skipped, go to next step.
    3. Do not repeat the steps.

    #Format#
    Return your answer with a valid JSON object in the following format: {{"thought": "your thinking process", "next_agent": "the selected tool name"}}. 
    Do not include any other text in your response, only the JSON object.  
    """
    chain = model_v3 | output_parser
    result = chain.invoke(template)
    print("rag_supervisor:", result["thought"])
    state.run_msg.append(result["thought"])
    return {
        "run_msg": state.run_msg,
        "next": result["next_agent"]
    }

#@timed_node("load_vectorstore")
def load_vectorstore(state: RAGState):
    # 如果用户未提供向量库路径，则询问并处理
    if not state.vectorstore:
        print("Vectorstore not provided. Please choose an option:")
        print("1. Vectorize new files")
        print("2. Use existing vectorstore")
        
        while True:    
            choice = input("Enter your choice (1 or 2): ")
            if choice == "1":
                print("Choose your content type:")
                print("a. Code")
                print("b. Notebook")
                print("c. Both")

                content_type = input("Enter your choice (a or b or c): ").lower()
                if content_type == "a":
                    path = get_file_or_folder_path()
                    print(f"Now saving code in: {path}...")
                    db_c = pycode_creator("file" if os.path.isfile(path) else "folder", path, "code_vectorstore",embedding_model=embedding_model)
                    state.vectorstore = "code_vectorstore"
                    break
                elif content_type == "b":
                    print(f"Now saving notebook in: {path}...")
                    db_n = notebook_creator_cell(path, include_outputs=False, max_output_length=20, remove_newline=True, db_name="notebook_vectorstore",embedding_model=embedding_model)
                    state.vectorstore = "notebook_vectorstore"
                    break
                elif content_type == "c":
                    c_path = get_file_or_folder_path()
                    print(f"Now saving code in: {c_path}...")
                    db_c = pycode_creator("file" if os.path.isfile(path) else "folder", path, "code_vectorstore",embedding_model=embedding_model)
                    n_path = get_file_or_folder_path()
                    print(f"Now saving notebook in: {n_path}...")
                    db_n = notebook_creator_cell(n_path, include_outputs=False, max_output_length=20, remove_newline=True, db_name="notebook_vectorstore",embedding_model=embedding_model)
                    state.vectorstore = ["code_vectorstore","notebook_vectorstore"]
                    break
                else:
                    print("无效选择，请重新运行程序")
                    continue

            elif choice == "2":
                state.vectorstore = get_vectorstore_path()
                break
            else:
                print("Invalid choice. Please enter 1 or 2.")
    
    if isinstance(state.vectorstore, str):
        if not state.vectorstore:
            raise ValueError("Provided vectorstore path cannot be empty.")
        if not os.path.exists(state.vectorstore):
            raise FileNotFoundError(f"Vectorstore path does not exist: {state.vectorstore}")
        
        name = Path(state.vectorstore).stem  # Extract filename without extension
        vectorstore = FAISS.load_local(state.vectorstore, embeddings=embedding_model, allow_dangerous_deserialization=True)
        state.vect_name = {name: vectorstore}
        print("Load vectorstore:",state.vect_name.keys())
        state.run_msg.append("Successfully loaded vectorstore.")
        return {
            "run_msg": state.run_msg,
            "vect_name": state.vect_name
        }

    elif isinstance(state.vectorstore, list):
        if not state.vectorstore:
            raise ValueError("提供的vectorstore列表不能为空。")
        
        state.vect_name = {}
        for path in state.vectorstore:
            if not isinstance(path, str):
                raise TypeError("列表中的元素必须为字符串类型。")
            if not path:
                raise ValueError("列表中的路径不能为空字符串。")
            if not os.path.exists(path):
                raise FileNotFoundError(f"Vectorstore path does not exist: {path}")
            
            name = Path(path).stem
            vectorstore = FAISS.load_local(path, embeddings=embedding_model, allow_dangerous_deserialization=True)
            state.vect_name[name] = vectorstore
        print("Load vectorstore:",state.vect_name.keys())
        state.run_msg.append("Successfully loaded vectorstore.")
        return {
            "run_msg": state.run_msg,
            "vect_name": state.vect_name
        }

    elif isinstance(state.vectorstore, dict):
        if not state.vectorstore:
            raise ValueError("提供的vectorstore字典不能为空。")
        
        state.vect_name = {}
        for name, path in state.vectorstore.items():
            if not isinstance(path, str):
                raise TypeError(f"{name} 的路径必须是字符串类型。")
            if not path:
                raise ValueError(f"{name} 的路径不能为空字符串。")
            if not os.path.exists(path):
                raise FileNotFoundError(f"Vectorstore path does not exist: {path}")

            vectorstore = FAISS.load_local(path, embeddings=embedding_model, allow_dangerous_deserialization=True)
            state.vect_name[name] = vectorstore
        print("Load vectorstore:",state.vect_name.keys())
        state.run_msg.append("Successfully loaded vectorstore.")
        return {
            "run_msg": state.run_msg,
            "vect_name": state.vect_name
        }
    
    else:
        raise TypeError("vectorstore必须是字符串、字符串列表或字典。")

def fuzzy_ratio_match(keyword, text, threshold=0):
    return fuzz.token_sort_ratio(keyword.lower(), text.lower()) >= threshold

#@timed_node("search_code")
def search_code(state: RAGState):
    vect_temp = f""" 
    #Task# 
    Follow the instruction to extract the name of the vectorstore that stores the code information from all the vectorstore.

    #Content#
    {state.vect_name.keys()}

    #Instruction# 
    Determine whether the name of the vectorstore is related to "code". If it is related, extract it; otherwise, do not extract it.

    #Example# 
    1)
    vectorstore = ["BRICK_code","python","jupyter_notebook","markdown","txt","py","notebook"]
    names = ["BRICK_code", "python", "py"] 
    2)
    vectorstore = ["jupyter_notebook","plot","format","png"]
    names = [] 

    #Output#
    Return the following JSON object: {{"names": "use list object to only store the name of the vectorstore that stores the code information "}}
    Do not include any other text in your response, only the JSON object.
    """
    chain = model_v3 | output_parser
    names = chain.invoke(vect_temp)

    # 检索+阈值过滤+排序
    result = []
    vs = names.get("names", [])
    if len(vs) > 0:
        key_temp = f""" 
        #Task# 
        Extract the keywords from the problem, such as function names, package names, and purposes.

        #Question#
        {state.query}

        #Example# 
        query = "Integration by BRICK.pp.complete_results() and Visualization by BRICK.pl.visualization()"
        keywords = ["Integration", "BRICK.pp.complete_results()", "BRICK", "BRICK.pp", "complete_results()", "complete_results", "Visualization", "BRICK.pl.visualization()", "BRICK.pl", "visualization()", "visualization"] 

        #Output#
        Return the following JSON object: {{"keywords": "use list object to only store the keywords in the query"}}
        Do not include any other text in your response, only the JSON object.
        """
        chain = model_v3 | output_parser
        keywords = chain.invoke(key_temp)
        keyword = keywords.get("keywords", [])
        print("Query keywords are",keyword)
        for name in vs:
            print("n: ",name)
            results = state.vect_name[name].similarity_search_with_score(state.query, k=state.search_k)
            filtered_sorted_results = sorted(
                [(doc, score) for doc, score in results if score >= state.score_threshold],
                key=lambda x: x[1],
                reverse=True
            )
            for i, (doc, score) in enumerate(filtered_sorted_results):
                print(f"[{i+1}] Score: {score:.4f} | Content: {doc.page_content}")
            #print("filtered_sorted_results",filtered_sorted_results)
            
            """ keyword_filtered_results = []
            for doc, score in filtered_sorted_results:
                print(f"\n文档预览: {doc.page_content[:150]}...")
                print(f"Score: {score}")
                matched = False  # 用来标记是否匹配到关键词

                print("Query keywords are",keyword)
                keyword = str(keyword)
                found_regex = re.search(re.escape(keyword), doc.page_content, re.IGNORECASE)
                fuzzy_score = fuzz.partial_ratio(keyword.lower(), doc.page_content.lower())

                print(f"检查关键词: {keyword} | 正则匹配: {bool(found_regex)} | 模糊分数: {fuzzy_score}")
                
                if found_regex or fuzzy_score >= 45:
                    matched = True
                    break  

                if matched:
                    keyword_filtered_results.append((doc, score)) """
            
            keyword_score_results = []
            for doc, _ in filtered_sorted_results:  
                max_score = max(
                    fuzz.partial_ratio(keyword.lower(), doc.page_content.lower())
                    for keyword in keywords
                )
                keyword_score_results.append((doc, max_score))

            #top_k_results = sorted(keyword_score_results, key=lambda x: x[1], reverse=True)[:1]

            if len(keyword_score_results) == 0:
                result.append([doc.page_content for doc, _ in filtered_sorted_results])
            else:
                max_score = max(score for _, score in keyword_score_results)
                top_k_results = [(doc, score) for doc, score in keyword_score_results if score == max_score]

                print("\n 根据关键词模糊匹配得分排序结果:")
                for i, (doc, fuzzy_score) in enumerate(top_k_results):
                    print(f"[{i+1}] Fuzzy Score: {fuzzy_score:.2f} | Content: {doc.page_content[:200]}...")

                result.append([doc.page_content for doc, _ in top_k_results])

            #print("keyword_filtered_results",keyword_filtered_results)
            #result.append([doc.page_content for doc, _ in keyword_filtered_results])
            #print("result",result)
        template = f"""
        #Task#
        Only based on the code context, select the approriate answer to generate a valid python code for the question.
        
        #Context#
        {result}
        
        #Question#
        {state.query}

        #Instruction#
        1. Use the most appropriate code retrieved for the response, and the output code must be exactly the same as the most suitable code. 
        2. The specific function content must be included, and only the function name cannot be answered.
        
        #Output#
        Return your answer in a valid str object.
        """
        chain = model_v3 | str_output_parser
        final_result = chain.invoke(template)
        print("find function:",final_result)
        state.run_msg.append("Successfully search code.")
        print("search_code:",state.run_msg[-1])
    else:
        final_result = ""
        state.run_msg.append("There is no vector library name related to the code, so skipping search code.")
        print("search_code:",state.run_msg[-1])
    return {
        "run_msg": state.run_msg,
        "code_output": final_result
    }

#@timed_node("search_notebook")
def search_notebook(state: RAGState):
    vect_temp = f""" 
    #Task# 
    Follow the instruction to extract the name of the vectorstore that stores the notebook information from all the vectorstore.

    #Content#
    {state.vect_name.keys()}

    #Instruction# 
    Determine whether the name of the vectorstore is related to "notebook". If it is related, extract it; otherwise, do not extract it.

    #Example# 
    1)
    vectorstore = ["BRICK_code","python","jupyter_notebook","markdown","txt","py","notebook"]
    names = ["jupyter_notebook", "markdown", "txt", "notebook"] 
    2)
    vectorstore = ["BRICK","BRICK_code","BRICK_code2","code","plot","format","png"]
    names = [] 

    #Output#
    Return the following JSON object: {{"names": "use list object to only store the name of the vectorstore that stores the notebook information"}}
    Do not include any other text in your response, only the JSON object.
    """
    chain = model_v3 | output_parser
    names = chain.invoke(vect_temp)

    # 检索+阈值过滤+排序
    result = []
    vs = names.get("names", [])
    if len(vs) > 0:
        key_temp = f""" 
        #Task# 
        Extract the keywords from the problem, such as function names, package names, and purposes.

        #Question#
        {state.query}

        #Example# 
        query = "Integration by BRICK.pp.complete_results() and Visualization by BRICK.pl.visualization()"
        keywords = ["Integration", "BRICK.pp.complete_results()", "BRICK", "BRICK.pp", "complete_results()", "complete_results", "Visualization", "BRICK.pl.visualization()", "BRICK.pl", "visualization()", "visualization"] 

        #Output#
        Return the following JSON object: {{"keywords": "use list object to only store the keywords in the query"}}
        Do not include any other text in your response, only the JSON object.
        """
        chain = model_v3 | output_parser
        keywords = chain.invoke(key_temp)
        keyword = keywords.get("keywords", [])
        print("Query keywords are",keyword)
        for name in vs:
            print("n: ",name)
            results = state.vect_name[name].similarity_search_with_score(state.query, k=state.search_k)
            filtered_sorted_results = sorted(
                [(doc, score) for doc, score in results if score >= state.score_threshold],
                key=lambda x: x[1],
                reverse=True
            )
            for i, (doc, score) in enumerate(filtered_sorted_results):
                print(f"[{i+1}] Score: {score:.4f} | Content: {doc.page_content[:200]}...")

            keyword_score_results = []
            for doc, _ in filtered_sorted_results:  
                max_score = max(
                    fuzz.partial_ratio(keyword.lower(), doc.page_content.lower())
                    for keyword in keywords
                )
                keyword_score_results.append((doc, max_score))

            if len(keyword_score_results) == 0:
                result.append([doc.page_content for doc, _ in filtered_sorted_results])
            else:
                max_score = max(score for _, score in keyword_score_results)
                top_k_results = [(doc, score) for doc, score in keyword_score_results if score == max_score]

                print("\n 根据关键词模糊匹配得分排序结果:")
                for i, (doc, fuzzy_score) in enumerate(top_k_results):
                    print(f"[{i+1}] Fuzzy Score: {fuzzy_score:.2f} | Content: {doc.page_content[:200]}...")

                result.append([doc.page_content for doc, _ in top_k_results])

        template = f"""
            #Task#
            Only based on the notebook context, select the approriate answer to generate a valid answer for the question.
            
            #Context#
            {result}
            
            #Question#
            {state.query}
            
            #Instruction# 
            1. Use the most appropriate context retrieved for the response, and the output context must be exactly the same as the most suitable context. 
            2. The output context needs to be complete. 

            #Output#
            Return your answer in a valid str object.
        """
        chain = model_v3 | str_output_parser
        final_result = chain.invoke(template)
        print("find notebook:",final_result)
        state.run_msg.append("Successfully search notebook.")
        print("search_notebook:",state.run_msg[-1])
    else:
        final_result = ""
        state.run_msg.append("There is no vector library name related to the notebook, so skipping search notebook.")
        print("search_notebook:",state.run_msg[-1])
    return {
        "run_msg": state.run_msg,
        "notebook_output": final_result
    }

def generate_final_answer(state: RAGState):
    template = f"""
    #Content#
    - query:{state.query}
    - code information:{state.code_output}
    - code execution examples:{state.notebook_output}
    
    #Instruction#
    1. Check the code information and the code execution examples are relevant or not.
        If relevant, complete the code context by using these two information, but avoid using repetitive parts.
            For example: you found 
                ```python
                import scanpy as sc
                adata = sc.read_h5ad(file_path)
                ``` 
                and 
                ```python
                def load_preprocessed_h5ad(file_path):
                adata = sc.read_h5ad(file_path)
                return adata
                ```
                The second piece of code employs the same code as the first one: adata = sc.read_h5ad(file_path). The difference between these two pieces of code lies in that the second one uses another load_preprocessed_h5ad function to wrap the sc.read_h5ad function.
                Thus, you can use the first code with real file_path.
        Else, check the query and select the best and relevant code to answer this query.
    2. Not only provide the function name, but also provide the function content.

    #Example#
    1) The code information provides a function `run_paga` that encapsulates the execution of `sc.tl.paga`, while the code execution example directly shows how to run `sc.tl.paga` on an `adata` object.
    To complete the code context without repetition, we can use the function `sc.tl.paga` to run the PAGA.
    
    2) The code information provides a function `sc.read_h5ad` that can load h5ad file, while the code execution example shows a function `load_preprocessed_h5ad` that encapsulates the execution of `sc.read_h5ad`
    To complete the code context without repetition, we can use the function `sc.read_h5ad` to load h5ad.

    #Format#
    Return your final answer in a valid str object.
    """
    chain = model_q | str_output_parser
    final_result = chain.invoke(template)
    print("final_result: ",final_result)
    return {"final_result": final_result}

def run_route(state: RAGState) -> str:
    if state.status == "FINISHED":
        return "end"
    else:
        if state.next == "rag_supervisor":
            return "rag_supervisor"
        elif state.next == "load_vectorstore":
            return "load_vectorstore"
        elif state.next == "search_code":
            return "search_code"
        elif state.next == "search_notebook":
            return "search_notebook"
        elif state.next == "generate_final_answer":
            return "generate_final_answer"
        return "Your next agent is not defined"

def run_route_func(state: RAGState) -> str:
    if state.status == "FINISHED":
        return "end"
    else:
        if state.next == "rag_supervisor":
            return "rag_supervisor"
        elif state.next == "load_vectorstore":
            return "load_vectorstore"
        elif state.next == "search_code":
            return "search_code"
        elif state.next == "generate_final_answer":
            return "generate_final_answer"
        return "Your next agent is not defined"

#@timed_node("RAG")
def RAG(query, vectorstore = None, search_k: int = 3, score_threshold: float = 0.6):
    print("query data:",query)
    state = {
        "query": query,
        "vectorstore": vectorstore,
        "search_k": search_k,
        "score_threshold": score_threshold
    }
    
    initial_state = RAGState(**state)

    builder = StateGraph(RAGState)
    builder.add_node("rag_supervisor", rag_supervisor)
    builder.add_node("load_vectorstore", load_vectorstore)
    builder.add_node("search_code", search_code)
    builder.add_node("search_notebook", search_notebook)
    builder.add_node("generate_final_answer", generate_final_answer)

    # Logic
    builder.add_edge(START, "rag_supervisor")
    builder.add_conditional_edges("rag_supervisor", run_route, {"load_vectorstore":"load_vectorstore", "search_code":"search_code", "search_notebook":"search_notebook", "generate_final_answer":"generate_final_answer"})
    builder.add_edge("load_vectorstore", "rag_supervisor")
    builder.add_edge("search_code", "rag_supervisor")
    builder.add_edge("search_notebook", "rag_supervisor")
    builder.add_edge("generate_final_answer", END)

    # Add
    graph = builder.compile()

    # View
    # display(Image(graph.get_graph().draw_mermaid_png()))
    config = {"configurable": {"thread_id": "1"},"recursion_limit": 200}
    answer = graph.invoke(initial_state, config)
    return answer["final_result"]

def RAG_func(query, vectorstore = None, search_k: int = 3, score_threshold: float = 0.6):
    print("query data:",query)
    
    # 检查query是否包含"BRICK"（不区分大小写）
    if "brick" not in query.lower():
        print("query not contain BRICK")
        return ""
    
    state = {
        "query": query,
        "vectorstore": vectorstore,
        "search_k": search_k,
        "score_threshold": score_threshold
    }
    
    initial_state = RAGState(**state)

    builder = StateGraph(RAGState)
    builder.add_node("rag_supervisor_func", rag_supervisor_func)
    builder.add_node("load_vectorstore", load_vectorstore)
    builder.add_node("search_code", search_code)
    builder.add_node("generate_final_answer", generate_final_answer)

    # Logic
    builder.add_edge(START, "rag_supervisor_func")
    builder.add_conditional_edges("rag_supervisor_func", run_route, {"load_vectorstore":"load_vectorstore", "search_code":"search_code",  "generate_final_answer":"generate_final_answer"})
    builder.add_edge("load_vectorstore", "rag_supervisor_func")
    builder.add_edge("search_code", "rag_supervisor_func")
    builder.add_edge("generate_final_answer", END)

    # Add
    graph = builder.compile()

    # View
    # display(Image(graph.get_graph().draw_mermaid_png()))
    config = {"configurable": {"thread_id": "1"},"recursion_limit": 200}
    answer = graph.invoke(initial_state, config)
    return answer["final_result"]

if __name__ == '__main__':
    a = RAG(query="Use BRICK.pp.complete_results to integrate KG-based insights with trajectory data.",vectorstore={"Code":"/home/lyt/checker_finallap/BRICK2/BRICK-main/vectorstore/BRICK_code4.faiss","Notebook":"/home/lyt/checker_finallap/BRICK/vectorstore/BRICK_notebook2.faiss"})
    print(a)
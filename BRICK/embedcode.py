import os
import ast
import glob
import json

#from langchain.document_loaders import NotebookLoader
from langchain_community.document_loaders import NotebookLoader
#from langchain.vectorstores import FAISS
from langchain_community.vectorstores import FAISS
from langchain_core.documents import Document

from ._settings import get_llm_embedding

def notebook_creator2(file_path, include_outputs=False, max_output_length=20, remove_newline=True, db_name="BRICK_notebook", embedding_model=None):
    """
    Vectorizes Jupyter Notebook (.ipynb) documents and stores them in a FAISS vector database.

    Args:
        file_path (str): The path to the Notebook file (supports wildcards, e.g., `./notebooks/*.ipynb`).
        include_outputs (bool): Whether to include the output of code cells, default is False.
        max_output_length (int): The maximum length of code cell outputs; outputs beyond this length will be truncated.
        remove_newline (bool): Whether to remove newline characters in the text to reduce redundancy during vectorization.
        db_name (str): The name of the vector database index file.
        embedding_model (any, optional): The embedding model used to convert the document into vector representations. Defaults to None (uses the initialized embedding model).

    Returns:
        FAISS: The FAISS database instance storing the vectorized Notebook documents.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The specified path does not exist: {file_path}")
    
    if embedding_model is None:
        embedding_model = get_llm_embedding()
    else:
        embedding_model = embedding_model

    batch_size = 10 
    all_text_embedding_pairs = []  # 用于收集所有的文本和嵌入对
    all_metadatas = []  # 用于收集所有的元数据
    documents = []
    files = glob.glob(file_path)
    for file in files:
        loader = NotebookLoader(file, include_outputs, max_output_length, remove_newline,)
        data = loader.load() 
        for d in data:
            documents.append(d)
    if type(embedding_model).__name__ == "OpenAI":
        texts = [doc.page_content for doc in documents]
        for i in range(0, len(texts), batch_size):
            batch = texts[i:i + batch_size]
            response = embedding_model.embeddings.create(input=batch, model=embedding_model.embeddings.model)
            embeddings = [item.embedding for item in response.data]
            text_embedding_pairs = list(zip(batch, embeddings))
            all_text_embedding_pairs.extend(text_embedding_pairs)
            all_metadatas.extend([doc.metadata for doc in documents[i:i + batch_size]])
        db = FAISS.from_embeddings(
            all_text_embedding_pairs,
            embedding_model,
            metadatas=all_metadatas
        )
    else:
        db = FAISS.from_documents(documents, embedding_model)
    db.save_local(db_name)
    return db

def notebook_creator(file_path, include_outputs=False, max_output_length=20, remove_newline=True, db_name="BRICK_notebook", embedding_model=None):
    """
    Vectorizes Jupyter Notebook (.ipynb) documents and stores them in a FAISS vector database.

    Args:
        file_path (str): The path to the Notebook files (supports wildcards, e.g., `./notebooks/*.ipynb`).
        include_outputs (bool): Whether to include the output of code cells, default is False.
        max_output_length (int): The maximum length of code cell outputs; outputs beyond this length will be truncated.
        remove_newline (bool): Whether to remove newline characters in the text to reduce redundancy during vectorization.
        db_name (str): The name of the vector database index file.
        embedding_model (any, optional): The embedding model used to convert the document into vector representations. Defaults to None (uses the initialized embedding model).

    Returns:
        FAISS: The FAISS database instance storing the vectorized Notebook documents.
    """
    if os.path.isdir(file_path):
        file_path = os.path.join(file_path, "*.ipynb")

    files = glob.glob(file_path)
    if not files:
        raise FileNotFoundError(f"No .ipynb files found in the specified path: {file_path}")
    
    if embedding_model is None:
        embedding_model = get_llm_embedding()

    documents = []
    files = glob.glob(file_path)
    for file in files:
        loader = NotebookLoader(file, include_outputs, max_output_length, remove_newline)
        data = loader.load() 
        documents.extend(data) 
    
    db = FAISS.from_documents(documents, embedding_model)
    db.save_local(db_name)
    return db

def notebook_creator_cell(file_path, include_outputs=False, max_output_length=20, remove_newline=True, db_name="BRICK_notebook", embedding_model=None):
    if os.path.isdir(file_path):
        file_path = os.path.join(file_path, "*.ipynb")

    files = glob.glob(file_path)
    if not files:
        raise FileNotFoundError(f"No .ipynb files found in the specified path: {file_path}")
    
    if embedding_model is None:
        embedding_model = get_llm_embedding()

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
                #all_functions.extend(extract_functions_from_py(file_path))
                all_functions.extend(extract_single_function(file_path))
    return all_functions

def pycode_creator2(data_type, data_path, db_name, embedding_model=None):
    """
    Parses Python code (either a single file or an entire folder), extracts functions and their docstrings, 
    and stores the results in a FAISS vector database.

    Args:
        data_type (str): The type of data to process. Possible values:
            - ``"file"``: A single Python file (.py).
            - ``"folder"``: A folder containing multiple Python files.
        data_path (str): The path to the Python file or folder.
        db_name (str): The path where the FAISS database will be stored.
        embedding_model (Any, optional): The embedding model used to convert extracted code snippets into vectors. 
            If not provided, the default initialized embedding model will be used.

    Returns:
        FAISS: A FAISS database instance that stores the vectorized Python code.
    """
    if data_type not in {"file", "folder"}:
        raise ValueError(f"Invalid data_type: {data_type}，data_type must be 'file' or 'folder'.")

    if not os.path.exists(data_path):
        raise FileNotFoundError(f"The specified path does not exist: {data_path}")

    if embedding_model is None:
        embedding_model = get_llm_embedding()
    else:
        embedding_model = embedding_model

    if data_type == "file":
        code_functions = extract_functions_from_py(data_path)
    elif data_type == "folder":
        code_functions = extract_from_folder(data_path)

    documents = [
        Document(
            page_content=f"Function: {func['name']}\nDocstring: {func['doc']}",
            metadata={"name": func["name"], "file": func["file"]}
        )
        for func in code_functions
    ]

    batch_size = 10 
    all_text_embedding_pairs = [] 
    all_metadatas = [] 

    # 如果直接使用OpenAI官方的model，则需要使用embeddings.create方法（只接受字符串）
    if type(embedding_model).__name__ == "OpenAI":
        texts = [doc.page_content for doc in documents]
        for i in range(0, len(texts), batch_size):
            batch = texts[i:i + batch_size] 
            response = embedding_model.embeddings.create(input=batch, model=embedding_model.embeddings.model) 
            embeddings = [item.embedding for item in response.data]  
            text_embedding_pairs = list(zip(batch, embeddings))  
            
            all_text_embedding_pairs.extend(text_embedding_pairs)
            all_metadatas.extend([doc.metadata for doc in documents[i:i + batch_size]])
        db = FAISS.from_embeddings(
            all_text_embedding_pairs,
            embedding_model.embeddings.model,
            metadatas=all_metadatas
        )   
    else:
        db = FAISS.from_documents(documents, embedding_model)
    db.save_local(db_name)
    
    return db

def pycode_creator(data_type, data_path, db_name, embedding_model=None):
    if data_type not in {"file", "folder"}:
        raise ValueError(f"Invalid data_type: {data_type}，data_type must be 'file' or 'folder'.")

    if not os.path.exists(data_path):
        raise FileNotFoundError(f"The specified path does not exist: {data_path}")

    if embedding_model is None:
        embedding_model = get_llm_embedding()
    else:
        embedding_model = embedding_model

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


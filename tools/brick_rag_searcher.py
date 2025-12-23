import os
from pathlib import Path
from langchain_community.vectorstores import FAISS
from typing import Literal, Union, List, Dict, Optional
from langgraph.graph import StateGraph, START, END
from pydantic import BaseModel
from IPython.display import Image, display
from rapidfuzz import fuzz
from langchain.tools import tool
from llm import embedding_model, basic_llm as model

class RAGState(BaseModel):
    query: str = ""
    vectorstore: Optional[Union[list[str],str,dict]] = None
    vect_name: Union[list[str],str,dict] = ""
    code_output: Union[list[str],str,dict] = ""
    notebook_output: Union[list[str],str,dict] = ""
    status: str = "START"
    search_k: int = 5
    score_threshold: float = 0.6
    run_msg: list = ["START"]
    final_result: str = ""
    next: Literal["load_vectorstore","search_code","search_notebook","generate_final_answer"] = "load_vectorstore"

def fuzzy_ratio_match(keyword, text, threshold=0):
    return fuzz.token_sort_ratio(keyword.lower(), text.lower()) >= threshold

def load_vectorstore(state: RAGState):
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
    chain = model | output_parser
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
        chain = model | output_parser
        keywords = chain.invoke(key_temp)
        keyword = keywords.get("keywords", [])
        #print("Query keywords are",keyword)
        for name in vs:
            #print("n: ",name)
            results = state.vect_name[name].similarity_search_with_score(state.query, k=state.search_k)
            filtered_sorted_results = sorted(
                [(doc, score) for doc, score in results if score >= state.score_threshold],
                key=lambda x: x[1],
                reverse=True
            )
            #for i, (doc, score) in enumerate(filtered_sorted_results):
                #print(f"[{i+1}] Score: {score:.4f} | Content: {doc.page_content}")
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

                #print("\n 根据关键词模糊匹配得分排序结果:")
                #for i, (doc, fuzzy_score) in enumerate(top_k_results):
                    #print(f"[{i+1}] Fuzzy Score: {fuzzy_score:.2f} | Content: {doc.page_content[:200]}...")

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
        chain = model | str_output_parser
        final_result = chain.invoke(template)
        #print("find function:",final_result)
        state.run_msg.append("Successfully search code.")
        #print("search_code:",state.run_msg[-1])
    else:
        final_result = ""
        state.run_msg.append("There is no vector library name related to the code, so skipping search code.")
        #print("search_code:",state.run_msg[-1])
    return {
        "run_msg": state.run_msg,
        "code_output": final_result
    }

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
    chain = model | output_parser
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
        chain = model | output_parser
        keywords = chain.invoke(key_temp)
        keyword = keywords.get("keywords", [])
        #print("Query keywords are",keyword)
        for name in vs:
            #print("n: ",name)
            results = state.vect_name[name].similarity_search_with_score(state.query, k=state.search_k)
            filtered_sorted_results = sorted(
                [(doc, score) for doc, score in results if score >= state.score_threshold],
                key=lambda x: x[1],
                reverse=True
            )
            #for i, (doc, score) in enumerate(filtered_sorted_results):
                #print(f"[{i+1}] Score: {score:.4f} | Content: {doc.page_content[:200]}...")

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

                #print("\n 根据关键词模糊匹配得分排序结果:")
                #for i, (doc, fuzzy_score) in enumerate(top_k_results):
                    #print(f"[{i+1}] Fuzzy Score: {fuzzy_score:.2f} | Content: {doc.page_content[:200]}...")

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
        chain = model | str_output_parser
        final_result = chain.invoke(template)
        #print("find notebook:",final_result)
        state.run_msg.append("Successfully search notebook.")
        #print("search_notebook:",state.run_msg[-1])
    else:
        final_result = ""
        state.run_msg.append("There is no vector library name related to the notebook, so skipping search notebook.")
        #print("search_notebook:",state.run_msg[-1])
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
    #print("final_result: ",final_result)
    return {"final_result": final_result}
    
def run_route(state: RAGState) -> str:
    routing_map = {
        "load_vectorstore": "load_vectorstore",
        "search_code": "search_code",
        "search_notebook": "search_notebook",
        "generate_final_answer": "generate_final_answer",
    }
    if state.status == "FINISHED":
        return "end"
    elif state.next in routing_map:
        return routing_map[state.next]
    else:
        raise ValueError(f"Invalid next step: {state.next}")

def RAG(query, vectorstore = None, search_k: int = 5, score_threshold: float = 0.6):
    state = {
        "query": query,
        "vectorstore": vectorstore,
        "search_k": search_k,
        "score_threshold": score_threshold
    }
    
    initial_state = RAGState(**state)

    builder = StateGraph(RAGState)
    builder.add_node("load_vectorstore", load_vectorstore)
    builder.add_node("search_code", search_code)
    builder.add_node("search_notebook", search_notebook)
    builder.add_node("generate_final_answer", generate_final_answer)

    # Logic
    builder.add_edge(START, "load_vectorstore")
    builder.add_edge("load_vectorstore", "search_code")
    builder.add_edge("search_code", "search_notebook")
    builder.add_edge("search_notebook", "generate_final_answer")
    builder.add_edge("generate_final_answer", END)

    # Add
    graph = builder.compile()

    # View
    # display(Image(graph.get_graph().draw_mermaid_png()))
    config = {"configurable": {"thread_id": "1"},"recursion_limit": 200}
    answer = graph.invoke(initial_state, config)
    return answer["final_result"]


@tool
def perform_rag_search_tool(
    query: str,
    vectorstore_paths: Union[str, List[str], Dict[str, str]], 
    search_k: int = 5,
    score_threshold: float = 0.6
) -> str:
    """
    Performs a Retrieval Augmented Generation (RAG) search for code and notebook examples
    based on a user query and returns a final answer. This tool loads vectorstores from provided paths
    before performing the search.

    Args:
        query (str): The user's question or problem statement.
        vectorstore_paths (Union[str, List[str], Dict[str, str]]):
            The path(s) or dictionary mapping names to paths of the local FAISS vectorstores.
            - If a single string, it's treated as a path to one vectorstore.
            - If a list of strings, each string is a path to a vectorstore.
            - If a dictionary, keys are logical names (e.g., "code_db", "notebook_docs")
              and values are the file paths to the corresponding local FAISS vectorstores.
              These vectorstores will be loaded internally by the RAG pipeline.
        search_k (int, optional): The number of top results to retrieve from each vectorstore. Defaults to 5.
        score_threshold (float, optional): The minimum similarity score for retrieved documents to be considered. Defaults to 0.6.

    Returns:
        str: The final answer, which is a valid Python code snippet or explanation
             generated by combining retrieved code and notebook information.
    """
    # 直接将 vectorstore_paths 传递给 RAG 函数，RAG 内部会调用 load_vectorstore 处理
    return RAG(query=query, vectorstore=vectorstore_paths, search_k=search_k, score_threshold=score_threshold)

@tool
def run_rag_pipeline_tool(
    query: str,
    vect_name: dict,
    search_k: int = 5,
    score_threshold: float = 0.6,
) -> dict:
    """
    Run RAG pipeline to answer a question using vectorstores that may contain code and notebook context.

    Args:
        query (str): The user's question or problem.
        vect_name (dict): A dictionary of named vectorstores, e.g. {"code": code_vs, "notebook": nb_vs}
        search_k (int): Number of top documents to retrieve per store.
        score_threshold (float): Threshold to filter similarity scores.

    Returns:
        dict: Final answer and intermediate outputs.
    """
    state = RAGState(
        query=query,
        vect_name=vect_name,
        search_k=search_k,
        score_threshold=score_threshold,
        run_msg=["START"]
    )
    
    # Step 1: Code Retrieval
    code_result = search_code(state)
    state.code_output = code_result["code_output"]
    state.run_msg = code_result["run_msg"]

    # Step 2: Notebook Retrieval
    notebook_result = search_notebook(state)
    state.notebook_output = notebook_result["notebook_output"]
    state.run_msg = notebook_result["run_msg"]

    # Step 3: Final Answer Generation
    answer_result = generate_final_answer(state)
    state.final_result = answer_result

    return {
        "final_result": state.final_result,
        "code_output": state.code_output,
        "notebook_output": state.notebook_output,
        "run_log": state.run_msg,
    }

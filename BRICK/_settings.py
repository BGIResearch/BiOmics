__version__:str = '0.0.1'
from neo4j import GraphDatabase
#from langchain.chat_models import AzureChatOpenAI
from langchain_community.chat_models import AzureChatOpenAI
from langchain_openai.chat_models import ChatOpenAI
#from langchain.embeddings.openai import OpenAIEmbeddings
from langchain_openai import OpenAIEmbeddings
from openai import OpenAI
from langchain_core.embeddings import Embeddings

# 全局 driver 实例
_driver_instance = None
_llm_instance = None
_llm_embedding_instance = None

def config(url, auth=None, encrypted=False):
    """
    Configures the database connection and initializes the global driver instance.

    Args:
        url (str): The database connection URL (e.g., "bolt://localhost:7687").
        auth (tuple, optional): A tuple containing authentication credentials (username, password). Defaults to None.
        encrypted (bool, optional): Whether to enable encrypted communication. Defaults to False.

    Returns:
        None
    """
    global _driver_instance
    _driver_instance = GraphDatabase.driver(url, auth=auth, encrypted=encrypted)
    print(f"Graph database has been configured and initialized successfully.")

    # if _driver_instance is None:
    #     _driver_instance = GraphDatabase.driver(url, auth=auth, encrypted=encrypted)
    # else:
    #     raise RuntimeError("Database has already been configured!")

def get_driver():
    """
    Retrieves the global driver instance.

    Returns:
        GraphDatabase.driver or None: The configured driver instance if available, otherwise None.
    """
    if _driver_instance is None:
        raise RuntimeError("Database has not been configured. Call BRICK.config(url) first.")
    return _driver_instance

def close():
    """
    Closes the global driver instance.
    """
    global _driver_instance
    if _driver_instance is not None:
        _driver_instance.close()
        _driver_instance = None

def config_llm(modeltype='AzureChatOpenAI', base_url=None, api_key=None, llm_params={'LLM':{}}):
    """
    Configures and initializes a globally unique LLM instance.

    Args:
        modeltype (str): The type of the model. Defaults to 'AzureChatOpenAI'.
        base_url (str, optional): The base URL of the service. Defaults to None.
        api_key (str, optional): The API key for the model. Defaults to None.
        llm_params (dict): Parameters for the model. Defaults to {'LLM': {}}.

    Returns:
        AzureChatOpenAI: An instance of the configured LLM.
    """
    global _llm_instance
    if modeltype == 'AzureChatOpenAI':
        try:
            """ if "model_name" in llm_params.keys():
                model_name = llm_params["model_name"] """
            _llm_instance = AzureChatOpenAI(base_url,api_key,**llm_params)
            print(f"LLM has been configured and initialized successfully.")
        except Exception as e:
            raise ValueError(f"Failed to initialize LLM: {e}")
    elif modeltype == 'ChatOpenAI':
        try:
            """ if "model_name" in llm_params.keys():
                model_name = llm_params["model_name"] """
            # print(llm_params[first_key]["llm_params"])
            _llm_instance = ChatOpenAI(base_url=base_url, openai_api_key=api_key, **llm_params)
            print(f"LLM has been configured and initialized successfully.")
        except Exception as e:
            raise ValueError(f"Failed to initialize LLM: {e}")


def get_llm():
    """
    Retrieves the globally configured LLM instance.

    Returns:
        object or None: The LLM instance if configured, otherwise None.
    """
    global _llm_instance
    if _llm_instance is None:
        raise RuntimeError("LLM has not been configured. Call BRICK.config_llm(base_url, llm_params) first.")
    return _llm_instance

def clear_llm_config():
    """
    Clears all LLM configurations.  
    This function is useful for debugging or resetting the configuration.
    """
    global _llm_instance
    _llm_instance.clear()

def config_llm_embedding(modeltype='OpenAIEmbeddings', base_url=None, api_key=None, llm_params={'LLM':{}}):
    """
    Configures and initializes a globally unique LLM instance.

    Args:
        modeltype (str, optional): The type of the model. Defaults to 'AzureChatOpenAI'.
        base_url (str, optional): The base URL of the service. Defaults to None.
        llm_params (dict, optional): Parameters for configuring the model. Defaults to {'LLM': {}}.

    Returns:
        AzureChatOpenAI: The initialized LLM instance.
    """
    global _llm_embedding_instance
    if modeltype == 'OpenAIEmbeddings':
        try:
            _llm_embedding_instance = OpenAIEmbeddings(base_url=base_url, openai_api_key=api_key, **llm_params)
            print(f"Embedding model has been configured and initialized successfully.")
        except Exception as e:
            raise ValueError(f"Failed to initialize LLM: {e}")
    elif modeltype == 'ChatOpenAI':
        try:
            _llm_embedding_instance = ChatOpenAI(base_url=base_url, openai_api_key=api_key, **llm_params)
            print(f"Embedding model has been configured and initialized successfully.")
        except Exception as e:
            raise ValueError(f"Failed to initialize LLM: {e}")
    elif modeltype == 'OpenAI':
        try:
            if "model" in llm_params:
                model = llm_params["model"]
            else:
                raise ValueError("model is not specified in llm_params, use llm_params={'model': 'model_name'} to specify it.")
            client = OpenAI(api_key=api_key, base_url=base_url)

            class _OpenAIEmbeddings(Embeddings):
                def embed_query(self, text: str):
                    response = client.embeddings.create(input=text, model=model)
                    return response.data[0].embedding

                def embed_documents(self, texts):
                    batch_size = 10  # OpenAI 限制最大 batch_size=10
                    all_embeddings = []

                    for i in range(0, len(texts), batch_size):
                        batch = texts[i : i + batch_size]  # 按 batch_size 取子集
                        response = client.embeddings.create(input=batch, model=model)
                        batch_embeddings = [item.embedding for item in response.data]
                        all_embeddings.extend(batch_embeddings)  # 累加所有结果

                    return all_embeddings
            
            _llm_embedding_instance = _OpenAIEmbeddings()
            print(f"Embedding model has been configured and initialized successfully.")
        except Exception as e:
            raise ValueError(f"Failed to initialize LLM: {e}")

def get_llm_embedding():
    """
    Retrieves the globally configured LLM instance.

    Returns:
        LLM instance or None: Returns the LLM instance if configured, otherwise returns None.
    """
    global _llm_embedding_instance
    if _llm_embedding_instance is None:
        raise RuntimeError("Embedding model has not been configured. Call BRICK.config_llm_embedding(modeltype='OpenAIEmbeddings', base_url=None, api_key=None, llm_params={'LLM':{}}) first.")
    return _llm_embedding_instance

def clear_llm_embedding():
    """
    Clears all configurations (used for debugging or resetting purposes).
    """
    global _llm_embedding_instance
    _llm_embedding_instance.clear()

def get_openai_embeddings():
    """
    Returns an Embeddings adapter instance as a function.
    """
    embedding_config = get_llm_embedding()  # 获取 OpenAI 配置信息
    model = embedding_config.embeddings.model
    client = OpenAI(api_key=embedding_config.api_key, base_url=embedding_config.base_url)

    class _OpenAIEmbeddings(Embeddings):
        # 适配 OpenAIEmbeddingsAdapter，确保 embedding_model 兼容 FAISS
        def embed_query(self, text: str):
            response = client.embeddings.create(input=text, model=model)
            return response.data[0].embedding

        def embed_documents(self, texts):
            response = client.embeddings.create(input=texts, model=model)
            return [item.embedding for item in response.data]

    return _OpenAIEmbeddings()


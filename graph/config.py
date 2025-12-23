import os
from dotenv import load_dotenv
import json
#from logger import setup_logger
from langchain_openai.chat_models import ChatOpenAI
from langchain_openai.embeddings import OpenAIEmbeddings
from langchain_core.output_parsers import JsonOutputParser, StrOutputParser
#import logging



#logger = setup_logger("config", 'logs/config.log', level=logging.DEBUG)

# 配置环境变量
load_dotenv(dotenv_path='/home/lyt/checker_finallap/brick_test_config.env')
api_key = os.getenv('API_KEY')
base_url = os.getenv('BASE_URL')
url = os.getenv('KG_URL')
auth = (os.getenv('KG_AUTH'), os.getenv('KG_PASS'))
dsr1_params = json.loads(os.getenv('DS_R1', '{}'))
dsv3_params = json.loads(os.getenv('DS_V3', '{}'))
qwen_params = json.loads(os.getenv('QWEN_MAX', '{}'))
emb_params = json.loads(os.getenv('EMBEDDING_MODEL', '{}'))


# 初始化模型，如果API密钥不存在则设为None
try:
    if api_key and base_url:
        model_v3 = ChatOpenAI(base_url=base_url, openai_api_key=api_key, **dsv3_params)
        #logger.info(f"""{dsv3_params["model_name"]} has been configured""")
        model_q = ChatOpenAI(base_url=base_url, openai_api_key=api_key, **qwen_params)
        #logger.info(f"""{qwen_params["model_name"]} has been configured""")
        # 对于嵌入模型，使用emb_params中的配置，避免重复参数
        embedding_model = OpenAIEmbeddings(**emb_params)
    else:
        model_v3 = None
        model_q = None
        embedding_model = None
        print("Warning: API_KEY or BASE_URL not found. Models will be set to None.")
except Exception as e:
    print(f"Warning: Failed to initialize models: {e}")
    model_v3 = None
    model_q = None
    embedding_model = None

output_parser = JsonOutputParser()
str_output_parser = StrOutputParser()
#logger.debug(f"{nodes}\n{relations}")

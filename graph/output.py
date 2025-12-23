from langchain_core.output_parsers import JsonOutputParser, StrOutputParser
from tenacity import retry, stop_after_attempt, wait_fixed
import json
from typing import Generator, Tuple
from langchain.schema import AIMessage

output_parser = JsonOutputParser()
str_output_parser = StrOutputParser()

def stream_chunks(text: str):
    for char in text:
        yield char

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))  # 最多重试3次，每次间隔3秒
def safe_invoke(chain, input_text):
    return chain.invoke(input_text)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3)) 
def safe_invoke_json(chain, input_text):
    result = chain.invoke(input_text)
    return output_parser.invoke(result)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))  
def safe_invoke_str(chain, input_text):
    result = chain.invoke(input_text)
    return str_output_parser.invoke(result)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_invoke_model(model, template):
    return model.invoke(template)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_invoke_model_json(model, template):
    return model.invoke(template)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_invoke_model_str(model, template):
    result = model.invoke(template)
    return str_output_parser.invoke(result)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_invoke_agent_json(agent, input_dict):
    result = agent.invoke(input_dict)
    # 获取最后一条AI消息内容并解析为JSON
    last_ai_message = next(
        (msg for msg in reversed(result["messages"]) 
         if isinstance(msg, AIMessage) or getattr(msg, 'role', None) == 'ai'),
        None
    )
    if last_ai_message and last_ai_message.content:
        try:
            return json.loads(last_ai_message.content)
        except json.JSONDecodeError:
            return {"error": "Invalid JSON format"}
    return {}

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_invoke_agent_str(agent, input_dict):
    result = agent.invoke(input_dict)
    # 获取最后一条AI消息内容
    last_ai_message = next(
        (msg for msg in reversed(result["messages"]) 
         if isinstance(msg, AIMessage) or getattr(msg, 'role', None) == 'ai'),
        None
    )
    return last_ai_message.content if last_ai_message else ""

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream(chain, input_text): 
    return chain.stream(input_text)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream_json(chain, input_text):
    stream = chain.stream(input_text)
    result_str = ""
    for chunk in stream:
        if chunk.content:
            print(chunk.content,end="" if chuck == stream[-1] else "\n", flush=True)
            result_str += chunk.content
    return output_parser.invoke(result_str)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream_str(chain, input_text):
    stream = chain.stream(input_text)
    result_str = ""
    for chunk in stream:
        if chunk.content:
            print(chunk.content,end="" if chuck == stream[-1] else "\n", flush=True)
            result_str += chunk.content
    return str_output_parser.invoke(result_str)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream_model(model, template):
    return model.stream(template)  

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream_model_json(model, template) -> Tuple[Generator[str, None, None], dict]:
    stream = model.stream(template)
    result = {"thought": [], "output": []}
    state = {"last_content": None}

    def generator():
        for chunk in stream:
            content = getattr(chunk, "content", None)
            if not content:
                continue

            state["last_content"] = content

            try:
                parsed = json.loads(content)
                
                if "thought" in parsed:
                    thought = parsed["thought"]
                    result["thought"].append(thought)
                    for part in stream_chunks(thought):
                        yield f"event: thought_stream\ndata: {json.dumps(part, ensure_ascii=False)}\n\n"
                
                if "output" in parsed:
                    output = parsed["output"]
                    result["output"].append(output)
                    for part in stream_chunks(output):
                        yield f"event: output_stream\ndata: {json.dumps(part, ensure_ascii=False)}\n\n"
            
            except json.JSONDecodeError as e:
                yield f"event: output_stream\ndata: {json.dumps('⚠️ JSON解析失败', ensure_ascii=False)}\n\n"
                continue
        
        # 最后统一返回 session 结束信息
        yield f"event: session_complete\ndata: {json.dumps({'complete': True, 'final_answer': ''.join(result['output'])}, ensure_ascii=False)}\n\n"

    return generator(), state, result

def safe_stream_model_json2(model, template):
    stream = model.stream(template) 
    result_str = ""
    for chunk in stream:
        if chunk.content:
            #print(chunk.content,end="" if chuck == stream[-1] else "\n", flush=True)
            result_str += chunk.content
    return output_parser.invoke(result_str)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream_model_str(model, template):
    stream = model.stream(template)  
    result_str = ""
    for chunk in stream:
        if chunk.content:
            #print(chunk.content,end="" if chuck == stream[-1] else "\n", flush=True)
            result_str += chunk.content
    return str_output_parser.invoke(result_str)

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream_agent_json(agent, input_dict) -> Tuple[Generator[str, None, None], dict]:
    stream = agent.stream(input_dict)
    result = {"thought": [], "output": []}
    state = {"last_content": None}

    def generator():
        for chunk in stream:
            if isinstance(chunk, dict) and 'messages' in chunk:
                ai_messages = [msg for msg in chunk['messages'] if isinstance(msg, AIMessage)]
                if not ai_messages:
                    continue

                last_msg = ai_messages[-1]
                content = last_msg.content
                if not content:
                    continue

                state["last_content"] = content
                
                try:
                    parsed = json.loads(content)

                    if "thought" in parsed:
                        thought = parsed["thought"]
                        print(f"[Thought] {thought}", flush=True)
                        result["thought"].append(thought)
                        for part in stream_chunks(thought):
                            yield f"event: thought_stream\ndata: {json.dumps(part, ensure_ascii=False)}\n\n"

                    if "output" in parsed:
                        output = parsed["output"]
                        print(f"[Output] {output}", flush=True)
                        result["output"].append(output)
                        for part in stream_chunks(output):
                            yield f"event: output_stream\ndata: {json.dumps(part, ensure_ascii=False)}\n\n"

                except json.JSONDecodeError as e:
                    print(f"[Error] Invalid JSON:\n{content}\n{e}", flush=True)
                    yield f"event: output_stream\ndata: {json.dumps('⚠️ JSON解析失败', ensure_ascii=False)}\n\n"

        yield f"event: session_complete\ndata: {json.dumps({'complete': True, 'final_answer': ''.join(result['output'])})}\n\n"

    return generator(), state, result

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def safe_stream_agent_json2(agent, input_dict):
    stream = agent.stream(input_dict)

    result = {
        "thought": [],
        "output": []
    }

    for chunk in stream:
        if isinstance(chunk, dict) and 'messages' in chunk:
            ai_messages = [msg for msg in chunk['messages'] if isinstance(msg, AIMessage)]
            if not ai_messages:
                continue

            last_msg = ai_messages[-1] 
            content = last_msg.content
   
            if not content:
                continue

            try:
                parsed = json.loads(content)
                thought = parsed.get("thought")
                output = parsed.get("output")

                if thought:
                    print(f"[Thought] {thought}", flush=True)
                    result["thought"].append(thought)
                if output:
                    print(f"[Output] {output}", flush=True)
                    result["output"].append(output)

            except json.JSONDecodeError as e:
                print(f"[Error] Invalid JSON:\n{content}\n{e}", flush=True)

    return content, result

@retry(stop=stop_after_attempt(3), wait=wait_fixed(3))
def stream_agent_response_to_frontend(agent, input_dict) -> Generator[str, None, None]:
    """
    多事件 SSE 推送 agent 输出（支持思考、输出、agent 切换等）。
    """
    stream = agent.stream(input_dict)

    current_agent = None
    buffer_thought = ""
    buffer_output = ""

    for chunk in stream:
        if isinstance(chunk, dict) and 'messages' in chunk:
            for msg in chunk['messages']:
                if isinstance(msg, AIMessage):
                    # 1. 如果 agent 切换了，发 agent_info 事件
                    if msg.name and msg.name != current_agent:
                        current_agent = msg.name
                        yield f"event: agent_info\ndata: {json.dumps({'agent_type': current_agent, 'requires_interaction': True})}\n\n"

                    # 2. 尝试从 content 中提取 JSON
                    content = msg.content
                    if not content:
                        continue

                    try:
                        parsed = json.loads(content)

                        # 流式输出 thought
                        if "thought" in parsed:
                            thought = parsed["thought"]
                            for chunk in stream_chunks(thought):
                                yield f"event: thought_stream\ndata: {json.dumps(chunk, ensure_ascii=False)}\n\n"
                            buffer_thought += thought

                        # 流式输出 output
                        if "output" in parsed:
                            output = parsed["output"]
                            for chunk in stream_chunks(output):
                                yield f"event: output_stream\ndata: {json.dumps(chunk, ensure_ascii=False)}\n\n"
                            buffer_output += output

                    except json.JSONDecodeError:
                        # 如果 JSON 解码失败，可以忽略或上报
                        yield f"event: output_stream\ndata: {json.dumps('⚠️ JSON解析失败', ensure_ascii=False)}\n\n"

    # 3. 单个 agent 结束
    yield f"event: agent_complete\ndata: {json.dumps({'complete': True})}\n\n"

    # 4. 整个任务结束（可以返回最终结果）
    yield f"event: session_complete\ndata: {json.dumps({'complete': True, 'final_answer': buffer_output}, ensure_ascii=False)}\n\n"

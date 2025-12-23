import time

_last_id = None


def build_config_id() -> str:
    """根据当前时间戳生成五位数字字符串 ID，尽量避免重复"""
    global _last_id
    # 毫秒级时间戳
    ts_ms = int(time.time() * 1000)
    candidate = ts_ms % 100000  # 取后 5 位

    # 若与上一次生成的相同（同一毫秒多次调用），则顺延一位
    if _last_id is not None and candidate == _last_id:
        candidate = (candidate + 1) % 100000

    _last_id = candidate
    return f"{candidate:05d}"
# ... existing code ...

if __name__ == "__main__":
    print(build_config_id())

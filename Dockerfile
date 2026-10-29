# ClinVar VEP Benchmark Audit -- reproducible runtime
FROM python:3.11-slim

# 1) 工作目录（相当于容器内的 cd /app）
WORKDIR /app

# 2) 先只拷依赖文件。这一步单独一层，后面改代码不会让 pip 重装（层缓存）
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 3) 再拷代码和数据
COPY . .

# 4) 默认命令：跑完整基线矩阵
CMD ["python", "evaluate.py", "--all"]

## BRICK Agent - User Guide

### Step 1: Clone the Repository

```bash
git clone -b master https://github.com/BGIResearch/BiOmics
```

### Step 2: Configure Environment Variables

Edit the `graph/brick_test_config.env` file and modify the following configuration:

```bash
# Change the project root path to your actual path
PROJECT_ROOT=/your/path/to/Biomics_agent
```

### Step 3: Create Conda Environment

Create a Python environment using the provided configuration file:

```bash
# Create environment using conda
conda env create -f biomics_environment.yml

# Activate the environment
conda activate biomics_agent
```

### Step 4: Start Docker Sandbox Container

#### 4.1 Download Docker Image

Download the pre-built image from Alibaba Cloud Container Registry:

```bash
# Download image (approximately 6.57GB)
docker pull crpi-b88i7r04wqgzpar4.cn-beijing.personal.cr.aliyuncs.com/biomics/biomics-agent:v6

# Optional: Tag the image
docker tag crpi-b88i7r04wqgzpar4.cn-beijing.personal.cr.aliyuncs.com/biomics/biomics-agent:v6 biomics_agent:v6
```

#### 4.2 Start the Sandbox Container

```bash
# Create a data directory under the project folder:
Biomics_agent/data

# Stop and remove old container (if exists)
docker rm -f my_code_sandbox

# Start new sandbox container
docker run -d \
  --name my_code_sandbox \
  --network="host" \
  -e PYTHONUNBUFFERED=1 \
  -v /your/path/to/Biomics_agent/data:/workspace/data \
  biomics_agent:v6 \
  jupyter kernelgateway \
  --KernelGatewayApp.ip=0.0.0.0 \
  --KernelGatewayApp.port=8888 \
  --KernelGatewayApp.auth_token="" \
  --JupyterWebsocketPersonality.list_kernels=True \
  --KernelManager.ip=0.0.0.0
```

**Notes:**
- Replace `/your/path/to/Biomics_agent/data` with the **absolute path** to your project's `data` directory
- The container will start Jupyter Kernel Gateway service on port `8888`
- Using `--network="host"` allows the container to share the host's network

**Verify container is running:**
```bash
# Check container status
docker ps | grep my_code_sandbox

# Test connection
curl http://127.0.0.1:8888/api
```

### Step 5: Run the Application

Start the NiceGUI Web application:

```bash
# Make sure you're in the project root directory
cd /your/path/to/Biomics_agent

# Activate the environment (if not already activated)
conda activate biomics_agent

# Start the application
python app_nicegui.py
```

Once the application starts, the browser will automatically open or you can visit:
```
http://localhost:8080
```

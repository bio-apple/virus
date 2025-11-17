#!/bin/bash

# 定义安装路径 (变量 1)
# 默认值: /staging/explify_china
BASE_DIR=${1:-/staging/explify_china}

# 定义应用端口 (变量 2)
# 默认值: 80
APP_PORT=${2:-80}

# 定义 Docker 镜像加速器地址 (变量 3)
# 默认值: https://docker-0.unsee.tech
DOCKER_MIRROR=${3:-https://docker-0.unsee.tech}

# 检查 BASE_DIR 是否设置
if [ -z "$BASE_DIR" ]; then
    echo "Error: BASE_DIR is not set. Please provide it as the first argument."
    exit 1
fi

echo "--- 脚本配置概览 ---"
echo "BASE_DIR: $BASE_DIR"
echo "APP_PORT: $APP_PORT"
echo "DOCKER_MIRROR: $DOCKER_MIRROR"
echo "----------------------"

# 创建主目录
mkdir -p "$BASE_DIR"
# 创建用户指定的子目录
mkdir -p "$BASE_DIR/customer_data"
mkdir -p "$BASE_DIR/Docker"
mkdir -p "$BASE_DIR/explify-databases"
mkdir -p "$BASE_DIR/Knowledge"
mkdir -p "$BASE_DIR/microbe"
mkdir -p "$BASE_DIR/script"

# 1. 安装 Java
echo "👉 1. Installing Java 11..."
yum install -y java-11-openjdk.x86_64 java-11-openjdk-devel.x86_64
alternatives --config java

# 2. 安装 LibreOffice
echo "👉 2. Installing LibreOffice RPMs..."
# 确保目录存在，尽管原始脚本中没有 mkdir，但 cd 命令要求目录存在
mkdir -p "$BASE_DIR/microbe/install/LibreOffice_7.5.3.2_Linux_x86-64_rpm/RPMS"
cd "$BASE_DIR/microbe/install/LibreOffice_7.5.3.2_Linux_x86-64_rpm/RPMS"
yum install -y *.rpm

# 3. 安装字体
echo "👉 3. Installing Chinese fonts..."
mkdir -p /usr/share/fonts/SourceHanSans
cp "$BASE_DIR/microbe/install/OTF/SimplifiedChinese/*" /usr/share/fonts/SourceHanSans
mkdir -p /usr/share/fonts/yueyuan
cp "$BASE_DIR/microbe/install/OTF/MFYueYuan_Noncommercial-Regular.otf" /usr/share/fonts/yueyuan
fc-cache -fv

# 4. 配置 Docker 镜像加速器和重启 Docker 服务
echo "👉 4. Configuring Docker registry mirror and restarting Docker..."

# 创建 Docker 配置文件目录
mkdir -p /etc/docker/

# 写入 daemon.json 配置（使用 DOCKER_MIRROR 变量）
tee /etc/docker/daemon.json <<EOF
{
    "registry-mirrors": [
        "$DOCKER_MIRROR"
    ]
}
EOF

# 重新加载 systemd 配置并重启 Docker 服务
systemctl daemon-reload && systemctl restart docker

# 4b. Docker (PostgreSQL) setup - 原始脚本中被注释的部分
echo "Setting up PostgreSQL 14.2 via Docker (uncommented)..."
#docker pull docker.io/postgres:14.2
docker run --name postgres14 -v "$BASE_DIR/microbe/data/postgres/data:/var/lib/postgresql/data" -e POSTGRES_PASSWORD=happy -d -p 5432:5432 docker.io/postgres:14.2
docker update --restart=always postgres14

# 5. 防火墙配置 (使用变量 APP_PORT)
echo "👉 5. Configuring firewall for port $APP_PORT/tcp..."
firewall-cmd --zone=public --add-port=$APP_PORT/tcp --permanent
firewall-cmd --reload

# 6. Systemd 服务文件创建
SERVICE_FILE="/usr/lib/systemd/system/IlluminaMicrobe.service"
echo "👉 6. Creating systemd service file: $SERVICE_FILE..."
cat <<EOF > "$SERVICE_FILE"
[Unit]
Description=IlluminaMicrobe
After=syslog.target network.target remote-fs.target nss-lookup.target
Requires=docker.service

[Service]
#LimitCORE=infinity
LimitNOFILE=100000
LimitNPROC=100000
Type=simple
WorkingDirectory=$BASE_DIR/microbe
ExecStart=$BASE_DIR/microbe/startup.sh
ExecReload=/bin/kill -s HUP \$MAINPID
ExecStop=/bin/kill -s QUIT \$MAINPID
PrivateTmp=true

[Install]
WantedBy=multi-user.target
EOF

# 7. 启用和启动服务
echo "👉 7. Enabling and starting services..."
systemctl enable docker.service
systemctl enable IlluminaMicrobe.service
systemctl start IlluminaMicrobe.service

echo "--- 脚本执行完成.---"
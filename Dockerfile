# Start with Python 3.10 image
FROM python:3.10

# Install C++ compilation tools
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    git \
    g++ \
    gdb \
    clang \
    make \
    ninja-build \
    cmake \
    autoconf \
    automake \
    libtool \
    valgrind \
    locales-all \
    dos2unix \
    rsync \
    tar \
    python3 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN python -m pip install -r requirements.txt

WORKDIR /app
COPY . .

CMD ["python", "MolecularDynamics.py"]
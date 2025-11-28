FROM ubuntu:24.04

# Install minimal dependencies
RUN apt-get update && apt-get install -y \
    git cmake build-essential clang libgmp-dev bash util-linux \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy source and scripts
COPY src/ src/
COPY include/ include/
COPY tests/ tests/
COPY scripts/ scripts/
COPY CMakeLists.txt .
COPY CMakePresets.json .

RUN chmod +x scripts/*.sh

# Install google benchmark
WORKDIR /app/benchmark
RUN git clone https://github.com/google/benchmark.git . \
    && cmake -E make_directory "build" \
    && cmake -E chdir "build" cmake -DBENCHMARK_DOWNLOAD_DEPENDENCIES=on -DCMAKE_BUILD_TYPE=Release ../ \
    && cmake --build "build" --config Release

# Default command: build + run benchmarks
WORKDIR /app
CMD ["bash", "-c", "scripts/build.sh && scripts/run_benchmarks.sh"]

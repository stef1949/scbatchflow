FROM python:3.12-slim

LABEL org.opencontainers.image.title="scbatchflow" \
      org.opencontainers.image.source="https://github.com/stef1949/scbatchflow" \
      org.opencontainers.image.licenses="MIT"

# fastqc for read QC; procps because Nextflow needs `ps` to collect task metrics
RUN apt-get update \
    && apt-get install -y --no-install-recommends fastqc procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/scbatchflow
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir . multiqc

# numba and matplotlib need writable caches when the container runs as a non-root user
ENV NUMBA_CACHE_DIR=/tmp/numba MPLCONFIGDIR=/tmp/matplotlib

CMD ["scbf", "--help"]

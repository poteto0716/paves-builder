# syntax=docker/dockerfile:1
# Ubuntu 24.04 matches the public Linux x86_64 wheel's runtime ABI.
FROM mambaorg/micromamba:2.9.0-ubuntu24.04@sha256:d586ee0cd12ffd035e36ffb28f245b0288de088f44a77a2e4e4d6feb670299e5

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        ffmpeg fonts-dejavu-core libdbus-1-3 libegl1 libgl1 libglib2.0-0t64 \
        libopengl0 libvulkan1 libx11-xcb1 libxcb-cursor0 libxi6 libxkbcommon0 \
        libxrender1 mesa-vulkan-drivers \
    && rm -rf /var/lib/apt/lists/*

COPY docker/jupyter_server_config.py /etc/jupyter/jupyter_server_config.py
COPY docker/healthcheck.py /usr/local/share/paves/healthcheck.py
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/environment.yml /tmp/environment.yml
# OpenMM has both CPU and CUDA platforms. The solver need not see a GPU at build time.
RUN CONDA_OVERRIDE_CUDA=12.9 micromamba install --yes --name base \
        --strict-channel-priority --file /tmp/environment.yml \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER dist/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl /tmp/
RUN python -m pip install --no-deps --no-cache-dir /tmp/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl \
    && python -m pip install --no-cache-dir 'ovito==3.16.1' \
    && rm /tmp/polypaves-0.4.0-cp311-cp311-linux_x86_64.whl \
    && python -m ipykernel install --sys-prefix --name python3 --display-name 'PAVES (Python 3.11)'

# Runtime audit tools used by the public distribution tests; no compiler is installed.
USER root
RUN apt-get update && apt-get install -y --no-install-recommends binutils file unzip \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Only public runtime assets are copied. Building does not require the private source checkout.
COPY --chown=$MAMBA_USER:$MAMBA_USER examples/ /workspace/examples/
COPY --chown=$MAMBA_USER:$MAMBA_USER external/ /workspace/external/
COPY --chown=$MAMBA_USER:$MAMBA_USER workflows/ /workspace/workflows/
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/ /workspace/docker/
COPY --chown=$MAMBA_USER:$MAMBA_USER tests/test_public_reactive_workflows.py /workspace/tests/test_public_reactive_workflows.py
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/audit_binary_wheel.sh /workspace/scripts/audit_binary_wheel.sh
COPY --chown=$MAMBA_USER:$MAMBA_USER LICENSE THIRD_PARTY_LICENSES.md /workspace/
COPY --chown=$MAMBA_USER:$MAMBA_USER LICENSES/ /workspace/LICENSES/

# OVITO selects its own headless Qt platform plugin, which also supports Vulkan.
ENV MPLBACKEND=Agg \
    XDG_CACHE_HOME=/tmp/paves-cache \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1
WORKDIR /workspace
EXPOSE 8888
HEALTHCHECK --interval=30s --timeout=10s --start-period=90s --retries=3 \
    CMD ["/opt/conda/bin/python", "/usr/local/share/paves/healthcheck.py"]
# Keep micromamba's entrypoint so `compose run ... python ...` activates the environment too.
CMD ["jupyter", "lab", "--config=/etc/jupyter/jupyter_server_config.py"]

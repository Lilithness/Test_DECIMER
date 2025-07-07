# ---- 1. base image ----------------------------------------------------------
FROM python:3.10-slim

# ---- 2. system libs ---------------------------------------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git build-essential poppler-utils \
        libgl1 libglib2.0-0 libsm6 && \
    rm -rf /var/lib/apt/lists/*

# Configure Git to use HTTPS without prompting
RUN git config --global url."https://github.com/".insteadOf git@github.com: && \
    git config --global advice.detachedHead false
# ---- 3. Python core deps ----------------------------------------------------
RUN pip install --no-cache-dir \
        "tensorflow==2.12.0"  "numpy==1.23.5" \
        opencv-python pillow pillow-heif pdf2image tqdm rdkit-pypi

# ---- 4. clone & editable-install the three DECIMER repos --------------------
WORKDIR /opt/decimer-src
RUN git clone https://github.com/Kohulan/DECIMER-Image-Segmentation.git && \
    git clone https://github.com/Iagea/DECIMER-Image-Classifier.git && \
    git clone https://github.com/Kohulan/DECIMER-Image_Transformer.git && \
    pip install --no-cache-dir -e DECIMER-Image-Segmentation \
                               -e DECIMER-Image-Classifier \
                               -e DECIMER-Image_Transformer
# ---- 5. copy your pipeline scripts -----------------------------------------
WORKDIR /app
COPY scripts/ /app/scripts/
COPY data/ /app/data/
COPY results/ /app/results
COPY test /app/test

# ---- 7. default entrypoint --------------------------------------------------
ENTRYPOINT ["bash", "/app/handler.sh"]


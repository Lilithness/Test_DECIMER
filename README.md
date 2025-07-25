# Test_DECIMER


# Patent-to-SMILES Pipeline (DECIMER)

End-to-end extraction of chemical structures from PDF patents using
DECIMER-Segmentation, DECIMER-Image-Classifier and the
DECIMER Image-Transformer (OCSR).

---

## 1. Quick start with Docker

```bash
git clone https://github.com/Lilithness/Test_DECIMER.git
cd Test_DECIMER

# build
docker build -t decimer-pipeline .

# run (data/ contains <patent>/pdf/*.pdf or page PNGs)
docker run --rm -it  --env-file test.env  -v ./data:/app/data   -v ./results:/app/results   decimer-pipeline /bin/bash


You will need env file, where PATENT should be either on pdf (just basename.pdf) or a txt file with \n separated list of pdfs
All pdfs should be in data/patents
Individual csv of results (per patent) will be saved in results/ folder 
```
## 2. Output

``` text
results/
 ├── <PATENT>_ocsr_results.csv          # Patent,SMILES,valid
 └── <PATENT>/filtered_segments/*.png
```

---

## 3. Folder layout


```text
data/
 └── patents/
     └── US11066404.pdf         
scripts/
 ├── pdf2png.py
 ├── segment_pngs.py
 ├── classify.py
 ├── ocsr_to_csv.py
 ├── run.sh
 └── handler.sh
```

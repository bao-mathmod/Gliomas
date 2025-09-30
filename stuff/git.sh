# 0) Go to your project root
cd /mnt/18T/chibao/gliomas/code

# 1) Init (ok if already initialized)
git init

# 2) Create .gitignore (CLOSES with GIT)
cat > .gitignore <<'GIT'
# OS / editors
.DS_Store
Thumbs.db
.vscode/
.idea/

# Python / R / notebooks
__pycache__/
*.pyc
.ipynb_checkpoints/
.Rproj.user/

# Environments, builds, logs
.venv/
env/
*.egg-info/
dist/
build/
logs/
*.log

# Nextflow / workflow caches
.work/
.nf-cache/
.singularity_cache/
.results/

# Large data / outputs (adjust to your layout)
data/
results/
GIT

# 3) Create README (CLOSE code fence ``` and MD)
cat > README.md <<'MD'
# Gliomas Pipeline

Reproducible pipelines (Nextflow) for downloading, QC, and processing sc/sn data.

## Quickstart
```bash
# create env, for example
# conda env create -f env.yml && conda activate glio

nextflow run main.nf -profile slurm

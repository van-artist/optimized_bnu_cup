chmod +x ./env_init.sh
chmod +x ./run_pipline.sh

conda init
source ~/miniconda3/etc/profile.d/conda.sh
conda activate bio

bash ./env_init.sh
bash ./run_pipline.sh

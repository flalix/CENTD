# CENTD
CENTD - Centro de Excelência para Descobertas de Alvos  - Instituto Butantan

### Portal em fase de teste

### Aplicativos

  - Resultados:
    - Biobanco
  - SARS-CoV-2 por país e tempo:
    - Entropia mutacional
    - Mapa de mutações
  - WebService
    - controle de temperatura


### Flask + Dash

Dash embeded in Flask:  
  - How to embed a Dash app into an existing Flask app

see: https://medium.com/@olegkomarov_77860/how-to-embed-a-dash-app-into-an-existing-flask-app-ea05d7a2210b


### Github CENTD

git rm --cached -r migrations/  
git rm --cached -r logs/  
git rm --cached -r venv/  

git add *  
git commit -m 'Dash embedded in Flask'  
git push  

### gitignore

git config --global core.excludesfile ./gitignore  

\# Cache
__pycache__

\# OS  
.DS_Store  

\# IDE  
.vscode  

\# Settings  
.envrc  

logs/  
migrations/  
venv/  

\*.db  

### Moving data to server

#### biobank screening - moving data to server

server='cetics@172.25.1.72:'  
root_centd=$server'/projects/web/CENTD/'  
root_screening=$server'/projects/colaboracoes/biobanco/screening/'  
echo $root_screening  

rsync -v \*.xls? $root_screening  

#### sarscov2-gisaid - moving data to server

  - server
cd /projects/colaboracoes/  
mkdir covid/  
cd covid  
mkdir fasta/  
cd fasta/  
mkdir sarscov2_202007  
cd sarscov2_202007  
mkdir msa_0713_202007  
cd msa_0713_202007  
mkdir trees  
mkdir protein  
cd protein  
mkdir html  
mkdir figure  
mkdir entropy  

  - local
root_local=/media/flalix/5c1ba0b4-f897-451c-9068-ac5e57194590/flalix/  
root_local_covid=$root_local'colaboracoes/covid/fasta/sarscov2_202007/'  
echo $root_local_covid  
cd  $root_local_covid  

root_fasta_covid=$server'/projects/colaboracoes/covid/fasta/sarscov2_202007/'  
echo $root_fasta_covid  

cd  $root_local_covid  
rsync -v metadata.tsv $root_fasta_covid  
cd msa_0713_202007/  
rsync -rv . $root_fasta_covid'msa_0713_202007/'  



### Install Locally

#--- install virtual environment
sudo apt-get install python3-venv

  - Anaconda3
source activate py_env36

  - Install virtual environment python (3.6 for microblog/centd)
python -m venv venv

  - Activate venv
source venv/bin/activate

   - Confirm python version
python --version

   - Don't forget to upgrade pip (for this version)
pip install --upgrade pip

   - install from requirement file
pip install -r requirements.txt   

pip freeze > requirements2.txt  

pip install -r requirements_bokeh.txt

pip freeze > requirements_microblog_bokeh.txt

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

### Server / Servidor

sudo snx -s 143.107.75.58 -u flavio.lichtenstein  
pwd: butantan@123  

terminator -p remotehost  

\#-- end --  
sudo snx -d  

# ssh cetics@172.25.1.72
terminator -p remotehost

cd /projects/web/CENTD

cat ~/.config/terminator/config

sudo cp ~/.config/terminator/config.centd   ~/.config/terminator/config
sudo cp ~/.config/terminator/config.bioinfo ~/.config/terminator/config

ssh -p 2200 flavio_lichtenstein@172.25.2.25

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

### Docker - installing CENTD

docker image list
docker container list
  - was active: microblog, elasticsearch, mysql
  - removing microblog

id=
docker stop $id
docker ps

rem - git checkout https://github.com/flalix/CENTD
docker build -t centd:latest .
docker image list

docker run --name centd -d -p 8000:5000 --rm centd:latest
docker ps

id=xxxxx
docker stop $id
docker ps

docker run --name mysql -d -e MYSQL_RANDOM_ROOT_PASSWORD=admin -e MYSQL_DATABASE=centd -e MYSQL_USER=root -e MYSQL_PASSWORD=admin mysql/mysql-server:5.7

docker ps

id=xxxxx
docker stop $id
docker ps


docker run --name centd -d -p 8000:5000 --rm -e SECRET_KEY=centdSecret#345 \
-e MAIL_SERVER=smtp.googlemail.com  -e MAIL_PORT=587  -e MAIL_USERNAME=flalix@gmail.com -e MAIL_USE_TLS=true  -e MAIL_PASSWORD=Flavi3131 \
--link mysql:dbserver -e DATABASE_URL=mysql+pymysql://centd:admin@dbserver/centd  \
centd:latest

docker ps

id=xxxxx
docker stop $id
docker ps

docker run --name elasticsearch -d --rm -e "discovery.type=single-node" docker.elastic.co/elasticsearch/elasticsearch-oss:6.1.1

docker run --name centd -d -p 8000:5000 --rm -e SECRET_KEY=centdSecret#345 \
-e MAIL_SERVER=smtp.googlemail.com  -e MAIL_PORT=587  -e MAIL_USERNAME=flalix@gmail.com -e MAIL_USE_TLS=true -e MAIL_PASSWORD=Flavi3131 \
--link mysql:dbserver -e DATABASE_URL=mysql+pymysql://root:admin@dbserver/centd \
--link elasticsearch:elasticsearch -e ELASTICSEARCH_URL=http://elasticsearch:9200 \
-v "$(pwd)"/../src:/src:ro \
centd:latest

docker ps

### MySql

#-- setup MySql database
mysql -u root -p
pwd: admin

CREATE USER 'centd'@'localhost' IDENTIFIED BY 'admin';
CREATE USER 'centd'@'172.17.0.3' IDENTIFIED BY 'admin';
GRANT ALL ON *.* TO 'centd'@'172.17.0.3';
flush privileges;

select User, Host from mysql.user;

# CREATE USER 'myuser'@'localhost' IDENTIFIED BY 'mypass';
# GRANT ALL ON *.* TO 'myuser'@'localhost';

mysql -u 'centd'@'172.17.0.3' -p


SELECT DATABASE();

SHOW TABLES;


### To run in Terminal

172.25.1.72:5000

  - see logs:
docker logs centd


### Run docker command

  - src: parallel to centd
  - colaboracoes: root
  - specific data have been copied to server /colaboracoes

  - config.py
MYSQL_USER='root'  
MYSQL_ROOT_PASSWORD='admin'  


root=/projects/web
root_prj=/projects

\#mysql/mysql-server  5.7   9c31a29b3f30
docker image rm -f 9c31a29b3f30

docker container list --all | grep mysql  
# dd53d2acb996    mysql/mysql-server:5.7     "/entrypoint.sh mysq…"   3 weeks ago     
docker container stop dd53d2acb996
docker container rm dd53d2acb996

docker run --name mysql -d -e MYSQL_RANDOM_ROOT_PASSWORD=admin -e MYSQL_DATABASE=centd -e MYSQL_USER=root -e MYSQL_PASSWORD=admin mysql/mysql-server:5.7


docker ps


docker run --name centd -d -p 8000:5000 --rm -e SECRET_KEY=centdSecret#345 \
-e MAIL_SERVER=smtp.googlemail.com  -e MAIL_PORT=587  -e MAIL_USERNAME=flalix@gmail.com -e MAIL_USE_TLS=true -e MAIL_PASSWORD=Flavi3131 \
--link mysql:dbserver -e DATABASE_URL=mysql+pymysql://root:admin@dbserver/centd \
--link elasticsearch:elasticsearch -e ELASTICSEARCH_URL=http://elasticsearch:9200 \
-v $(pwd)"../src:/home/src" \
-v $(pwd)"../../colaboracoes/covid/fasta/sarscov2_202007:/colaboracoes/covid/fasta/sarscov2_202007" \
-v $(pwd)"../../colaboracoes/biobanco/screening:/colaboracoes/biobanco/screening" \
centd:latest

docker ps
id=

docker cp boot.sh $id:/home/centd
docker cp config.py $id:/home/centd

\#-- faltaram os templates
docker cp ../../colaboracoes/covid/templates $id:/colaboracoes/covid  
\#-- faltaram algumas tabelas como condrocitos  
docker cp ../../colaboracoes/biobanco/screening/ $id:/colaboracoes/biobanco  

../../colaboracoes/biobanco/screening/condrocytes_paola_final.xls
../../colaboracoes/biobanco/screening/condrocytes_paola_final.xls

pip install xlrd

import pandas as pd
fname = '../../colaboracoes/biobanco/screening/condrocytes_paola_final.xls'
import os
os.path.exists(fname)
pd.read_excel(fname)


docker exec -it $id bash

source venv/bin/activate

pip uninstall biopython
pip install biopython==1.76
pip install email-validator
pip install ipython



docker exec -it $id source boot.sh


# TMUX

\#-- https://tmuxcheatsheet.com/
\#-- https://tmuxguide.readthedocs.io/en/latest/tmux/tmux.html

tmux ls

tmux new -s mysession

tmux kill-session -a -t mysession

CTRL+b d

tmux attach-session


# source boot.sh

source venv/bin/activate

export FLASK_APP=centd.py
export FLASK_DEBUG=0
export MAIL_SERVER=localhost
export MAIL_PORT=8025

flask db upgrade
flask translate compile
exec gunicorn -b :5000 --access-logfile - --error-logfile - centd:app





### copying files to container
id=76d29ef7a5c9
docker exec -it $id bash

docker cp  ../src $id:/home/
docker cp  ../../colaboracoes/covid/fasta/sarscov2_202007 $id:/colaboracoes/covid/fasta/
docker cp  ../../colaboracoes/biobanco/screening $id:/colaboracoes/biobanco/
docker cp boot.sh $id:/home/centd

docker exec -it $id bash

docker exec -it $id "ls -ls"  ???  

docker logs centd


export FLASK_APP=centd.py
export FLASK_DEBUG=0
export MAIL_SERVER=localhost
export MAIL_PORT=8025


### Run manually - the container

  - run manually
docker ps  
id=
docker exec -it $id bash   

  - inside the container:
centd@6c40fdbe834f:~$  
source venv/bin/activate  
flask db upgrade  
flask translate compile  
exec gunicorn -b :5000 --access-logfile - --error-logfile - centd:app  

  -- run driver externally
docker exec -it $id bash /source/driver.sh




### Stop container

docker ps
docker container stop xxxxxxx


### Dockerfile

FROM python:3.8-slim  

LABEL MAINTAINER="Flavio Lichtenstein <flalix@gmail.com>"  

\# RUN adduser -D centd  
RUN   adduser --system centd  

WORKDIR /home/centd  

\# https://dev.to/faizanbashir/building-python-data-science-container-usingdocker-3f8p  
\# Building Python Data Science Container using Docker  Faizan Bashir 20 de jan. de 2019  
\# SOFTWARE PACKAGES  
\#   * musl: standard C library  
\#   * lib6-compat: compatibility libraries for glibc  
\#   * linux-headers: commonly needed, and an unusual package name from Alpine.  
\#   * build-base: used so we include the basic development packages (gcc)  
\#   * bash: so we can access /bin/bash  
\#   * git: to ease up clones of repos  
\#   * ca-certificates: for SSL verification during Pip and easy_install  
\#   * freetype: library used to render text onto bitmaps, and provides support font-related operations  
\#   * libgfortran: contains a Fortran shared library, needed to run Fortran  
\#   * libgcc: contains shared code that would be inefficient to duplicate every time as well as auxiliary helper routines and runtime support  
\#   * libstdc++: The GNU Standard C++ Library. This package contains an additional runtime library for C++ programs built with the GNU compiler  
\#   * openblas: open source implementation of the BLAS(Basic Linear Algebra Subprograms) API with many hand-crafted optimizations for specific processor types  
\#   * tcl: scripting language  
\#   * tk: GUI toolkit for the Tcl scripting language  
\#   * libssl1.0: SSL shared libraries  

ENV PACKAGES="\  
    dumb-init \  
    musl \  
    libc6-compat \  
    linux-headers \  
    build-base \  
    bash \  
    git \  
    ca-certificates \  
    freetype \  
    libgfortran \  
    libgcc \  
    libstdc++ \  
    openblas \  
    tcl \  
    tk "  

\# apt-cache search openblas  
\# libblas-test - Basic Linear Algebra Subroutines 3, testing programs  
\# libopenblas-base - Optimized BLAS (linear algebra) library based on GotoBLAS2  
\# libopenblas-dev - Optimized BLAS (linear algebra) library based on GotoBLAS2  

RUN apt-get -y update && apt-get upgrade  
\# RUN apt-get add libssl1.0  
RUN apt-get install -y python3 python3-dev build-essential  
RUN apt-get install -y musl-dev libc6-dev  
RUN apt-get install -y libopenblas-dev libopenblas-base libfreetype6-dev pkg-config  
RUN apt-get install -y dumb-init musl libc6-dev  
RUN apt-get install -y bash git ca-certificates freetype2-demos  
RUN apt-get install -y gfortran-7 aptitude libstdc++6 libopenblas-dev tcl  
\# RUN apt-get install -y linux-headers-$(uname -r)  


COPY requirements.txt requirements.txt  
RUN python -m venv venv  
RUN venv/bin/pip install --upgrade pip  
RUN venv/bin/pip install -r requirements.txt  
RUN venv/bin/pip install gunicorn pymysql  

COPY app app  
COPY migrations migrations  
COPY centd.py config.py boot.sh ./  
RUN chmod a+x boot.sh  

ENV FLASK_APP cented.py  

RUN chown -R centd:centd ./  
USER centd  

EXPOSE 5000  
ENTRYPOINT ["./boot.sh"]  


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

  - SRC (python - my libs)  
root_local=/media/flalix/5c1ba0b4-f897-451c-9068-ac5e57194590/flalix/   
root_src=$root_local'python/src/'  
echo $root_src  
cd $root_src  

root_src_server=$server'/projects/web/src/'   
rsync -v * $root_src_server  



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



diego di456
catarina cat678
olga olga864
anamarisa ana0567
isadora isa444
edu edu459

usuario: isadora
senha: isa444


Iniciando Site de Resultados de Experimentos do CENTD
tem que estar no I Butantan, ou via VPN
http://172.25.1.72:8000/
usuario:
pwd:

entre com login
na tela inicial: Screening/Resultados de Agosto/2020 <- clique aqui
escolha os experimentos
tem um combo acima para mudar de Anti para Pró
Use o botão de voltar do browser para voltar

em mutações do SARS-CoV-2 as datas estão como mes/dia/ano
ainda em testes ... dados da Australia estão errados, estamos corrigindo


ainda em testes
pode comentar .. e sugerir
abr
Flavio


em mutações as datas estão como mes/dia/ano


veja se acessa, dentro do IBu, ou via VPN

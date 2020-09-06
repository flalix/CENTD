source activate py_env36
source venv/bin/activate
export FLASK_APP=centd.py

export FLASK_DEBUG=0

export MAIL_SERVER=localhost
export MAIL_PORT=8025

echo $MAIL_SERVER  $MAIL_PORT
echo $FLASK_DEBUG
echo $FLASK_APP
echo done

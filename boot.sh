#!/bin/sh
# this script is used to boot a Docker container
source venv/bin/activate

export FLASK_APP=centd.py
export FLASK_DEBUG=0
export MAIL_SERVER=localhost
export MAIL_PORT=8025

while true; do
    flask db upgrade
    if [[ "$?" == "0" ]]; then
        break
    fi
    echo Deploy command failed, retrying in 5 secs...
    sleep 5
done
flask translate compile
exec gunicorn -b :5000 --access-logfile - --error-logfile - centd:app

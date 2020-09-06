FROM python:3.6-alpine

RUN adduser -D centd

WORKDIR /home/centd

COPY requirements.txt requirements.txt
RUN python -m venv venv
RUN venv/bin/pip install -r requirements.txt
RUN venv/bin/pip install gunicorn pymysql

COPY app app
COPY migrations migrations
COPY centd.py config.py boot.sh ./
RUN chmod a+x boot.sh

ENV FLASK_APP centd.py

RUN chown -R centd:centd ./
USER centd

EXPOSE 5000
ENTRYPOINT ["./boot.sh"]

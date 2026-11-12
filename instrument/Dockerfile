FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    INSTRUMENT_PORT=9000

WORKDIR /srv

COPY server.py ./

RUN useradd --create-home --uid 10002 --shell /usr/sbin/nologin instrument

USER instrument

EXPOSE 9000

CMD ["python3", "server.py"]

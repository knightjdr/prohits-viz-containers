FROM python:3

WORKDIR /app

COPY ./utilities/python/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY ./utilities/python/ .

RUN mkdir /files
WORKDIR /files

ENTRYPOINT ["python"]
FROM ubuntu:20.04
RUN apt update
RUN apt install -y curl cron
COPY bagent.sh /app/
RUN /app/bagent.sh


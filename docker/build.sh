#!/usr/bin/env bash
docker build -t wot .
docker tag wot regevlab/wot
docker push regevlab/wot

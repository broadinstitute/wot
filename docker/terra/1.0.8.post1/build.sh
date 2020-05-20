#!/usr/bin/env bash

docker build -t wot-terra-1.0.8.post1 .
docker tag wot-terra-1.0.8.post1 klarmancellobservatory/wot-terra:1.0.8.post1
docker push klarmancellobservatory/wot-terra:1.0.8.post1
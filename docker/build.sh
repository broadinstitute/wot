#!/usr/bin/env bash
export WOT_VERSION=1.0.2
docker build --build-arg WOT_VERSION=$WOT_VERSION -t wot-$WOT_VERSION .
docker tag wot regevlab/wot-$WOT_VERSION
docker push regevlab/wot-$WOT_VERSION

#!/bin/bash

for x in *; do mv "$x" "${x%???}"; done


#!/bin/bash

# modify the below to account for how many characters we need to remove. layer10000_XX is to become layer10000
for x in layer*; do mv "$x" "${x%??}"; done


#!/bin/bash

set -e
errors=0

# Run unit tests
python fastqe-convert/fastqe-convert_test.py || {
    echo "'python python/fastqe-convert/fastqe-convert_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E fastqe-convert/*.py || {
    echo 'pylint -E fastqe-convert/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"

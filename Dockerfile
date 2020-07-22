FROM python:3.7.3-stretch
WORKDIR /fastqe-convert
COPY . .

# Install the python package (and executable)
RUN pip3 install .

# Override some of the dependencies with the hard-coded versions
RUN pip3 install -r requirements-dev.txt

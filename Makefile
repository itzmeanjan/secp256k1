PYTHON = python3

all: testing

testing: ecdsa/*.py field/*.py point/*.py
	$(PYTHON) -m pytest -v --cache-clear

clean:
	find . -name __pycache__ -o -name .pytest* -o -name .benchmarks | xargs rm -rf

format:
	find . -name '*.py' | xargs $(PYTHON) -m black

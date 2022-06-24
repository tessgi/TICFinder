.PHONY: all clean pytest flake8 black isort

CMD:=poetry run
PYMODULE:=src
TESTS:=tests

# Run all the checks which do not change files
all: pytest flake8 isort black

# Run the unit tests using `pytest`
pytest:
	$(CMD) pytest $(PYMODULE) $(TESTS) --remote-data

# Lint the code using `flake8`
flake8:
	$(CMD) flake8 $(PYMODULE) $(TESTS)

# Automatically format the code using `black`
black:
	$(CMD) black $(PYMODULE) $(TESTS)

# Order the imports using `isort`
isort:
	$(CMD) isort $(PYMODULE) $(TESTS)

# Serve docs
serve:
	$(CMD) mkdocs serve

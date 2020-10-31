.PHONY: install release dist test
 
install:
	pip install -e .

release:
ifndef VERSION
	$(error "Usage: make release VERSION=0.1.1")
endif
	git checkout master
	git pull
	echo "__version__ = '$(VERSION)'" > abnumber/__version__.py
	git add abnumber/__version__.py
	git commit -m "Set version to $(VERSION)"
	git push
	@echo "Create a new release version on: https://github.com/prihoda/abnumber/releases"

dist:
	python setup.py sdist bdist_wheel

test: unit-test jupyter-test

jupyter-test: examples/AbNumber_getting_started.ipynb
	papermill $< > /dev/null

unit-test:
	pytest test

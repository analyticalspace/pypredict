SHELL=/bin/bash

all: venv
	@if [[ -e requirements.txt ]] ; then source venv/bin/activate && pip install -r requirements.txt; fi
	@source venv/bin/activate && ./setup.py build install
	@echo
	@echo 'To activate virtualenv:'
	@echo 'source venv/bin/activate'
	@echo

venv:
	@virtualenv -p python3 venv

clean:
	$(RM) -r venv

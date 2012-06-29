ifndef HELPERS_MK
HELPERS_MK=helpers.mk

clean:
	rm -f bootstrap-prefix.sh
	rm -rvf tmp/

help:
	@cat README.txt
	@echo
	@echo "Available actions (targets):"
	@echo "----------------------------"
	@./utils/list_make_targets.sh
	@echo
	@echo "To see which commands will be executed by each action, use:"
	@echo "\$$ make -n action"
	@echo

endif


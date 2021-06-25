SUBDIRS := $(subst /,,$(dir $(shell ls */Makefile)))

SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))

.PHONY: all checkfail copytodoc clean $(SUBDIRS) $(SUBDIRSCLEAN)

all: $(SUBDIRS)
	scripts/check

checkfail:
	scripts/check --fail

copytodoc:
	scripts/copytodoc

clean: $(SUBDIRSCLEAN)


$(SUBDIRS):
	$(MAKE) -C $@

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean


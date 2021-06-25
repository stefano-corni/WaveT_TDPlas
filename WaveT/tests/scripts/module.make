SUBDIRS := $(subst /Makefile,,$(wildcard */Makefile))
SUBDIRSCLEAN := $(addsuffix .clean,$(SUBDIRS))

.PHONY: all clean $(SUBDIRS) $(SUBDIRSCLEAN)

all: $(SUBDIRS)

clean: $(SUBDIRSCLEAN)

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBDIRSCLEAN): %.clean:
	$(MAKE) -C $* clean


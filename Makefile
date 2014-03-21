# This makefile builds the command-line version of Morpheus for a unix
# environment. Make sure Mono's xbuild program is on your path.

uc = $(if $(filter "$@","release"),"Release","Debug")

.PHONY: all
all: release debug

.PHONY: release debug
release debug:
	mkdir -p build
	rm -rf build/$@
	xbuild /p:Configuration=$(uc) "Morpheus/Morpheus (MonoDevelop).sln"
	cp -rp "Morpheus/Morpheus (command line)/bin/$(uc)" "build/$@"

.PHONY: clean
clean:
	rm -rf "external/CommandLine/CommandLine/bin/"
	rm -rf "external/CommandLine/CommandLine/obj/"
	rm -rf "Morpheus/bin/"
	rm -rf "Morpheus/obj/"
	rm -rf "Morpheus/Morpheus (command line)/bin/"
	rm -rf "Morpheus/Morpheus (command line)/obj/"
	rm -rf "build"


SOURCES = eopd.c\
          Makefile COPYRIGHT.txt LICENSE.txt README.md

all: build/eopd

clean:
	rm -rf build
	rm -rf dist

build/eopd: eopd.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

sources: dist/eopd-sources.zip dist/eopd-sources.tar.gz

dist/eopd-sources.zip: $(SOURCES)
	mkdir -p dist
	zip dist/eopd-sources $(SOURCES)

dist/eopd-sources.tar.gz: $(SOURCES)
	mkdir -p dist
	tar czf dist/eopd-sources.tar.gz $(SOURCES)

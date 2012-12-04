ifndef PREFIX
    PREFIX := $(PWD)
endif

ifndef CONFIGURE
    CONFIGURE := --with-sockets
endif

ifndef GA
    GA := ga-5-1
endif

# http://www.emsl.pnl.gov/docs/global/download/ga-5-1.tgz

ga: $(GA).tgz
	tar -xzf $(GA).tgz
	patch -p0 < $(GA).patch
	cd $(GA) && \
           ./configure --prefix=$(PREFIX) $(CONFIGURE) CFLAGS="$(CFLAGS) -g" FFLAGS="$(FFLAGS) -g" && \
           make clean && make install

clean:
	rm -fr $(GA)
	rm -fr bin/ include/ lib/

$(GA).tgz:
	wget http://www.emsl.pnl.gov/docs/global/download/$(GA).tgz


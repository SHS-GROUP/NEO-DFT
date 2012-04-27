check-syntax:
	g++ -Wall -fsyntax-only -S ${CHK_SOURCES} $(AM_CPPFLAGS)
